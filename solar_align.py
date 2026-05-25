"""
Solar image alignment tools — Python port of J. Chae's IDL routines.

Ports:
    ch_shift_sub.pro    -> ch_shift_sub()
    ch_alignoffset.pro  -> ch_alignoffset()
    align_hri.pro       -> align_fits_directory()

Differences vs the IDL originals (intentional):
    * Boundary handling: out-of-bounds pixels are filled with NaN (or a chosen
      cval) instead of being wrapped circularly. Output array shape is identical
      to input.
    * Reference policy: ALL frames are aligned to the FIRST frame (fixed
      reference), not a sliding reference. This avoids drift accumulation.
    * Multi-pass refinement re-applies the *accumulated* offset to the original
      data instead of chaining interpolations, so we incur only one
      interpolation per output frame.

Author: ported for Song, 2026.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple, Union

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.fft import fft2, ifft2
from scipy.ndimage import shift as ndi_shift

try:
    import sunpy.map  # noqa: F401
    HAS_SUNPY = True
except Exception:
    HAS_SUNPY = False


# ---------------------------------------------------------------------------
# Layer 1 — subpixel shift (non-cyclic)
# ---------------------------------------------------------------------------

def ch_shift_sub(image: np.ndarray,
                 x0: float,
                 y0: float,
                 order: int = 3,
                 cval: float = np.nan) -> np.ndarray:
    """
    Shift a 2-D image with subpixel accuracy.

    Out-of-bounds pixels are filled with ``cval`` (default NaN) — NOT wrapped.
    Output shape == input shape.

    Parameters
    ----------
    image : (ny, nx) ndarray
        Input image. NumPy convention: axis 0 = y (row), axis 1 = x (col).
    x0, y0 : float
        Shift in pixels. ``x0 > 0`` moves content toward larger x.
    order : int, default 3
        Spline interpolation order. 3 ≈ IDL ``cubic=-0.5``; use 1 for fastest.
    cval : float, default NaN
        Fill value for out-of-bounds.
    """
    if image.ndim != 2:
        raise ValueError(f"ch_shift_sub expects 2-D array, got shape {image.shape}")
    # scipy.ndimage.shift takes shift = (axis0, axis1) = (y, x)
    return ndi_shift(image, shift=(y0, x0), order=order,
                     mode='constant', cval=cval)


# ---------------------------------------------------------------------------
# Layer 2 — FFT cross-correlation alignment
# ---------------------------------------------------------------------------

def _gaussian_apodization(ny: int, nx: int) -> np.ndarray:
    """Gaussian apodization window with sigma = N/6 (matches IDL routine).

    The window is multiplied onto BOTH image and reference, so we take sqrt
    of the gaussian so that the effective window in the cross-correlation
    integrand is a single full-amplitude gaussian.
    """
    sigma_x = nx / 6.0
    sigma_y = ny / 6.0
    yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float64)
    xx -= nx / 2.0
    yy -= ny / 2.0
    win = np.exp(-0.5 * ((xx / sigma_x) ** 2 + (yy / sigma_y) ** 2))
    return np.sqrt(win)


def ch_alignoffset(image: np.ndarray,
                   reference: np.ndarray,
                   apodize: bool = True,
                   pad_to_pow2: bool = False,
                   compute_cor: bool = True
                   ) -> Tuple[Tuple[float, float], float]:
    """
    FFT-based cross-correlation alignment with subpixel parabolic refinement.

    Convention: the returned ``(dx, dy)`` is the offset of ``image`` relative
    to ``reference``. Calling ``ch_shift_sub(image, -dx, -dy)`` brings the
    image into alignment with the reference.

    Parameters
    ----------
    image, reference : 2-D ndarray, same shape
        Inputs may contain NaN (treated as zero after mean removal).
    apodize : bool, default True
        Apply Gaussian window to suppress FFT edge artifacts.
    pad_to_pow2 : bool, default False
        Resample to power-of-2 dims before FFT (legacy IDL behaviour).
        Modern FFT libs handle arbitrary sizes well, so default is off.
    compute_cor : bool, default True
        Compute the normalized cross-correlation at the recovered offset.

    Returns
    -------
    (dx, dy) : tuple of float
        Subpixel offset in *original* pixel units.
    cor : float
        Normalized cross-correlation coefficient (NaN if ``compute_cor=False``).
    """
    if image.shape != reference.shape:
        raise ValueError(f"Shape mismatch: {image.shape} vs {reference.shape}")

    ny_orig, nx_orig = image.shape

    img1 = np.asarray(image, dtype=np.float64)
    ref1 = np.asarray(reference, dtype=np.float64)
    img1 = img1 - np.nanmean(img1)
    ref1 = ref1 - np.nanmean(ref1)
    img1 = np.nan_to_num(img1, nan=0.0)
    ref1 = np.nan_to_num(ref1, nan=0.0)

    # Optional resample to power-of-2 dims
    if pad_to_pow2:
        from scipy.ndimage import zoom
        nx = 1 << (nx_orig - 1).bit_length()
        ny = 1 << (ny_orig - 1).bit_length()
        if nx < nx_orig:
            nx *= 2
        if ny < ny_orig:
            ny *= 2
        if (nx, ny) != (nx_orig, ny_orig):
            img1 = zoom(img1, (ny / ny_orig, nx / nx_orig), order=3)
            ref1 = zoom(ref1, (ny / ny_orig, nx / nx_orig), order=3)
    else:
        ny, nx = ny_orig, nx_orig

    # Apodize
    if apodize:
        win = _gaussian_apodization(ny, nx)
        img1 = img1 * win
        ref1 = ref1 * win

    # Cross-correlation via FFT theorem
    F_img = fft2(img1)
    F_ref = fft2(ref1)
    cor_map = np.real(ifft2(F_img * np.conj(F_ref)))

    # Integer-pixel peak — numpy index order is (row, col) = (y, x)
    y0_int, x0_int = np.unravel_index(np.argmax(cor_map), cor_map.shape)
    if x0_int > nx // 2:
        x0_int -= nx
    if y0_int > ny // 2:
        y0_int -= ny

    # Subpixel parabolic refinement on the 3x3 neighborhood of the peak
    cor_roll = np.roll(cor_map, shift=(-y0_int + 1, -x0_int + 1), axis=(0, 1))
    cc = cor_roll[:3, :3]  # cc[row(y), col(x)]
    denom_x = cc[1, 2] + cc[1, 0] - 2.0 * cc[1, 1]
    denom_y = cc[2, 1] + cc[0, 1] - 2.0 * cc[1, 1]
    x1 = 0.5 * (cc[1, 0] - cc[1, 2]) / denom_x if abs(denom_x) > 1e-30 else 0.0
    y1 = 0.5 * (cc[0, 1] - cc[2, 1]) / denom_y if abs(denom_y) > 1e-30 else 0.0

    x_off = (x0_int + x1) * nx_orig / nx
    y_off = (y0_int + y1) * ny_orig / ny

    cor_val = np.nan
    if compute_cor:
        img_aligned = ch_shift_sub(img1, -(x0_int + x1), -(y0_int + y1),
                                   order=3, cval=0.0)
        # Valid mask: pixels where both arrays have non-trivial signal
        a = np.nan_to_num(img_aligned, nan=0.0)
        b = ref1
        denom = np.sqrt(np.sum(a * a) * np.sum(b * b))
        cor_val = float(np.sum(a * b) / denom) if denom > 0 else 0.0

    return (float(x_off), float(y_off)), float(cor_val)


# ---------------------------------------------------------------------------
# Layer 3 — batch driver
# ---------------------------------------------------------------------------

def _read_fits(path: Path, use_sunpy: bool = True):
    """Return ``(data_array, metadata_or_map)``."""
    if use_sunpy and HAS_SUNPY:
        m = sunpy.map.Map(str(path))
        return np.asarray(m.data), m
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                return np.asarray(hdu.data), hdu.header.copy()
    raise IOError(f"No image HDU found in {path}")


def _read_fits_shape(path: Path):
    """Header-only shape read (fast, doesn't load pixel data)."""
    with fits.open(path) as hdul:
        for hdu in hdul:
            naxis = hdu.header.get("NAXIS", 0)
            if naxis >= 2:
                n1 = hdu.header.get("NAXIS1")
                n2 = hdu.header.get("NAXIS2")
                if n1 and n2:
                    return (int(n2), int(n1))  # numpy (ny, nx)
    return None


def _center_crop(data: np.ndarray, target_shape) -> np.ndarray:
    """Center-crop a 2-D array to target_shape (ny, nx). No-op if already that size."""
    if data.shape == target_shape:
        return data
    cy, cx = data.shape
    ty, tx = target_shape
    if ty > cy or tx > cx:
        raise ValueError(f"Cannot crop {data.shape} -> {target_shape} (target larger)")
    y0 = (cy - ty) // 2
    x0 = (cx - tx) // 2
    return data[y0:y0 + ty, x0:x0 + tx]


def scan_fits_shapes(input_dir: Union[str, Path],
                     pattern: str = "*.fits",
                     verbose: bool = True) -> dict:
    """
    Scan a directory and report the (ny, nx) shape of every FITS file.

    Useful before calling ``align_fits_directory`` when you suspect a mix of
    image sizes. Reads only headers (fast).

    Returns
    -------
    dict mapping shape tuple -> list of filenames.
    """
    from collections import defaultdict
    input_dir = Path(input_dir)
    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {input_dir}")
    groups = defaultdict(list)
    for f in files:
        s = _read_fits_shape(f)
        if s is not None:
            groups[s].append(f.name)
    if verbose:
        print(f"Scanned {len(files)} files in {input_dir}")
        print(f"Found {len(groups)} distinct shape(s):\n")
        for sh, names in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            print(f"  {sh}: {len(names)} files")
            for n in names[:3]:
                print(f"    {n}")
            if len(names) > 3:
                print(f"    ... and {len(names)-3} more (last: {names[-1]})")
            print()
    return dict(groups)


def _write_fits(path: Path, data: np.ndarray, meta,
                use_sunpy: bool = True, overwrite: bool = True) -> None:
    """Write FITS preserving as much metadata as available."""
    is_map = (HAS_SUNPY and hasattr(meta, "meta")
              and not isinstance(meta, (dict, fits.Header)))
    if use_sunpy and is_map:
        new_map = sunpy.map.Map(data, meta.meta)
        new_map.save(str(path), overwrite=overwrite)
        return
    if isinstance(meta, fits.Header):
        header = meta
    elif isinstance(meta, dict):
        header = fits.Header(meta)
    else:
        header = fits.Header(dict(meta.meta) if hasattr(meta, "meta") else {})
    fits.writeto(str(path), data, header, overwrite=overwrite)


def _get_time_str(meta) -> str:
    if hasattr(meta, "date") and not isinstance(meta, (dict, fits.Header)):
        try:
            return str(meta.date)
        except Exception:
            pass
    if isinstance(meta, (dict, fits.Header)):
        for key in ("DATE-OBS", "DATE_OBS", "DATEOBS", "T_OBS"):
            if key in meta:
                return str(meta[key])
    return ""


def align_fits_directory(input_dir: Union[str, Path],
                         output_dir: Union[str, Path],
                         pattern: str = "*.fits",
                         *,
                         use_sunpy: bool = True,
                         cval: float = np.nan,
                         cor_thresholds: Tuple[float, ...] = (0.85, 0.95),
                         interp_order: int = 3,
                         apodize: bool = True,
                         auto_crop: bool = True,
                         max_shape_diff: float = 0.05,
                         overwrite: bool = True,
                         verbose: bool = True
                         ) -> pd.DataFrame:
    """
    Align all FITS files in ``input_dir`` to the FIRST file; write results to
    ``output_dir``; also save an offset summary (CSV + NPZ).

    Parameters
    ----------
    input_dir, output_dir : str or Path
    auto_crop : bool, default True
        If files in the directory have slightly different shapes (e.g. 864x864
        vs 865x865 from different preprocessing batches), center-crop all
        frames down to the common minimum shape. If False, mismatched frames
        are skipped (legacy behavior).
    max_shape_diff : float, default 0.05
        Maximum relative shape difference allowed when auto-cropping (5%).
        Bigger discrepancies abort with an error, since they probably mean
        the files belong to different datasets and should not be aligned.
    pattern : str
        Glob pattern (e.g. ``'*.fits'``, ``'solo_L2_eui-hri*.fits'``).
    use_sunpy : bool
        Use sunpy.map.Map for I/O when available (preserves WCS).
    cval : float
        Fill value for out-of-bounds after the shift (NaN by default).
    cor_thresholds : tuple of float
        Re-align if the running cor drops below each successive threshold.
        IDL default: (0.85, 0.95) — up to 3 passes total.
    interp_order : int
        Spline order for the final shift (3 = cubic).
    apodize : bool
        Apply Gaussian apodization inside the FFT alignment.
    overwrite : bool
        Overwrite existing output files.
    verbose : bool
        Print per-frame progress.

    Returns
    -------
    df : pandas.DataFrame with columns
        filename, time, x_offset, y_offset, correlation, n_iter
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern} in {input_dir}")

    if use_sunpy and not HAS_SUNPY:
        if verbose:
            print("sunpy unavailable; falling back to astropy.io.fits")
        use_sunpy = False

    if verbose:
        print(f"Found {len(files)} files; reference = {files[0].name}")

    # --- Shape scan (header-only, fast) and determine common target shape ---
    target_shape = None
    if auto_crop:
        from collections import defaultdict
        shape_groups = defaultdict(list)
        for f in files:
            s = _read_fits_shape(f)
            if s is not None:
                shape_groups[s].append(f.name)
        if not shape_groups:
            raise IOError("No 2-D image HDUs found in any file")
        shapes = list(shape_groups.keys())
        ny_min = min(s[0] for s in shapes)
        nx_min = min(s[1] for s in shapes)
        ny_max = max(s[0] for s in shapes)
        nx_max = max(s[1] for s in shapes)
        target_shape = (ny_min, nx_min)
        # Sanity check: if discrepancy is huge, this isn't a preprocessing
        # quirk — files probably belong to different datasets.
        rel_diff = max((ny_max - ny_min) / ny_min,
                       (nx_max - nx_min) / nx_min)
        if rel_diff > max_shape_diff:
            # Build a detailed breakdown so user can see which files are which.
            breakdown_lines = ["Shape distribution:"]
            for sh, names in sorted(shape_groups.items(),
                                    key=lambda kv: -len(kv[1])):
                breakdown_lines.append(
                    f"  {sh}: {len(names)} files  "
                    f"(e.g. {names[0]}{', ...' if len(names) > 1 else ''})"
                )
            raise ValueError(
                f"Shape discrepancy too large to auto-crop: min={target_shape}, "
                f"max=({ny_max}, {nx_max}), relative diff {rel_diff:.1%} > "
                f"max_shape_diff={max_shape_diff:.1%}.\n"
                + "\n".join(breakdown_lines) +
                "\nOptions:\n"
                "  1. Filter the `pattern` argument to one consistent group\n"
                "  2. Move outlier files out of the input directory\n"
                "  3. Pass max_shape_diff=... to force auto-crop anyway\n"
                "  4. Call scan_fits_shapes(input_dir) to inspect first"
            )
        n_to_crop = sum(len(names) for sh, names in shape_groups.items()
                        if sh != target_shape)
        if verbose and n_to_crop:
            print(f"Shape scan: target={target_shape}, will center-crop "
                  f"{n_to_crop}/{len(files)} files (diff up to "
                  f"{ny_max-ny_min},{nx_max-nx_min} px)")

    # --- Frame 0: reference ---
    ref_data, ref_meta = _read_fits(files[0], use_sunpy=use_sunpy)
    if ref_data.ndim != 2:
        raise ValueError(f"Reference must be 2-D, got shape {ref_data.shape}")
    if auto_crop:
        ref_data = _center_crop(ref_data, target_shape)
    ref_for_corr = np.asarray(ref_data, dtype=np.float64)

    _write_fits(output_dir / files[0].name, ref_data, ref_meta,
                use_sunpy=use_sunpy, overwrite=overwrite)

    records = [{
        "filename": files[0].name,
        "time": _get_time_str(ref_meta),
        "x_offset": 0.0,
        "y_offset": 0.0,
        "correlation": 1.0,
        "n_iter": 0,
    }]

    # --- Subsequent frames ---
    for i, fpath in enumerate(files[1:], start=1):
        data, meta = _read_fits(fpath, use_sunpy=use_sunpy)

        if auto_crop and data.shape != ref_data.shape:
            try:
                data = _center_crop(data, ref_data.shape)
            except ValueError as e:
                if verbose:
                    print(f"  WARN: {fpath.name}: {e}, skipping")
                continue
        elif data.shape != ref_data.shape:
            if verbose:
                print(f"  WARN: {fpath.name} shape {data.shape} != "
                      f"reference {ref_data.shape}, skipping")
            continue

        data64 = np.asarray(data, dtype=np.float64)

        # Pass 1 always runs
        (dx, dy), cor = ch_alignoffset(data64, ref_for_corr, apodize=apodize)
        total_dx, total_dy = dx, dy
        n_iter = 1

        # Additional passes: each gated by its threshold
        for thresh in cor_thresholds:
            if cor >= thresh:
                continue
            # Re-shift the ORIGINAL data with the cumulative offset, then
            # measure the residual. This keeps only ONE interpolation chain.
            current = ch_shift_sub(data64, -total_dx, -total_dy,
                                   order=interp_order, cval=cval)
            (dx, dy), cor = ch_alignoffset(current, ref_for_corr,
                                           apodize=apodize)
            total_dx += dx
            total_dy += dy
            n_iter += 1

        # Apply the final accumulated shift once to the original data
        aligned = ch_shift_sub(data64, -total_dx, -total_dy,
                               order=interp_order, cval=cval)
        # Preserve original dtype if it was integer/float32 etc.
        aligned_out = aligned.astype(data.dtype, copy=False)

        _write_fits(output_dir / fpath.name, aligned_out, meta,
                    use_sunpy=use_sunpy, overwrite=overwrite)

        records.append({
            "filename": fpath.name,
            "time": _get_time_str(meta),
            "x_offset": total_dx,
            "y_offset": total_dy,
            "correlation": cor,
            "n_iter": n_iter,
        })

        if verbose:
            print(f"  [{i:4d}/{len(files)-1}] {fpath.name}: "
                  f"dx={total_dx:+8.4f}  dy={total_dy:+8.4f}  "
                  f"cor={cor:.4f}  iter={n_iter}")

    # --- Offset summary ---
    df = pd.DataFrame(records)
    df.to_csv(output_dir / "offsets_summary.csv", index=False)
    # Use fixed-width unicode for strings (avoids np.load(allow_pickle=True))
    np.savez(output_dir / "offsets_summary.npz",
             filenames=np.array(df["filename"].tolist(), dtype="U256"),
             times=np.array(df["time"].astype(str).tolist(), dtype="U64"),
             x_offset=df["x_offset"].to_numpy(dtype=np.float64),
             y_offset=df["y_offset"].to_numpy(dtype=np.float64),
             correlation=df["correlation"].to_numpy(dtype=np.float64),
             n_iter=df["n_iter"].to_numpy(dtype=np.int32))

    if verbose:
        print(f"\nWrote {len(df)} aligned FITS to {output_dir}")
        print(f"Summary: {output_dir/'offsets_summary.csv'} (+ .npz)")

    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Align solar FITS images using FFT cross-correlation.")
    p.add_argument("input_dir", help="Directory of input FITS files")
    p.add_argument("output_dir", help="Directory for aligned FITS files")
    p.add_argument("--pattern", default="*.fits", help="Glob pattern")
    p.add_argument("--no-sunpy", action="store_true",
                   help="Disable sunpy I/O (use astropy.io.fits)")
    p.add_argument("--cval", type=float, default=float('nan'),
                   help="Fill value for out-of-bounds (default: NaN)")
    p.add_argument("--no-apodize", action="store_true",
                   help="Disable Gaussian apodization")
    args = p.parse_args()

    align_fits_directory(
        args.input_dir, args.output_dir,
        pattern=args.pattern,
        use_sunpy=not args.no_sunpy,
        cval=args.cval,
        apodize=not args.no_apodize,
    )


if __name__ == "__main__":
    main()
