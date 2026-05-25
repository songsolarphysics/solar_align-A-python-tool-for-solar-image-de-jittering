"""
Render PNGs from raw and aligned HRI FITS frames for side-by-side animation.

- Reads all FITS from INPUT_DIR  -> writes PNGs to PNG_RAW_DIR
- Reads all FITS from OUTPUT_DIR -> writes PNGs to PNG_ALIGNED_DIR
- Uses ONE shared color scale + log stretch so both animations are comparable
  and don't flicker frame-to-frame.
- NaN-safe (aligned frames have NaN edges from non-cyclic shift).
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')  # headless rendering
import matplotlib.pyplot as plt
from astropy.io import fits


# ============================ CONFIG ============================
INPUT_DIR  = Path("/Users/tansong/Downloads/hri_data")     # raw frames
OUTPUT_DIR = Path("/Users/tansong/Downloads/hri_data_de")  # aligned frames

# Where to put the PNGs (parallel folders next to the FITS dirs)
PNG_RAW_DIR     = INPUT_DIR.parent  / "hri_pngs_raw"
PNG_ALIGNED_DIR = OUTPUT_DIR.parent / "hri_pngs_aligned"

PATTERN     = "*.fits"
USE_LOG     = True               # log10 stretch (good for EUV)
PERCENTILE  = (1.0, 99.5)        # robust intensity clipping
SAMPLE_N    = 8                  # frames sampled to determine global scale
DPI         = 100
FIGSIZE     = (7, 7)             # inches, MUST be constant across frames
SHOW_TITLE  = True               # overlay filename + DATE-OBS
CMAP_NAME   = 'sdoaia171'        # AIA 171 LUT; fallback to 'inferno' if missing
# ================================================================


def _setup_cmap():
    """Try sunpy's AIA colormap; fall back if unavailable. NaN -> black."""
    cmap_name = CMAP_NAME
    try:
        import sunpy.visualization.colormaps  # noqa: F401 (registers 'sdoaia171')
    except Exception:
        cmap_name = 'inferno'
        print(f"  (sunpy colormaps unavailable, using '{cmap_name}')")
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color='black')   # NaN regions render as black
    return cmap


def _read_fits_data(path):
    """Return (data, date_obs_string) for the first 2-D image HDU."""
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = np.asarray(hdu.data, dtype=np.float64)
                date_obs = hdu.header.get('DATE-OBS',
                            hdu.header.get('DATE_OBS', ''))
                return data, date_obs
    raise IOError(f"No image HDU in {path}")


def _stretch(data):
    """Apply log10 stretch if requested; NaN-safe."""
    if not USE_LOG:
        return data
    out = np.full_like(data, np.nan)
    pos = data > 0
    out[pos] = np.log10(data[pos])
    return out


def compute_global_scale(fits_dir):
    """Robust intensity range computed from SAMPLE_N evenly-spaced frames."""
    files = sorted(fits_dir.glob(PATTERN))
    if not files:
        raise FileNotFoundError(f"No {PATTERN} in {fits_dir}")
    idxs = np.linspace(0, len(files) - 1,
                       min(SAMPLE_N, len(files))).astype(int)
    lows, highs = [], []
    for i in idxs:
        data, _ = _read_fits_data(files[i])
        s = _stretch(data)
        lo, hi = np.nanpercentile(s, PERCENTILE)
        lows.append(lo)
        highs.append(hi)
    return float(np.nanmin(lows)), float(np.nanmax(highs))


def render_directory(fits_dir, png_dir, vmin, vmax, cmap, label=""):
    """Render every FITS in fits_dir -> PNG of fixed size/scale in png_dir."""
    files = sorted(fits_dir.glob(PATTERN))
    png_dir.mkdir(exist_ok=True, parents=True)
    print(f"\n{label}  {fits_dir}  ->  {png_dir}  ({len(files)} files)")

    for i, f in enumerate(files):
        data, date_obs = _read_fits_data(f)
        img = _stretch(data)

        fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
        ax.imshow(img, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
                  interpolation='nearest')
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

        if SHOW_TITLE:
            ax.text(0.02, 0.98,
                    f"{f.name}\n{date_obs}",
                    transform=ax.transAxes, ha='left', va='top',
                    color='white', fontsize=9,
                    bbox=dict(facecolor='black', alpha=0.4, pad=2,
                              edgecolor='none'))

        out = png_dir / f"{f.stem}.png"
        fig.savefig(out, dpi=DPI, bbox_inches='tight', pad_inches=0,
                    facecolor='black')
        plt.close(fig)

        if (i + 1) % 20 == 0 or i == len(files) - 1:
            print(f"    [{i+1:4d}/{len(files)}]  {f.name}")

    return len(files)


def main():
    cmap = _setup_cmap()

    # Use the RAW directory to set the color scale, so the aligned frames
    # are displayed with the exact same stretch.
    vmin, vmax = compute_global_scale(INPUT_DIR)
    print(f"\nGlobal display scale  (log={USE_LOG}):")
    print(f"  vmin = {vmin:.4f}")
    print(f"  vmax = {vmax:.4f}")

    n_raw = render_directory(INPUT_DIR,  PNG_RAW_DIR,     vmin, vmax, cmap,
                             label="RAW    ")
    n_ali = render_directory(OUTPUT_DIR, PNG_ALIGNED_DIR, vmin, vmax, cmap,
                             label="ALIGNED")

    print(f"\nDone:  {n_raw} raw + {n_ali} aligned PNGs")
    print(f"\nTo build comparison movies with ffmpeg:")
    print(f"  ffmpeg -y -framerate 12 -pattern_type glob "
          f"-i '{PNG_RAW_DIR}/*.png' "
          f"-c:v libx264 -pix_fmt yuv420p -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' "
          f"raw.mp4")
    print(f"  ffmpeg -y -framerate 12 -pattern_type glob "
          f"-i '{PNG_ALIGNED_DIR}/*.png' "
          f"-c:v libx264 -pix_fmt yuv420p -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' "
          f"aligned.mp4")
    print(f"\nOr side-by-side single video:")
    print(f"  ffmpeg -y -framerate 12 -pattern_type glob "
          f"-i '{PNG_RAW_DIR}/*.png' "
          f"-framerate 12 -pattern_type glob "
          f"-i '{PNG_ALIGNED_DIR}/*.png' "
          f"-filter_complex hstack -c:v libx264 -pix_fmt yuv420p compare.mp4")


if __name__ == "__main__":
    main()
