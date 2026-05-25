# solar_align

**[English]** | [中文](README_zh.md)

> A Python toolkit for de-jittering solar image sequences.
> A modern port of J. Chae's IDL routines `CH_ALIGNOFFSET` and `CH_SHIFT_SUB`,
> with extensions for non-cyclic shifts, multi-shape inputs, and offset summaries.
> Validated on ground-based NVST/Hα and space-based Solar Orbiter EUI/HRI data.

![python](https://img.shields.io/badge/python-%E2%89%A53.8-blue)
![license](https://img.shields.io/badge/license-MIT-green)
![status](https://img.shields.io/badge/status-tested%20on%20NVST%20%2B%20EUI/HRI-success)

---

## What it does

Given a directory of FITS image frames affected by frame-to-frame jitter
(atmospheric seeing residuals, telescope mount vibration, spacecraft pointing
instabilities, etc.), `solar_align` aligns every frame to a common reference
using FFT-based cross-correlation with subpixel refinement, and writes:

- aligned FITS files (same filenames, same shapes, no cyclic wrap-around)
- a per-frame offset summary (`offsets_summary.csv` and `.npz`)

The algorithmic core (FFT cross-correlation + 3×3 parabolic subpixel fit +
Gaussian apodization + multi-pass refinement) is unchanged from
Chae's IDL implementation, which has been used in the solar physics community
since ~1999.

## Features

- **Subpixel accuracy** (≈0.05–0.15 px single-pass, ≈0.03 px after two iterations)
- **Non-cyclic shifts** — out-of-bounds pixels filled with `NaN` (or user-chosen value), not wrapped
- **Single-interpolation guarantee** — multi-iteration refinement re-applies the cumulative offset to the *original* data, avoiding cumulative blurring
- **Auto-handling of shape-inconsistent inputs** — center-crops to a common minimum shape with safety thresholds
- **Sunpy-native I/O** preserves WCS and instrument metadata
- **Pure Python**, no compiled dependencies

## Installation

```bash
pip install scipy astropy pandas sunpy
```

Then clone or download this repository, or simply drop `solar_align.py` next
to your script / notebook.

## Quick start

```python
from solar_align import align_fits_directory

df = align_fits_directory(
    input_dir  = '/path/to/raw_fits',
    output_dir = '/path/to/aligned_fits',
    pattern    = '*.fits',
)
```

Produces:

```
/path/to/aligned_fits/
    ├── frame_001.fits          ← aligned, same filename, same shape
    ├── frame_002.fits
    ├── ...
    ├── offsets_summary.csv     ← per-frame offsets, human-readable
    └── offsets_summary.npz     ← per-frame offsets, machine-readable
```

The first file in `input_dir` (alphabetically) is used as the fixed reference.
Every subsequent frame is aligned to it.

### Inspecting the result

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/path/to/aligned_fits/offsets_summary.csv')
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(df['x_offset'], '.-', label='dx')
ax.plot(df['y_offset'], '.-', label='dy')
ax.set_xlabel('Frame'); ax.set_ylabel('Offset (px)')
ax.legend()
```

The CSV contains:

| column | meaning |
|---|---|
| `filename` | output filename |
| `time` | observation time from header |
| `x_offset`, `y_offset` | measured shift in pixels |
| `correlation` | normalized cross-correlation after alignment (1.0 = perfect) |
| `n_iter` | number of refinement passes applied |

## Repository layout

```
solar_align/
├── solar_align.py          ← core module
├── fits_to_png.py          ← render aligned FITS to PNGs for inspection
├── make_compare_gif.py     ← build raw-vs-aligned comparison GIFs
├── example_usage.py        ← end-to-end example script
├── test_solar_align.py     ← unit + integration tests on synthetic data
├── README.md               ← this file
└── README_zh.md            ← Chinese version
```

## Tool API summary

| Function | Purpose |
|---|---|
| `align_fits_directory(input_dir, output_dir, ...)` | batch alignment of a directory |
| `ch_alignoffset(image, reference)` | measure offset between two arrays |
| `ch_shift_sub(image, x0, y0)` | shift an array with subpixel accuracy (non-cyclic) |
| `scan_fits_shapes(input_dir)` | report the (ny, nx) distribution of files in a directory |

Full signatures and parameters are documented in the source.

## Tested on

| Dataset | Wavelength | Telescope | Notes |
|---|---|---|---|
| NVST | Hα 6562.8 Å | ground-based | sustained seeing-induced jitter |
| Solar Orbiter EUI/HRI | EUV 174 Å | space-based | high-cadence (1–3 s), platform jitter |
| Synthetic | — | — | offset recovery accuracy quantified |

Synthetic-data recovery accuracy: single-pass `≤ 0.15 px`, two-pass `≤ 0.05 px`
on 30-frame sequences with ~1 px rms input jitter.

## Failure modes (brief)

The fixed-reference + FFT-correlation strategy has known limits:

- Long sequences where the Sun itself evolves (correlation drifts down) → use segmented alignment
- Very large jitter (>50 px) → parabolic fit becomes unreliable, do a coarse pre-alignment first
- Low SNR (cor < 0.7) → spatial smoothing or temporal averaging before alignment
- Bright transient events dominating the field → mask the dynamic region (not yet implemented; planned)
- Poor reference frame (bad seeing, cloud) → manually choose a good reference

A detailed discussion of failure modes and recommended mitigations is in the
accompanying article (Chinese, see WeChat post linked in `README_zh.md`).

## Reusing the offset CSV

The per-frame offsets are themselves physical data. Common secondary analyses:

- Power spectrum of `x_offset` / `y_offset` → identify periodic disturbance sources
- Long-term drift trend → thermal vs mechanical
- Cross-correlation between two instruments' jitter series → common disturbance events
- `correlation` column → quality weight for time-series analysis

## Credits

The algorithmic core (FFT cross-correlation, 3×3 parabolic subpixel fit, Gaussian
apodization, multi-pass refinement) is unchanged from
**J. Chae's** IDL routines `CH_ALIGNOFFSET.pro` and `CH_SHIFT_SUB.pro`
(Seoul National University), originally written around 1999 and incorporated
into the FISS data pipeline. This Python port is a re-implementation and
modernization; all credit for the underlying method goes to Chae.

If you use this tool in published work, please acknowledge both the Python
port and the original IDL routines.

## License

MIT License. See `LICENSE` for details.

## Contact

Issues, feature requests, and pull requests are welcome via GitHub.

Author: Song Tan, Leibniz Institute for Astrophysics Potsdam (AIP)
Contact: [your email here]
