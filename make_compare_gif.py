"""
Build a side-by-side comparison GIF from two PNG directories.

Reads paired PNGs (by sorted filename) from PNG_RAW_DIR and PNG_ALIGNED_DIR,
stitches each pair horizontally with a labeled banner, saves an animated GIF.
"""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont


# ============================ CONFIG ============================
PNG_RAW_DIR     = Path("/Users/tansong/Downloads/hri_pngs_raw")
PNG_ALIGNED_DIR = Path("/Users/tansong/Downloads/hri_pngs_aligned")
OUTPUT_GIF      = Path("/Users/tansong/Downloads/compare.gif")

DURATION_MS = 100      # ms per frame (100 = 10 fps)
LOOP        = 0        # 0 = loop forever
SCALE       = 0.6      # 0.6 = shrink to 60% to keep GIF size manageable
BANNER_H    = 28       # px height of top label banner
LABEL_LEFT  = "RAW"
LABEL_RIGHT = "ALIGNED"
GAP_PX      = 4        # black gap between the two panels
OPTIMIZE    = True     # palette optimize (smaller files, slower save)
# ================================================================


def _load_font(size=18):
    """Find a usable bold font on the system; fall back to PIL default."""
    candidates = [
        "/System/Library/Fonts/Helvetica.ttc",                # macOS
        "/System/Library/Fonts/HelveticaNeue.ttc",            # macOS
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",  # Linux
        "C:\\Windows\\Fonts\\arialbd.ttf",                    # Windows
        "C:\\Windows\\Fonts\\arial.ttf",                      # Windows
    ]
    for p in candidates:
        if Path(p).exists():
            try:
                return ImageFont.truetype(p, size)
            except Exception:
                continue
    return ImageFont.load_default()


def _combine_pair(left_path, right_path, font):
    """Open one PNG from each side, scale, add banner, return composite RGB image."""
    left  = Image.open(left_path).convert("RGB")
    right = Image.open(right_path).convert("RGB")

    if SCALE != 1.0:
        new_w = int(left.width * SCALE)
        new_h = int(left.height * SCALE)
        left  = left.resize((new_w, new_h),  Image.LANCZOS)
        right = right.resize((new_w, new_h), Image.LANCZOS)

    # Force both panels to same height (in case of off-by-one sizes)
    h = min(left.height, right.height)
    if left.height != h:
        left  = left.crop((0, 0, left.width, h))
    if right.height != h:
        right = right.crop((0, 0, right.width, h))

    total_w = left.width + GAP_PX + right.width
    total_h = BANNER_H + h
    canvas = Image.new("RGB", (total_w, total_h), color="black")

    # Paste images below the banner
    canvas.paste(left,  (0,                     BANNER_H))
    canvas.paste(right, (left.width + GAP_PX,   BANNER_H))

    # Draw banner labels
    draw = ImageDraw.Draw(canvas)
    # Left label
    draw.text((10, 4), LABEL_LEFT, fill="white", font=font)
    # Right label
    draw.text((left.width + GAP_PX + 10, 4), LABEL_RIGHT, fill="white", font=font)
    # Frame name (right-aligned, useful for debugging)
    name = left_path.stem
    try:
        bbox = draw.textbbox((0, 0), name, font=font)
        tw = bbox[2] - bbox[0]
    except AttributeError:  # older PIL
        tw = draw.textlength(name, font=font) if hasattr(draw, 'textlength') else len(name) * 8
    draw.text((total_w - tw - 10, 4), name, fill="gray", font=font)

    return canvas


def main():
    raw_files     = sorted(PNG_RAW_DIR.glob("*.png"))
    aligned_files = sorted(PNG_ALIGNED_DIR.glob("*.png"))

    if not raw_files:
        raise FileNotFoundError(f"No PNGs in {PNG_RAW_DIR}")
    if not aligned_files:
        raise FileNotFoundError(f"No PNGs in {PNG_ALIGNED_DIR}")

    # Match by filename (stem) so order can't drift
    raw_map     = {p.name: p for p in raw_files}
    aligned_map = {p.name: p for p in aligned_files}
    common = sorted(set(raw_map) & set(aligned_map))

    n_skipped = len(raw_files) + len(aligned_files) - 2 * len(common)
    print(f"raw: {len(raw_files)}, aligned: {len(aligned_files)}, "
          f"common: {len(common)}  (skipped {n_skipped} unmatched)")

    font = _load_font(size=18)

    frames = []
    for i, name in enumerate(common):
        frame = _combine_pair(raw_map[name], aligned_map[name], font)
        frames.append(frame)
        if (i + 1) % 20 == 0 or i == len(common) - 1:
            print(f"  composed [{i+1:4d}/{len(common)}]")

    OUTPUT_GIF.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nSaving {len(frames)} frames -> {OUTPUT_GIF}")
    frames[0].save(
        OUTPUT_GIF,
        save_all=True,
        append_images=frames[1:],
        duration=DURATION_MS,
        loop=LOOP,
        optimize=OPTIMIZE,
        disposal=2,           # restore-to-background between frames
    )
    size_mb = OUTPUT_GIF.stat().st_size / (1024 * 1024)
    w, h = frames[0].size
    print(f"GIF size: {w}x{h}, {size_mb:.2f} MB, {DURATION_MS} ms/frame "
          f"({1000/DURATION_MS:.1f} fps)")


if __name__ == "__main__":
    main()
