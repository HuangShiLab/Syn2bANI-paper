#!/usr/bin/env python3
"""Compose a grid of PNG panels into one figure.

Each input is a path. Layout is specified as rows separated by '--row'.
Example: python compose_grid.py a.png b.png --row c.png d.png -o out.png
"""
import argparse
from pathlib import Path
from PIL import Image


def parse_layout(args):
    """Return list of rows, each row is list of paths."""
    rows = [[]]
    for item in args:
        if item == 'ROW':
            rows.append([])
        else:
            rows[-1].append(Path(item))
    return [r for r in rows if r]


def resize_to_height(img, h):
    if img.height == h:
        return img
    ratio = h / img.height
    return img.resize((int(img.width * ratio), h), Image.LANCZOS)


def compose(rows, out, pad=20, bg=(255, 255, 255)):
    imgs = [[Image.open(p).convert('RGB') for p in row] for row in rows]
    # For each row, normalize heights to the smallest height in that row
    norm = []
    for row in imgs:
        h = min(i.height for i in row)
        norm.append([resize_to_height(i, h) for i in row])
    row_heights = [max(i.height for i in row) for row in norm]
    row_widths = [sum(i.width for i in row) + pad * (len(row) - 1) for row in norm]
    total_w = max(row_widths)
    total_h = sum(row_heights) + pad * (len(norm) - 1)
    canvas = Image.new('RGB', (total_w, total_h), bg)
    y = 0
    for row, h in zip(norm, row_heights):
        x = 0
        for img in row:
            canvas.paste(img, (x, y))
            x += img.width + pad
        y += h + pad
    canvas.save(out, dpi=(300, 300))
    print(f'Wrote {out} ({total_w}x{total_h})')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('items', nargs='+', help='image paths and --row separators')
    parser.add_argument('-o', '--out', required=True, type=Path)
    args = parser.parse_args()
    rows = parse_layout(args.items)
    compose(rows, args.out)


if __name__ == '__main__':
    main()
