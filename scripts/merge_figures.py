#!/usr/bin/env python3
"""Merge two PNG figures side-by-side or top-bottom into one composite."""
import argparse
from pathlib import Path
from PIL import Image


def merge(paths, out, direction='horizontal', pad=20, bg=(255, 255, 255)):
    imgs = [Image.open(p).convert('RGB') for p in paths]
    widths = [i.width for i in imgs]
    heights = [i.height for i in imgs]
    if direction == 'horizontal':
        total_w = sum(widths) + pad * (len(imgs) - 1)
        total_h = max(heights)
        new = Image.new('RGB', (total_w, total_h), bg)
        x = 0
        for i in imgs:
            y = (total_h - i.height) // 2
            new.paste(i, (x, y))
            x += i.width + pad
    else:
        total_w = max(widths)
        total_h = sum(heights) + pad * (len(imgs) - 1)
        new = Image.new('RGB', (total_w, total_h), bg)
        y = 0
        for i in imgs:
            x = (total_w - i.width) // 2
            new.paste(i, (x, y))
            y += i.height + pad
    new.save(out, dpi=(300, 300))
    print(f'Wrote {out}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', nargs='+', type=Path)
    parser.add_argument('-o', '--out', required=True, type=Path)
    parser.add_argument('-d', '--direction', choices=['horizontal', 'vertical'], default='horizontal')
    args = parser.parse_args()
    merge(args.inputs, args.out, args.direction)


if __name__ == '__main__':
    main()
