#!/usr/bin/env python3
"""Generate the packaged Extended Data Figure 2.

Input: ``paper/figures/main/fig6_syntracker_structure.png``.
Output: ``paper/figures/extended_data/extended_data_figure2.pdf`` and ``.png``.

The source panel is retained unchanged; panel labels are added in fixed
normalized coordinates so that the caption and rendered panels correspond.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "paper/figures/main/fig6_syntracker_structure.png"
OUTPUT = ROOT / "paper/figures/extended_data/extended_data_figure2"

LABELS = {
    "a": (0.017, 0.976),
    "b": (0.517, 0.976),
    "c": (0.017, 0.505),
    "d": (0.517, 0.505),
}


def main() -> None:
    image = mpimg.imread(SOURCE)
    height, width = image.shape[:2]
    dpi = 300.0
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.imshow(image, interpolation="none")
    for label, (x, y) in LABELS.items():
        ax.text(
            x, y, label, transform=ax.transAxes,
            fontsize=13, fontweight="bold", ha="left", va="top",
        )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT.with_suffix(".pdf"), dpi=dpi)
    fig.savefig(OUTPUT.with_suffix(".png"), dpi=dpi)
    plt.close(fig)
    print(f"Generated {OUTPUT.with_suffix('.pdf')} and {OUTPUT.with_suffix('.png')}")


if __name__ == "__main__":
    main()
