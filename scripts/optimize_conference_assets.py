#!/usr/bin/env python3
"""Regenerate presentation assets from vector figures using Arial and no outer border."""
from __future__ import annotations
from pathlib import Path
import shutil
import subprocess
from PIL import Image, ImageOps

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "conference_assets"
MAPPING = {
    "figure1.pdf": "paper/figures/main/fig1_algorithm.pdf",
    "figure2.pdf": "paper/figures/main/fig2_simulations.pdf",
    "figure3.pdf": "paper/figures/main/fig3_gtdb_r207_benchmark.pdf",
    "figure4.pdf": "paper/figures/main/fig5_sv_validation.pdf",
    "figure5.pdf": "paper/figures/main/fig7_b_longum_abfa.pdf",
    "figure6.pdf": "paper/figures/main/fig8_cagpai_h_pylori.pdf",
    "extended_data_figure1.pdf": "paper/figures/main/fig4_efficiency.pdf",
    "extended_data_figure3.pdf": "results/census/hrom/hypermodes/fig_hrom_hypermodes.pdf",
    "extended_data_figure4.pdf": "results/census/oapgc/figures/fig_oapgc_site_axes.pdf",
    "extended_data_figure5.pdf": "results/census/gtdb/figures/fig_gtdb_ani_sv_census.pdf",
}


def optimize_pdf(pdf: Path, out: Path, dpi: int = 400) -> None:
    tmp = OUT / f".{out.stem}.raw"
    subprocess.run(
        ["pdftoppm", "-png", "-singlefile", "-r", str(dpi), str(pdf), str(tmp)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    raster = OUT / f".{out.stem}.raw.png"
    with Image.open(raster) as im:
        rgb = im.convert("RGB")
        # Remove surrounding white space so no page border or blank frame remains.
        gray = ImageOps.grayscale(rgb)
        content = ImageOps.invert(gray).point(lambda p: 255 if p > 6 else 0)
        bbox = content.getbbox()
        if bbox:
            margin = int(dpi * 0.035)  # ~0.035 inch of white margin
            left = max(0, bbox[0] - margin)
            top = max(0, bbox[1] - margin)
            right = min(rgb.width, bbox[2] + margin)
            bottom = min(rgb.height, bbox[3] + margin)
            rgb = rgb.crop((left, top, right, bottom))
        rgb.save(out, format="PNG", optimize=True, dpi=(dpi, dpi))
    raster.unlink(missing_ok=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for out_name, source in MAPPING.items():
        source_path = ROOT / source
        if not source_path.exists():
            raise FileNotFoundError(source_path)
        out = OUT / out_name.replace(".pdf", ".png")
        optimize_pdf(source_path, out, dpi=400)
        print(f"optimized {source} -> {out}")
        # Keep a lossless vector-bearing PDF beside the PNG for presentation export.
        shutil.copyfile(source_path, OUT / out_name)


if __name__ == "__main__":
    main()
