"""Shared input and export helpers for fixed-layout manuscript figures."""
from __future__ import annotations

import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def input_directory(ctx, figure_name):
    return ctx.path(figure_name + "_source", "source_data/" + figure_name)


def verify_inputs(directory, names):
    missing = [directory / name for name in names if not (directory / name).is_file()]
    if missing:
        raise FileNotFoundError("Missing figure inputs: " + ", ".join(str(path) for path in missing))
    return [directory / name for name in names]


def save_layout(fig, ctx):
    ctx.output.mkdir(parents=True, exist_ok=True)
    for suffix in ctx.formats:
        options = {}
        if suffix in ("tif", "tiff"):
            options["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(ctx.output / f"{ctx.figure_id}.{suffix}", dpi=ctx.dpi, **options)
    plt.close(fig)
    print(f"{ctx.figure_id}: saved to {ctx.output}")


COLORS = {"SSI-12": "#007FAF", "GHCD SPEI-12": "#DF5C4B",
          "ERA5 SPEI-12": "#4B7888", "JRA-3Q SPEI-12": "#EC8A1C"}
NAMES = {"SSI-12": "SSI", "GHCD SPEI-12": "GHCD",
         "ERA5 SPEI-12": "ERA5", "JRA-3Q SPEI-12": "JRA-3Q"}
MAP_BOUNDS = np.array([-150, -75, -25, -5, 5, 25, 75, 150], dtype=float)
MAP_COLORS = ["#543005", "#A66A20", "#DFC27D", "#F5F5F5",
              "#80CDC1", "#238F83", "#004C43"]


def p_label(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.4f}" if p < 0.01 else f"p = {p:.3f}"


def signed(value, digits):
    return f"{value:+.{digits}f}".replace("-", "−")


