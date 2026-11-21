"""Manuscript figure S1."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main():
    from drought.cli import figure_context
    ctx = figure_context("S1")
    from drought.style import setup_style
    setup_style()
    from drought.map_figures import render_ocean_land_relative_humidity
    render_ocean_land_relative_humidity(ctx)


if __name__ == "__main__":
    main()
