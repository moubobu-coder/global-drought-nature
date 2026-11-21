"""Plot precipitation trends from GPCC, MSWEP and GPCP."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from drought.cli import figure_context


def main():
    ctx = figure_context("ED8")
    from drought.map_figures import render_precipitation_trends
    render_precipitation_trends(ctx)


if __name__ == "__main__":
    main()
