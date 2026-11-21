"""Manuscript figure ED5."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main():
    from drought.cli import figure_context
    ctx = figure_context("ED5")
    from drought.style import setup_style
    setup_style()
    from drought.map_figures import render_climatological_demand_spei_trends
    render_climatological_demand_spei_trends(ctx)


if __name__ == "__main__":
    main()
