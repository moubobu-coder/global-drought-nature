"""Plot the relative-humidity and wind-speed station networks."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from drought.cli import figure_context


def main():
    ctx = figure_context("S4")
    from drought.style import setup_style
    from drought.map_figures import render_station_networks
    setup_style()
    render_station_networks(ctx)


if __name__ == "__main__":
    main()
