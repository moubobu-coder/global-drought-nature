"""Draw manuscript display ED9."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from drought.cli import figure_context


def main():
    ctx = figure_context("ED9")
    from drought.extended_data09 import render
    render(ctx)


if __name__ == "__main__":
    main()
