"""Draw manuscript display ED7."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from drought.cli import figure_context


def main():
    ctx = figure_context("ED7")
    from drought.extended_data07 import render
    render(ctx)


if __name__ == "__main__":
    main()
