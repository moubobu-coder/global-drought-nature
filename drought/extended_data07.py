#!/usr/bin/env python3
"""Plot basin-network SSI and SPEI and the effect of replacing ERA5 precipitation."""
from __future__ import annotations

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from .manuscript_layout import input_directory, read_csv, save_layout, verify_inputs


SERIES = ("SSI", "ERA5", "GHCD", "MIX")
LABELS = ("SSI", "ERA5 SPEI", "GHCD SPEI", "ERA5 SPEI with GHCD precipitation")
COLORS = ("#3C3C3C", "#255C93", "#267455", "#AD510F")
BLOCKS = (3, 5, 7)
YEARS = np.arange(1981, 2020)
POINT_UNIT = "index units per decade"


def render(ctx):
    directory = input_directory(ctx, "figED07")
    annual, bootstrap = verify_inputs(directory, (
        "source_data_panel_a.csv", "source_data_panel_b.csv",
    ))
    annual_rows = read_csv(annual)
    values, trend_lines = {}, {}
    for key, label in zip(SERIES, LABELS):
        rows = sorted((row for row in annual_rows if row["series"] == label),
                      key=lambda row: int(row["year"]))
        if [int(row["year"]) for row in rows] != YEARS.tolist():
            raise ValueError(f"Expected 39 observations for {label}.")
        values[key] = np.array([float(row["annual_index"]) for row in rows])
        trend_lines[key] = np.array([float(row["sen_trend_line"]) for row in rows])
        if not np.isfinite([values[key], trend_lines[key]]).all():
            raise ValueError(f"Non-finite observation or Sen trend for {label}.")
    intervals = read_csv(bootstrap)
    if [int(row["block_length_years"]) for row in intervals] != list(BLOCKS):
        raise ValueError("Expected paired intervals for 3-, 5- and 7-year blocks.")
    estimate = float(intervals[0]["gap_reduction_index_units_per_decade"])
    for row in intervals:
        low = float(row["lower_95_index_units_per_decade"])
        high = float(row["upper_95_index_units_per_decade"])
        if not (0 < low < estimate < high and int(row["replicates"]) == 5000
                and float(row["gap_reduction_index_units_per_decade"]) == estimate):
            raise ValueError("Unexpected bootstrap interval or replicate count.")
    if ctx.check_inputs:
        print(f"{ctx.figure_id}: inputs checked (156 annual values and three intervals)")
        return
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "font.size": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.65,
        "axes.unicode_minus": True,
        "legend.frameon": False,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })
    width_mm = 183.0
    height_mm = 75.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))
    grid = fig.add_gridspec(1, 2, left=0.088, right=0.981, bottom=0.20, top=0.77,
                           width_ratios=(1.75, 1), wspace=0.44)
    left, right = (fig.add_subplot(grid[0, column]) for column in range(2))

    left.axhline(0, color="#BABABA", linestyle=(0, (1, 2.4)), linewidth=0.7, zorder=0)
    annual_handles = []
    for name, label, color in zip(SERIES, LABELS, COLORS):
        annual_handle, = left.plot(YEARS, values[name], color=color,
                                   linewidth=1.0 if name == "MIX" else 0.8,
                                   alpha=0.88, label=label, zorder=2)
        annual_handles.append(annual_handle)
        left.plot(YEARS, trend_lines[name], color=color, linewidth=1.15 if name == "MIX" else 0.95,
                  linestyle=(0, (4.5, 2.5)), zorder=3)
    legend_handles = [annual_handles[index] for index in (0, 2, 1, 3)]
    left.legend(handles=legend_handles, ncols=2, loc="lower left", bbox_to_anchor=(0, 1.055),
                borderaxespad=0, handlelength=1.55, handletextpad=0.55, columnspacing=1.3)
    left.set(xlim=(1980, 2020), xticks=[1981, 1990, 2000, 2010, 2019],
             ylim=(-0.60, 0.60), yticks=[-0.6, -0.3, 0, 0.3, 0.6],
             xlabel="Year", ylabel="Basin-network mean index")

    right.axvline(0, color="#8A8A8A", linestyle=(0, (2, 2)), linewidth=0.7, zorder=0)
    positions = np.arange(len(intervals))[::-1]
    for position, row in zip(positions, intervals):
        low, high = float(row["lower_95_index_units_per_decade"]), float(row["upper_95_index_units_per_decade"])
        right.errorbar(
            estimate, position, xerr=[[estimate - low], [high - estimate]],
            fmt="o", color=COLORS[-1], markersize=4.4, markeredgecolor="white",
            markeredgewidth=0.55, ecolor=COLORS[-1], elinewidth=1.2,
            capsize=2.4, capthick=0.95, zorder=3,
        )
    right.set(
        xlim=(-0.008, 0.157), xticks=[0, 0.05, 0.10, 0.15],
        ylim=(-0.60, 2.60), yticks=positions, yticklabels=[str(n) for n in BLOCKS],
        xlabel=f"Trend-gap reduction\n({POINT_UNIT})", ylabel="Bootstrap block (years)",
    )
    right.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    offset = mpl.transforms.ScaledTranslation(-0.23, 0.15, fig.dpi_scale_trans)
    for letter, axis in zip("ab", (left, right)):
        axis.text(0, 1, letter, transform=axis.transAxes + offset,
                  ha="left", va="bottom", fontweight="bold", fontsize=8)
        axis.tick_params(length=2.7, width=0.65, pad=2.6)

    actual_min = min(float(series.min()) for series in (*values.values(), *trend_lines.values()))
    actual_max = max(float(series.max()) for series in (*values.values(), *trend_lines.values()))
    if not left.get_ylim()[0] < actual_min <= actual_max < left.get_ylim()[1]:
        raise ValueError(f"Panel a would clip data: min {actual_min}, max {actual_max}.")

    save_layout(fig, ctx)
