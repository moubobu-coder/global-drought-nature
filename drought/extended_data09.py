"""Extended Data Fig. 9 with African water-storage and GHCD SPEI time series."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
from shapely.geometry import box

from .manuscript_layout import input_directory, read_csv, read_json, save_layout, verify_inputs

WIDTH_MM, HEIGHT_MM = 190, 190
PRECIP_COLORS = {"GHCD": "#0072B2", "GPCC": "#CC79A7"}
SERIES_COLORS = {"SPEI12": "#DF5C4B", "TWSA": "#007F73"}
MAP_BOUNDS = np.array([-150, -75, -25, -5, 5, 25, 75, 150], float)
MAP_COLORS = ["#543005", "#A66A20", "#DFC27D", "#F5F5F5", "#80CDC1", "#238F83", "#004C43"]


def rect_mm(left, bottom, width, height):
    return [left / WIDTH_MM, bottom / HEIGHT_MM, width / WIDTH_MM, height / HEIGHT_MM]


def header(fig, ax, letter, label, left, bottom):
    ax.text((left - 9.9) / WIDTH_MM, bottom / HEIGHT_MM, letter,
            transform=fig.transFigure, ha="left", va="bottom", fontsize=10,
            fontweight="bold", clip_on=False)
    ax.text((left + 33.3) / WIDTH_MM, bottom / HEIGHT_MM, label,
            transform=fig.transFigure, ha="center", va="bottom", fontsize=8.3,
            clip_on=False)


def load_series(path):
    rows = [row for row in read_csv(path) if row["region"] == "Africa"]
    arrays = {}
    for variable in ("SPEI12", "TWSA"):
        take = sorted((r for r in rows if r["variable"] == variable), key=lambda r: int(r["year"]))
        assert len(take) == 23 and [int(r["year"]) for r in take] == list(range(2002, 2025))
        arrays[variable] = {key: np.array([float(r[key]) if r[key] else np.nan for r in take])
                            for key in ("year", "mean", "p25", "p75", "sen_line")}
        dat = arrays[variable]
        good = np.isfinite(dat["mean"])
        assert np.all(dat["p25"][good] <= dat["p75"][good])
        assert good.sum() == (23 if variable == "SPEI12" else 21)
        if variable == "TWSA":
            assert np.array_equal(dat["year"][~good], [2017, 2018])
    return arrays


def signed(x, digits):
    return f"{x:+.{digits}f}".replace("-", "−")


def p_label(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.4f}" if p < 0.01 else f"p = {p:.3f}"


def nice_step(value):
    assert np.isfinite(value) and value > 0
    magnitude = 10 ** np.floor(np.log10(value))
    return next(s * magnitude for s in (1, 2, 2.5, 5, 10) if s * magnitude >= value)


def render(ctx):
    directory = input_directory(ctx, "figED09")
    annual, statistics, china_path, storage_map = verify_inputs(directory, (
        "Africa_annual_SPEI_TWS.csv", "Africa_trend_statistics.json",
        "china_precipitation_original.csv", "africa_TWSA_Sen_common.nc",
    ))
    series = load_series(annual)
    stats = read_json(statistics)
    china = read_csv(china_path)
    if len(china) != 124:
        raise ValueError("Expected 124 annual precipitation entries from China.")
    with xr.open_dataset(storage_map) as dataset:
        lon, lat = dataset.lon.values, dataset.lat.values
        slope = dataset.sen_slope_mm_decade.where(dataset.valid_mask > 0).values
    if int(np.isfinite(slope).sum()) != 27313:
        raise ValueError("Unexpected number of finite African CSR trend grid cells.")
    if ctx.check_inputs:
        print(f"{ctx.figure_id}: inputs checked (124 precipitation entries, 27313 CSR cells)")
        return

    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8.3,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.65, "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    a = fig.add_axes(rect_mm(18, 134.2, 66.6, 41.4), label="a")
    b = fig.add_axes(rect_mm(108, 134.2, 66.6, 41.4), label="b")
    c = fig.add_axes(rect_mm(18, 41.4, 66.6, 66.6), projection=ccrs.PlateCarree(), label="c")
    d = fig.add_axes(rect_mm(108, 41.4, 66.6, 66.6), label="d")
    d_tws = d.twinx()
    d_tws.set_label("d_tws_axis")
    header(fig, a, "a", "Precipitation over China", 18, 180.46)
    header(fig, b, "b", "Precipitation over China: recent period", 108, 180.46)
    header(fig, c, "c", "Water storage trends in Africa", 18, 113.04)
    header(fig, d, "d", "Water storage and 12-month SPEI", 108, 113.04)

    for ax, panel in ((a, "a"), (b, "b")):
        for product in ("GHCD", "GPCC"):
            rows = sorted([r for r in china if r["panel"] == panel and r["product"] == product], key=lambda r: int(r["year"]))
            years = np.array([int(r["year"]) for r in rows])
            values = np.array([float(r["precipitation_mm_year"]) for r in rows])
            assert np.array_equal(years, np.arange(1973, 2025) if panel == "a" else np.arange(2015, 2025))
            ax.plot(years, values, color=PRECIP_COLORS[product], linewidth=1.3, label=product)
        ax.set(ylim=(600, 1250), xlabel="Year")
        ax.set_yticks([600, 800, 1000, 1200])
        ax.set_ylabel("Precipitation (mm year$^{-1}$)", fontsize=7.5, labelpad=4)
        ax.set_xlim((1973, 2024) if panel == "a" else (2015, 2024))
        ax.set_xticks([1980, 1990, 2000, 2010, 2020] if panel == "a" else [2015, 2018, 2021, 2024])
        ax.legend(loc="upper left", frameon=False, fontsize=7, ncol=2,
                  borderaxespad=0.4, handlelength=2.1, columnspacing=1.4)

    region_box = box(-20, -40, 60, 40)
    land = [g.intersection(region_box) for g in cfeature.LAND.with_scale("50m").geometries() if g.intersects(region_box)]
    coast = [g.intersection(region_box) for g in cfeature.COASTLINE.with_scale("50m").geometries() if g.intersects(region_box)]
    c.add_geometries(land, crs=ccrs.PlateCarree(), facecolor="#E7E7E7", edgecolor="none", zorder=0)
    mesh = c.pcolormesh(lon, lat, np.ma.masked_invalid(slope), transform=ccrs.PlateCarree(),
                       cmap=ListedColormap(MAP_COLORS), norm=BoundaryNorm(MAP_BOUNDS, 7),
                       shading="auto", rasterized=True, zorder=1)
    c.set_extent([-20, 60, -40, 40], crs=ccrs.PlateCarree())
    c.add_geometries(coast, crs=ccrs.PlateCarree(), facecolor="none", edgecolor="#454545", linewidth=0.4, zorder=2)
    c.spines["geo"].set_linewidth(0.6)
    gl = c.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, xlocs=[-20, 0, 20, 40, 60],
                     ylocs=[-40, -20, 0, 20, 40], linewidth=0, x_inline=False, y_inline=False)
    gl.top_labels, gl.right_labels = False, False
    gl.xlabel_style = {"size": 6.5, "color": "0.3"}
    gl.ylabel_style = {"size": 6.5, "color": "0.3"}
    gl.xpadding, gl.ypadding = 3, 3
    cb_ax = fig.add_axes(rect_mm(18, 24.3, 66.6, 2.7), label="<colorbar>")
    cb = fig.colorbar(mesh, cax=cb_ax, orientation="horizontal", extend="both", ticks=MAP_BOUNDS, spacing="uniform")
    cb.ax.tick_params(labelsize=6.5, length=2, pad=2)
    cb.outline.set_linewidth(0.55)
    cb.set_label("Water storage trend (mm decade$^{-1}$)", fontsize=7.5, labelpad=4)

    # Different physical units remain on labelled symmetric axes.
    range_coverage = {}
    for variable, ax in (("SPEI12", d), ("TWSA", d_tws)):
        dat, color = series[variable], SERIES_COLORS[variable]
        values = np.concatenate([dat["mean"], dat["p25"], dat["p75"], dat["sen_line"]])
        limit = 2.0 if variable == "SPEI12" else 200.0
        assert np.nanmax(np.abs(values)) < limit
        range_coverage[variable] = {"minimum": float(np.nanmin(values)),
                                    "maximum": float(np.nanmax(values)), "fully_in_range": True}
        ax.set_ylim(-limit, limit)
        ax.set_yticks(np.linspace(-limit, limit, 5))
        ax.fill_between(dat["year"], dat["p25"], dat["p75"], color=color, alpha=0.14, linewidth=0, zorder=1)
        ax.plot(dat["year"], dat["mean"], color=color, linewidth=1.2,
                marker="o" if variable == "SPEI12" else "s", markersize=2.0,
                markeredgewidth=0, zorder=3)
        ax.plot(dat["year"], dat["sen_line"], color=color, linewidth=1.0, linestyle=(0, (5, 3)), zorder=2)
    d_tws.spines["right"].set_visible(True)
    d_tws.spines["left"].set_visible(False)
    d_tws.patch.set_visible(False)
    d_tws.set_zorder(d.get_zorder() + 1)
    d_tws.set_ylabel("Water storage anomaly (mm)", color=SERIES_COLORS["TWSA"], labelpad=5)
    d_tws.tick_params(axis="y", colors=SERIES_COLORS["TWSA"])
    d.set(xlim=(2001, 2025), xlabel="Year")
    d.set_xticks([2005, 2010, 2015, 2020, 2025])
    d.set_ylabel("12-month SPEI (z units)", color=SERIES_COLORS["SPEI12"], labelpad=4)
    d.tick_params(axis="y", colors=SERIES_COLORS["SPEI12"])
    handles = [Line2D([], [], color=SERIES_COLORS[v], lw=1.2, label=label)
               for v, label in (("SPEI12", "GHCD SPEI"), ("TWSA", "GRACE TWS"))]
    for variable, y in (("SPEI12", 0.885), ("TWSA", 0.835)):
        info = stats[variable]
        units = "decade$^{-1}$" if variable == "SPEI12" else "mm decade$^{-1}$"
        label = f"{signed(info['sen_per_decade'], 4 if variable == 'SPEI12' else 2)} {units}, {p_label(info['MK_p'])}"
        d.text(0.04, y, label, transform=d.transAxes, color=SERIES_COLORS[variable],
               fontsize=7.2, ha="left", va="top", zorder=6)
    styles = [Line2D([], [], color="0.4", lw=1, linestyle=(0, (5, 3)), label="Sen trend"),
              Patch(facecolor="0.7", edgecolor="none", alpha=0.35, label="Area-weighted 25th–75th percentile")]
    # The series key is inside d; the line/band styles sit below the plot.
    series_legend = d.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.01, 0.995),
                             frameon=False, ncol=2, fontsize=6.7, handlelength=1.8, columnspacing=1.0)
    fig.legend(handles=styles, loc="upper center", bbox_to_anchor=(141.3 / WIDTH_MM, 26 / HEIGHT_MM),
               ncol=1, frameon=False, fontsize=6.7, handlelength=2.3, labelspacing=0.65)

    save_layout(fig, ctx)
