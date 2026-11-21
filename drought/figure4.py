"""Plot hemispheric SPEI and water storage and basin drought indices.

The basin source and CSR trend map are retained. Regional summaries are read
from a separately calculated fixed common domain; no SPEI is refitted here.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
from shapely.geometry import box

from .manuscript_layout import (
    COLORS, MAP_BOUNDS, MAP_COLORS, NAMES, input_directory, p_label,
    read_csv, read_json, save_layout, signed, verify_inputs,
)

WIDTH_MM = 200
HEIGHT_MM = 165
LEFT_X, LEFT_W, RIGHT_X, RIGHT_W = 0.075, 0.40, 0.555, 0.40
REGIONAL_COLORS = {
    ("SPEI12", "NH"): "#DF5C4B", ("SPEI12", "SH"): "#DF5C4B",
    ("TWSA", "NH"): "#007F73", ("TWSA", "SH"): "#007F73",
}


def rect_mm(left, bottom, width, height):
    return [left / WIDTH_MM, bottom / HEIGHT_MM, width / WIDTH_MM, height / HEIGHT_MM]


def header(fig, ax, letter, label, left, width, y):
    ax.text(left - 0.045, y, letter, transform=fig.transFigure,
            ha="left", va="bottom", fontsize=10, fontweight="bold", clip_on=False)
    return ax.text(left + width / 2, y, label, transform=fig.transFigure,
                   ha="center", va="bottom", fontsize=8.3, clip_on=False)


def regional_inputs(annual_path):
    source_rows = read_csv(annual_path)
    series = {}
    for variable in ("SPEI12", "TWSA"):
        for region in ("NH", "SH"):
            selected = sorted([r for r in source_rows if r["variable"] == variable and r["region"] == region],
                              key=lambda r: int(r["year"]))
            value = {k: np.array([float(r[k]) if r[k] else np.nan for r in selected])
                     for k in ("year", "mean", "p25", "p75", "sen_line")}
            assert np.array_equal(value["year"], np.arange(2002, 2025))
            complete = np.isfinite(value["mean"])
            assert np.all(value["p25"][complete] <= value["p75"][complete])
            assert complete.sum() == (23 if variable == "SPEI12" else 21)
            assert np.isfinite(value["sen_line"]).all()
            if variable == "TWSA":
                assert np.array_equal(value["year"][~complete], [2017, 2018])
            series[variable, region] = value
    return series


def render(ctx):
    directory = input_directory(ctx, "fig04")
    basin_path, hemisphere_path, chart_path, stats_path = verify_inputs(directory, (
        "basin_SSI_SPEI_annual.csv", "hemisphere_GHCD_SPEI_TWS_annual.csv",
        "CSR_TWSA_Sen_2002_2024_common_land.nc", "figure_statistics.json",
    ))
    statistics = read_json(stats_path)
    basin_stats = statistics["basin"]
    regional_stats = statistics["regional"]
    regional = regional_inputs(hemisphere_path)
    years = np.arange(1981, 2020)
    basin_rows = read_csv(basin_path)
    source = {"series": {}}
    for key, name in NAMES.items():
        rows = sorted((row for row in basin_rows if row["series"] == name),
                      key=lambda row: int(row["year"]))
        if [int(row["year"]) for row in rows] != years.tolist():
            raise ValueError(f"Expected 39 basin-network observations for {name}.")
        source["series"][key] = {
            label: np.array([float(row[column]) for row in rows])
            for label, column in (("mean", "mean"), ("p25", "weighted_p25"), ("p75", "weighted_p75"))
        }
        basin_stats[key]["sen_line"] = np.array([float(row["sen_line"]) for row in rows])
    with xr.open_dataset(chart_path) as dataset:
        lon, lat = dataset.lon.values, dataset.lat.values
        slope = dataset.sen_slope.where(dataset.valid_mask > 0).values
    if int(np.isfinite(slope).sum()) != 160310:
        raise ValueError("Unexpected number of finite CSR trend grid cells.")

    if ctx.check_inputs:
        print(f"{ctx.figure_id}: inputs checked (39 basin years, two hemispheres, 160310 CSR cells)")
        return
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8.3,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.65, "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    a = fig.add_axes(rect_mm(15, 17, 80, 37), label="c")
    b = fig.add_axes(rect_mm(111, 17, 80, 37), label="d")
    c = fig.add_axes(rect_mm(111, 81, 80, 70), label="b")
    c.set_axis_off()
    map_projection = ccrs.EckertIV(central_longitude=0)
    map_y60 = map_projection.transform_point(0, 60, ccrs.PlateCarree())[1]
    natural_map_height_mm = 80 * 2 * map_y60 / np.ptp(map_projection.x_limits)
    c_map = fig.add_axes(rect_mm(111, 103, 80, natural_map_height_mm),
                         projection=map_projection, label="<map_b>")
    # Containers describe the upper modules; their visible data axes are also
    # measured and checked independently after the final layout draw.
    d = fig.add_axes(rect_mm(15, 81, 80, 70), label="a")
    d.set_axis_off()
    header(fig, d, "a", "Hemispheric SPEI and water storage", LEFT_X, LEFT_W, 154 / HEIGHT_MM)
    c_title = header(fig, c, "b", "Water storage trends", RIGHT_X, RIGHT_W, 154 / HEIGHT_MM)
    header(fig, a, "c", "12-month SSI", LEFT_X, LEFT_W, 66 / HEIGHT_MM)
    header(fig, b, "d", "12-month SPEI", RIGHT_X, RIGHT_W, 66 / HEIGHT_MM)
    for ax, keys in ((a, ["SSI-12"]), (b, ["GHCD SPEI-12", "ERA5 SPEI-12", "JRA-3Q SPEI-12"])):
        ax.axhline(0, color="0.55", linestyle=(0, (3, 2)), linewidth=0.65, zorder=0)
        for key in keys:
            value = source["series"][key]
            ax.fill_between(years, value["p25"], value["p75"], color=COLORS[key],
                            alpha=0.14, linewidth=0, zorder=0)
            ax.plot(years, value["mean"], color=COLORS[key], linewidth=1.2, zorder=3)
            ax.plot(years, basin_stats[key]["sen_line"], color=COLORS[key],
                    linewidth=1, linestyle=(0, (5, 3)), zorder=2)
        text_ys = [61.0] if len(keys) == 1 else [64.3, 61.0, 57.7]
        for key, y in zip(keys, text_ys):
            info = basin_stats[key]
            centre = 55 if ax is a else 151
            ax.text(centre / WIDTH_MM, y / HEIGHT_MM,
                    f"{NAMES[key]}: {signed(info['sen_per_decade'], 4)} decade$^{{-1}}$, {p_label(info['MK_p'])}",
                    transform=fig.transFigure, color=COLORS[key], fontsize=7.2, ha="center", va="top", clip_on=False)
        bound = 1.5 if ax is a else 1.0
        ax.set(xlim=(1980, 2020), ylim=(-bound, bound), xlabel="Year",
               ylabel="Area-weighted index (z units)")
        ax.set_xticks(np.arange(1980, 2021, 10))
        ax.set_yticks(np.arange(-bound, bound + 0.01, 0.5))
        ax.yaxis.labelpad = 4
    c_map.set_global()
    c_map.set_anchor("N")
    visible_geo_box = box(-180, -60, 180, 60)
    land_geometry = [value.intersection(visible_geo_box)
                     for value in cfeature.LAND.with_scale("110m").geometries()
                     if value.intersects(visible_geo_box)]
    coast_geometry = [value.intersection(visible_geo_box)
                      for value in cfeature.COASTLINE.with_scale("110m").geometries()
                      if value.intersects(visible_geo_box)]
    c_map.add_geometries(land_geometry, crs=ccrs.PlateCarree(), facecolor="#E7E7E7", edgecolor="none", zorder=0)
    cmap, norm = ListedColormap(MAP_COLORS, name="storage_change"), BoundaryNorm(MAP_BOUNDS, len(MAP_COLORS))
    mesh = c_map.pcolormesh(lon, lat, np.ma.masked_invalid(slope), cmap=cmap, norm=norm,
                        transform=ccrs.PlateCarree(), shading="auto", rasterized=True, zorder=1)
    y60 = c_map.projection.transform_point(0, 60, ccrs.PlateCarree())[1]
    c_map.set_xlim(*c_map.projection.x_limits)
    c_map.set_ylim(-y60, y60)
    c_map.add_geometries(coast_geometry, crs=ccrs.PlateCarree(), facecolor="none",
                         edgecolor="#454545", linewidth=0.35, zorder=2)
    c_map.spines["geo"].set_linewidth(0.6)
    gl = c_map.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, xlocs=[-120, 0, 120],
                     ylocs=[-50, 0, 50], linewidth=0, x_inline=False, y_inline=False)
    gl.top_labels, gl.right_labels, gl.left_labels = False, False, False
    gl.xlabel_style = {"size": 6.2, "color": "0.3"}
    gl.xpadding = 3
    for latitude, label in [(-50, "50°S"), (0, "0°"), (50, "50°N")]:
        projected_y = c_map.projection.transform_point(0, latitude, ccrs.PlateCarree())[1]
        c_map.text(-0.025, 0.5 + projected_y / (2 * y60), label, transform=c_map.transAxes,
               ha="right", va="center", fontsize=6.2, color="0.3", clip_on=False)
    cb_ax = fig.add_axes(rect_mm(114.6, 91.3, 72.8, 2.55), label="<colorbar>")
    cb = fig.colorbar(mesh, cax=cb_ax, orientation="horizontal", extend="both", ticks=MAP_BOUNDS, spacing="uniform")
    cb.ax.tick_params(labelsize=6.5, length=2, pad=2)
    cb.outline.set_linewidth(0.55)
    cb.set_label("Water storage trend (mm decade$^{-1}$)", fontsize=7.5, labelpad=4)
    handles = [Line2D([], [], color="#DF5C4B", linewidth=1.2, label="GHCD SPEI"),
               Line2D([], [], color="#007F73", linewidth=1.2, label="TWS")]
    indicator_legend = d.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, (149.5 - 81) / 70),
             frameon=False, ncol=2, fontsize=6.8, columnspacing=1.6,
             handlelength=2.2, borderaxespad=0)
    split_axes, twin_axes = [], []
    spei_max, tws_max = 0.8, 30.0
    outside_means = {}
    for region, bottom, region_title_y, text_top in (("NH", 117, 139.5, 143.8),
                                                    ("SH", 81, 104.5, 107.8)):
        axis = fig.add_axes(rect_mm(15, bottom, 80, 22), label="<region_" + region + ">")
        twin = axis.twinx()
        twin.set_label("<twin_TWSA_" + region + ">")
        twin.spines["right"].set_visible(True)
        twin.spines["top"].set_visible(False)
        twin.patch.set_visible(False)
        axis.tick_params(axis="y", labelsize=6.5, pad=2)
        twin.tick_params(axis="y", labelsize=6.5, pad=2)
        for variable, target in (("SPEI12", axis), ("TWSA", twin)):
            value = regional[variable, region]
            color = REGIONAL_COLORS[variable, region]
            target.plot(value["year"], value["mean"], color=color, linewidth=1.15,
                        clip_on=False, zorder=3)
            target.plot(value["year"], value["sen_line"], color=color, linewidth=0.9,
                        linestyle=(0, (4, 2.5)), zorder=2)
            bound = spei_max if variable == "SPEI12" else tws_max
            assert np.nanmax(np.abs(value["sen_line"])) < bound
            above, below = value["mean"] > bound, value["mean"] < -bound
            outside_means[region + "_" + variable] = {
                "years": value["year"][above | below].astype(int).tolist(),
                "values": value["mean"][above | below].tolist(),
            }
        axis.set(xlim=(2001, 2025), ylim=(-spei_max, spei_max))
        twin.set_ylim(-tws_max, tws_max)
        axis.set_yticks([-0.8, -0.4, 0, 0.4, 0.8])
        twin.set_yticks([-30, -15, 0, 15, 30])
        axis.set_xticks([2005, 2010, 2015, 2020, 2025])
        axis.axhline(0, color="0.6", linewidth=0.6, linestyle=(0, (3, 2)), zorder=1)
        if region == "NH":
            axis.tick_params(axis="x", labelbottom=False)
        else:
            axis.set_xlabel("Year")
        name = "Northern Hemisphere" if region == "NH" else "Southern Hemisphere"
        axis.text(55 / WIDTH_MM, region_title_y / HEIGHT_MM, name,
                  transform=fig.transFigure, ha="center", va="bottom", fontsize=7.4)
        for i, variable in enumerate(("SPEI12", "TWSA")):
            info = regional_stats["regions"][region][variable]
            name = "SPEI" if variable == "SPEI12" else "TWS"
            label = f"{name}: {signed(info['sen_per_decade'], 3 if variable == 'SPEI12' else 2)}, {p_label(info['MK_p'])}"
            axis.text((15 + i * 41) / WIDTH_MM, text_top / HEIGHT_MM, label,
                      transform=fig.transFigure, color=REGIONAL_COLORS[variable, region],
                      fontsize=7.2, va="bottom", ha="left")
        split_axes.append(axis)
        twin_axes.append(twin)
    left_axis_title = d.text(5 / WIDTH_MM, 110 / HEIGHT_MM, "GHCD SPEI (z units)", transform=fig.transFigure,
           ha="center", va="center", rotation=90, fontsize=7.5)
    right_axis_title = d.text(102 / WIDTH_MM, 110 / HEIGHT_MM, "Water storage anomaly (mm)", transform=fig.transFigure,
           ha="center", va="center", rotation=90, fontsize=7.5)
    style = [Line2D([], [], color="0.4", linewidth=1, linestyle=(0, (5, 3)), label="Sen trend (per decade)"),
             Patch(facecolor="0.7", alpha=0.35, edgecolor="none", label="Area-weighted 25th–75th percentile (c, d)")]
    fig.legend(handles=style, loc="lower center", bbox_to_anchor=(0.52, 0.012), ncol=2,
               frameon=False, fontsize=7, handlelength=2.6, columnspacing=1.8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def tight_group(artists):
        boxes = [artist.get_tightbbox(renderer) if hasattr(artist, "get_tightbbox")
                 else artist.get_window_extent(renderer) for artist in artists]
        return Bbox.union([value for value in boxes if value is not None])

    d_content = split_axes + twin_axes + [indicator_legend, left_axis_title, right_axis_title]
    d_content_box = tight_group(d_content)
    c_content_box = tight_group([c_map, cb_ax])
    map_shift_mm = ((d_content_box.y0 + d_content_box.y1) -
                    (c_content_box.y0 + c_content_box.y1)) * 0.5 * 25.4 / fig.dpi
    for target in (c_map, cb_ax):
        position = target.get_position()
        target.set_position([position.x0, position.y0 + map_shift_mm / HEIGHT_MM,
                             position.width, position.height])
    fig.canvas.draw()
    c_title.set_y(c_map.get_position().y1 + 3 / HEIGHT_MM)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    c_final_box, d_final_box = tight_group([c_map, cb_ax]), tight_group(d_content)
    centre_error_mm = abs((c_final_box.y0 + c_final_box.y1) -
                          (d_final_box.y0 + d_final_box.y1)) * 0.5 * 25.4 / fig.dpi
    assert centre_error_mm < 0.01
    save_layout(fig, ctx)
