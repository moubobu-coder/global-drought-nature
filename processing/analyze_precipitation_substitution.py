#!/usr/bin/env python
"""Compare Sen slopes of basin-network SSI and SPEI annual series.

The input CSV includes the SPEI series from the precipitation replacement
experiment using GHCD precipitation and ERA5 atmospheric evaporative demand.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Iterable

import numpy as np


SERIES = ("SSI", "GHCD", "ERA5", "JRA3Q", "MIX")
SPEI = ("GHCD", "ERA5", "JRA3Q", "MIX")
EXPECTED_YEARS = np.arange(1981, 2020, dtype=int)
BLOCK_LENGTHS = (3, 5, 7)
BASE_SEED = 20261017


def sen_slope_per_decade(years: np.ndarray, values: np.ndarray) -> float:
    """Median of all pairwise slopes, expressed per decade."""
    years = np.asarray(years, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if years.ndim != 1 or values.shape != years.shape or len(years) < 2:
        raise ValueError("Sen slope requires matching 1D arrays with at least 2 values")
    if not np.isfinite(years).all() or not np.isfinite(values).all():
        raise ValueError("Sen slope inputs must be finite")
    i, j = np.triu_indices(len(years), k=1)
    dx = years[j] - years[i]
    if np.any(dx == 0):
        raise ValueError("Years must be unique")
    return float(np.median((values[j] - values[i]) / dx) * 10.0)


def read_csv(path: Path, required: Iterable[str], *, exact_expected_years: bool = True) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"CSV has no data rows: {path}")
    try:
        years = np.array([int(row["year"]) for row in rows], dtype=int)
        series = {
            name: np.array([float(row[name]) for row in rows], dtype=np.float64)
            for name in required if name != "year"
        }
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Non-numeric or empty required value in {path}") from exc
    if len(np.unique(years)) != len(years):
        raise ValueError(f"Duplicate year in {path}")
    order = np.argsort(years)
    years = years[order]
    series = {name: values[order] for name, values in series.items()}
    if exact_expected_years and not np.array_equal(years, EXPECTED_YEARS):
        raise ValueError(f"{path} must contain every year 1981–2019 exactly once")
    if not exact_expected_years:
        if len(years) < 1 or np.any(np.diff(years) != 1):
            raise ValueError(f"{path} must contain one contiguous annual sequence")
        if years[0] > EXPECTED_YEARS[0] or years[-1] < EXPECTED_YEARS[-1]:
            raise ValueError(f"{path} must cover every year 1981–2019; additional contiguous years are allowed")
    if any(not np.isfinite(values).all() for values in series.values()):
        raise ValueError(f"Non-finite required value in {path}")
    return years, series


def window_table(years: np.ndarray, data: dict[str, np.ndarray]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    n = len(years)
    for start in range(n):
        for stop in range(start + 30, n + 1):
            wyears = years[start:stop]
            slopes = {name: sen_slope_per_decade(wyears, data[name][start:stop]) for name in SERIES}
            gaps = {name: abs(slopes[name] - slopes["SSI"]) for name in SPEI}
            row: dict[str, object] = {
                "start_year": int(wyears[0]), "end_year": int(wyears[-1]), "n_years": len(wyears),
                "sen_SSI_per_decade": slopes["SSI"],
            }
            for name in SPEI:
                row[f"sen_{name}_per_decade"] = slopes[name]
                row[f"gap_{name}"] = gaps[name]
            row["absolute_gap_decrease_MIX_vs_ERA5"] = gaps["ERA5"] - gaps["MIX"]
            row["percentage_gap_decrease_MIX_vs_ERA5"] = (
                None if gaps["ERA5"] == 0 else 100.0 * (1.0 - gaps["MIX"] / gaps["ERA5"])
            )
            rows.append(row)
    return rows


def detrended_correlations(years: np.ndarray, data: dict[str, np.ndarray]) -> dict[str, float]:
    x = years.astype(np.float64)
    x -= x.mean()
    residuals = {}
    for name, values in data.items():
        design = np.column_stack((np.ones(len(x)), x))
        beta, *_ = np.linalg.lstsq(design, values, rcond=None)
        residuals[name] = values - design @ beta
    return {name: float(np.corrcoef(residuals["SSI"], residuals[name])[0, 1]) for name in SPEI}


def deterministic_seed(block_length: int) -> int:
    # Fixed documented derivation, independent of Python's randomized hash().
    return BASE_SEED + block_length * 1009


def paired_circular_bootstrap(
    years: np.ndarray,
    data: dict[str, np.ndarray],
    reps: int,
    block_length: int,
    seed: int,
    chunk_size: int = 128,
) -> dict[str, np.ndarray]:
    """Joint circular moving-block residual bootstrap, processed in chunks."""
    names = list(SERIES)
    y = np.column_stack([data[name] for name in names])
    centered_year = years.astype(np.float64) - float(np.mean(years))
    design = np.column_stack((np.ones(len(years)), centered_year))
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ beta
    residual = y - fitted
    n, k = y.shape
    i, j = np.triu_indices(n, k=1)
    slopes_out = np.empty((reps, k), dtype=np.float64)
    rng = np.random.default_rng(seed)
    completed = 0
    while completed < reps:
        batch = min(chunk_size, reps - completed)
        starts = rng.integers(0, n, size=(batch, math.ceil(n / block_length)))
        offsets = np.arange(block_length, dtype=int)
        idx = (starts[:, :, None] + offsets[None, None, :]) % n
        idx = idx.reshape(batch, -1)[:, :n]
        sampled = residual[idx, :]
        synthetic = fitted[None, :, :] + sampled
        dy = synthetic[:, j, :] - synthetic[:, i, :]
        slopes_out[completed:completed + batch, :] = np.median(dy / (years[j] - years[i])[None, :, None], axis=1) * 10.0
        completed += batch

    ssi_index = names.index("SSI")
    result: dict[str, np.ndarray] = {f"sen_{name}_per_decade": slopes_out[:, idx] for idx, name in enumerate(names)}
    for name in SPEI:
        idx = names.index(name)
        signed = slopes_out[:, idx] - slopes_out[:, ssi_index]
        gap = np.abs(signed)
        result[f"signed_gap_{name}"] = signed
        result[f"gap_{name}"] = gap
    era_gap = result["gap_ERA5"]
    mix_gap = result["gap_MIX"]
    result["absolute_gap_decrease_MIX_vs_ERA5"] = era_gap - mix_gap
    ratio = np.full(reps, np.nan, dtype=np.float64)
    nonzero = era_gap != 0
    ratio[nonzero] = 1.0 - mix_gap[nonzero] / era_gap[nonzero]
    result["percentage_gap_decrease_MIX_vs_ERA5"] = ratio * 100.0
    return result


def percentile_interval(values: np.ndarray) -> dict[str, object]:
    valid = values[np.isfinite(values)]
    if not len(valid):
        return {"lower_2.5": None, "median": None, "upper_97.5": None,
                "valid_replicates": 0, "non_finite_replicates": int(values.size)}
    lo, med, hi = np.percentile(valid, [2.5, 50, 97.5])
    return {"lower_2.5": float(lo), "median": float(med), "upper_97.5": float(hi),
            "valid_replicates": int(len(valid)), "non_finite_replicates": int(values.size - len(valid))}


def bootstrap_summary(data: dict[str, np.ndarray], years: np.ndarray, reps: int) -> list[dict[str, object]]:
    output = []
    for block_length in BLOCK_LENGTHS:
        seed = deterministic_seed(block_length)
        draws = paired_circular_bootstrap(years, data, reps, block_length, seed)
        era5_gap_draws = draws["gap_ERA5"]
        denominator_near_zero_count = int(np.count_nonzero(era5_gap_draws < 1e-6))
        for metric, values in draws.items():
            output.append({"block_length_years": block_length, "replicates": reps, "seed": seed,
                           "metric": metric,
                           "denominator_near_zero_count_lt_1e-6": denominator_near_zero_count if metric == "percentage_gap_decrease_MIX_vs_ERA5" else None,
                           "near_zero_threshold_units": "absolute ERA5-SPEI to SSI gap < 1e-6 index units/decade; reported only, not filtered" if metric == "percentage_gap_decrease_MIX_vs_ERA5" else None,
                           **percentile_interval(values)})
    return output


def precipitation_summary(path: Path) -> dict[str, object]:
    years_all, p_all = read_csv(path, ("year", "P_GHCD", "P_ERA5"), exact_expected_years=False)
    mask = (years_all >= EXPECTED_YEARS[0]) & (years_all <= EXPECTED_YEARS[-1])
    years = years_all[mask]
    p = {name: values[mask] for name, values in p_all.items()}
    if not np.array_equal(years, EXPECTED_YEARS):
        raise ValueError("Selected precipitation records do not exactly cover 1981–2019")
    diff = p["P_ERA5"] - p["P_GHCD"]
    def mean_period(a: int, b: int) -> float:
        mask = (years >= a) & (years <= b)
        return float(np.mean(diff[mask]))
    return {
        "input_csv": str(path.resolve()),
        "input_years_preserved_in_source_file": [int(years_all[0]), int(years_all[-1])],
        "analysis_years_selected": [1981, 2019],
        "difference_definition": "P_ERA5 - P_GHCD",
        "full_period": {"start_year": 1981, "end_year": 2019, "mean_difference": float(np.mean(diff)),
                        "sen_difference_per_decade": sen_slope_per_decade(years, diff),
                        "sen_P_GHCD_per_decade": sen_slope_per_decade(years, p["P_GHCD"]),
                        "sen_P_ERA5_per_decade": sen_slope_per_decade(years, p["P_ERA5"])},
        "predeclared_period_means": {
            "1981-1996": {"mean_P_GHCD": float(np.mean(p["P_GHCD"][:16])), "mean_P_ERA5": float(np.mean(p["P_ERA5"][:16])), "mean_P_ERA5_minus_P_GHCD": mean_period(1981, 1996)},
            "2004-2019": {"mean_P_GHCD": float(np.mean(p["P_GHCD"][-16:])), "mean_P_ERA5": float(np.mean(p["P_ERA5"][-16:])), "mean_P_ERA5_minus_P_GHCD": mean_period(2004, 2019)},
        },
        "transition_interval_retained_in_full_period": [1997, 2003],
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_checks(years: np.ndarray, data: dict[str, np.ndarray], windows: list[dict[str, object]]) -> dict[str, object]:
    from scipy.stats import theilslopes
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from drought.stats import ktaub
    samples = [("SSI", 0, 39), ("GHCD", 0, 30), ("ERA5", 5, 39), ("MIX", 0, 39)]
    checks = []
    for name, a, b in samples:
        own = sen_slope_per_decade(years[a:b], data[name][a:b])
        scipy_slope = float(theilslopes(data[name][a:b], years[a:b]).slope * 10.0)
        ktaub_slope = float(ktaub(np.column_stack((years[a:b], data[name][a:b]))).sen * 10.0)
        if not np.isclose(own, scipy_slope, rtol=1e-12, atol=1e-12):
            raise AssertionError(f"Sen parity failure for {name} slice {a}:{b}: {own} vs {scipy_slope}")
        if not np.isclose(own, ktaub_slope, rtol=1e-12, atol=1e-12):
            raise AssertionError(f"KTAUB Sen parity failure for {name} slice {a}:{b}: {own} vs {ktaub_slope}")
        checks.append({"series": name, "slice_start_year": int(years[a]), "slice_end_year": int(years[b-1]),
                       "sen_own_per_decade": own, "scipy_theilslopes_per_decade": scipy_slope,
                       "project_ktaub_port_per_decade": ktaub_slope,
                       "absolute_difference": abs(own - scipy_slope)})
    if len(windows) != 55:
        raise AssertionError(f"Expected 55 contiguous >=30-year windows; found {len(windows)}")
    return {"scipy_theilslopes_parity": checks, "window_count": len(windows), "window_count_expected": 55}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, required=True, help="Annual global SSI/SPEI CSV with year,SSI,GHCD,ERA5,JRA3Q,MIX")
    parser.add_argument("--global-p-csv", type=Path, help="Optional annual global precipitation CSV with year,P_GHCD,P_ERA5")
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_output/precipitation_substitution"))
    parser.add_argument("--bootstrap-reps", type=int, default=5000)
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of existing output files")
    args = parser.parse_args()
    if args.bootstrap_reps < 1:
        parser.error("--bootstrap-reps must be positive")
    input_path = args.input_csv.resolve(strict=True)
    years, data = read_csv(input_path, ("year", *SERIES))
    windows = window_table(years, data)
    checks = run_checks(years, data, windows)

    slopes = {name: sen_slope_per_decade(years, data[name]) for name in SERIES}
    gaps = {name: abs(slopes[name] - slopes["SSI"]) for name in SPEI}
    signed_gaps = {name: slopes[name] - slopes["SSI"] for name in SPEI}
    era_gap = gaps["ERA5"]
    mix_gap = gaps["MIX"]
    summary: dict[str, object] = {
        "input_csv": str(input_path),
        "years": [int(years[0]), int(years[-1])],
        "n_years": int(len(years)),
        "sen_slopes_per_decade": slopes,
        "signed_trend_difference_s_SPEI_minus_s_SSI": signed_gaps,
        "absolute_trend_gap": gaps,
        "mix_vs_era5": {
            "signed_gap_change": signed_gaps["MIX"] - signed_gaps["ERA5"],
            "absolute_gap_decrease_primary": era_gap - mix_gap,
            "percentage_gap_decrease_secondary": None if era_gap == 0 else 100.0 * (1.0 - mix_gap / era_gap),
            "percentage_gap_decrease_units": "%",
            "percentage_denominator": "absolute ERA5-SPEI to SSI Sen-trend gap; undefined if zero",
        },
        "ghcd_vs_era5": {
            "absolute_gap_decrease": era_gap - gaps["GHCD"],
            "percentage_gap_decrease": None if era_gap == 0 else 100.0 * (1.0 - gaps["GHCD"] / era_gap),
            "percentage_gap_decrease_units": "%",
            "percentage_denominator": "absolute ERA5-SPEI to SSI Sen-trend gap",
        },
        "all_contiguous_windows": {"minimum_years": 30, "count": len(windows), "overlapping_windows_are_not_independent": True},
        "detrended_pearson_correlations_secondary": detrended_correlations(years, data),
        "conditional_temporal_uncertainty": {
            "method": "paired circular moving-block residual bootstrap around separate OLS intercept+fixed-year trends",
            "block_lengths_years": list(BLOCK_LENGTHS), "replicates_per_block_length": args.bootstrap_reps,
            "base_seed": BASE_SEED,
            "seed_derivation": "seed = 20261017 + block_length_years * 1009",
            "shared_resampling_indices_across_series": True,
            "interval": "percentiles 2.5, 50, 97.5",
            "percentage_gap_denominator_near_zero_threshold": 1e-6,
            "percentage_gap_denominator_near_zero_rule": "count ERA5 bootstrap gaps < 1e-6 index units per decade; report only and do not filter replicates",
            "percentage_gap_non_finite_replicates": "counted per block length in the bootstrap intervals CSV; non-finite percentage ratios are excluded from percentile calculation but are not otherwise filtered",
            "scope": "conditional temporal uncertainty only; excludes spatial/SSI fitting uncertainty; no bootstrap p-values",
        },
        "numerical_checks": checks,
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = [args.output_dir / "global_swap_statistics.json", args.output_dir / "all_contiguous_windows_ge30.csv",
             args.output_dir / "paired_residual_block_bootstrap_intervals.csv"]
    if args.global_p_csv:
        paths.append(args.output_dir / "global_precipitation_statistics.json")
    if not args.overwrite:
        existing = [str(p) for p in paths if p.exists()]
        if existing:
            raise FileExistsError("Output already exists; pass --overwrite to replace: " + ", ".join(existing))

    bootstrap_rows = bootstrap_summary(data, years, args.bootstrap_reps)
    summary["bootstrap_percentage_diagnostics_by_block_length"] = [
        {
            "block_length_years": int(row["block_length_years"]),
            "replicates": int(row["replicates"]),
            "era5_gap_denominator_near_zero_count_lt_1e-6": int(row["denominator_near_zero_count_lt_1e-6"]),
            "percentage_non_finite_replicate_count": int(row["non_finite_replicates"]),
            "near_zero_replicates_reported_but_not_filtered": True,
        }
        for row in bootstrap_rows if row["metric"] == "percentage_gap_decrease_MIX_vs_ERA5"
    ]
    precip = precipitation_summary(args.global_p_csv.resolve(strict=True)) if args.global_p_csv else None
    with paths[0].open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
    write_csv(paths[1], windows)
    write_csv(paths[2], bootstrap_rows)
    if precip is not None:
        with paths[3].open("w", encoding="utf-8") as stream:
            json.dump(precip, stream, ensure_ascii=False, indent=2, allow_nan=False)
    print(json.dumps({"status": "success", "output_files": [str(p.resolve()) for p in paths]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
