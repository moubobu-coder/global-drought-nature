# Global drought analysis and plotting code

This repository contains analysis and plotting code for the global drought study. Use Python 3.12 and install the dependencies in `requirements.txt`.

```sh
python -m pip install -r requirements.txt
python figures/fig04.py --input-root /path/to/inputs --output-dir figure_output --dpi 600
python figures/figED07.py --input-root /path/to/inputs --check-inputs
```

Each figure command accepts `--input-root`, `--output-dir`, `--config`, `--formats`, `--dpi` and `--check-inputs`. A JSON configuration can override input paths in its `paths` object. Relative paths are resolved under the input root. Cartopy uses Natural Earth boundaries for maps.

## Figure commands

| Manuscript figure | Script | Content |
| --- | --- | --- |
| Fig. 1 | `figures/fig01.py` | Annual meteorological drivers |
| Fig. 2 | `figures/fig02.py` | SPEI trends, differences and global annual means |
| Fig. 3 | `figures/fig03.py` | Drought-area trends |
| Fig. 4 | `figures/fig04.py` | Hemispheric and global water storage, and basin SSI and SPEI |
| Extended Data Figs. 1–4 | `figures/figED01.py` through `figures/figED04.py` | Driver trends and continental climate series |
| Extended Data Fig. 5 | `figures/figED05.py` | SPEI trends with climatological atmospheric evaporative demand |
| Extended Data Fig. 6 | `figures/figED06.py` | Continental drought-area trends |
| Extended Data Fig. 7 | `figures/figED07.py` | Basin SSI, three SPEI series and the precipitation replacement comparison |
| Extended Data Fig. 8 | `figures/figED08.py` | GPCC, MSWEP and GPCP precipitation trends |
| Extended Data Fig. 9 | `figures/figED09.py` | China precipitation and African water storage and SPEI |
| Supplementary Fig. 1 | `figures/figS01.py` | Land and ocean relative humidity |
| Supplementary Fig. 2 | `figures/figS02.py` | Continental annual SPEI |
| Supplementary Fig. 3 | `figures/figS03.py` | Winter Niño 3.4–SPEI correlations |
| Supplementary Fig. 4 | `figures/figS04.py` | Relative-humidity and wind-speed station networks |

## Data sources

- [GHCD meteorological fields and SPEI](https://doi.org/10.11888/Atmos.tpdc.303673)
- [ERA5](https://cds.climate.copernicus.eu/datasets) and [JRA-3Q](https://www.data.jma.go.jp/jra/html/JRA-3Q/index_en.html)
- [GPCC](https://opendata.dwd.de/climate_environment/GPCC/), [MSWEP](https://www.gloh2o.org/mswep/) and [GPCP v2.3](https://doi.org/10.7289/V56971M6)
- [CERES EBAF 4.2.1](https://doi.org/10.5067/TERRA-AQUA-NOAA20/CERES/EBAF_L3B004.2.1) and [SSRIH](https://doi.org/10.6084/m9.figshare.21625079.v1)
- [Niño 3.4 index](https://psl.noaa.gov/data/timeseries/month/Nino34/)
- [CSR GRACE and GRACE-FO Release 6.3](https://www2.csr.utexas.edu/grace/RL06_mascons.html)

## Prepared figure inputs

Use the study's prepared series, masks, trend grids and statistics for the figure commands. The following paths are relative to the input root.

| Figure | Input directory or key | Required files |
| --- | --- | --- |
| Fig. 4 | `source_data/fig04` (`fig04_source`) | `basin_SSI_SPEI_annual.csv`, `hemisphere_GHCD_SPEI_TWS_annual.csv`, `CSR_TWSA_Sen_2002_2024_common_land.nc`, `figure_statistics.json` |
| Extended Data Fig. 7 | `source_data/figED07` (`figED07_source`) | `source_data_panel_a.csv`, `source_data_panel_b.csv` |
| Extended Data Fig. 8 | `hydroclimate_trends` (`hydroclimate_trend_dir`), `fig/data` (`figure_statistics_dir`), `forcing_mask` | `GPCC_precip_sen_mk_1981_2024.nc`, `MSWEP_precip_sen_mk_1981_2024.nc`, `GPCP_precip_sen_mk_1981_2024.nc`, corresponding `*_sig_Dr_Wet_1981_2024.mat` files and `mask/Forcing_GL_0-359.75_90--90_0.25deg.mat` |
| Extended Data Fig. 9 | `source_data/figED09` (`figED09_source`) | `Africa_annual_SPEI_TWS.csv`, `Africa_trend_statistics.json`, `china_precipitation_original.csv`, `africa_TWSA_Sen_common.nc` |
| Supplementary Fig. 4 | `station_rh` and `station_wind` | Prepared GHCD station coordinate tables for RH (7,275 stations) and WS (5,515 stations) |

An Extended Data Fig. 7 configuration can be written as follows.

```json
{"paths": {"figED07_source": "/path/to/prepared/figED07"}}
```

Panel a reads 1981–2019 annual values and Sen trends for SSI, ERA5 SPEI, GHCD SPEI, and ERA5 SPEI with GHCD precipitation. Panel b reads the trend-gap reduction and paired bootstrap intervals. The precipitation replacement uses ERA5 atmospheric evaporative demand.

Each precipitation trend NetCDF contains `lon`, `lat`, `sen_slope` and `mk_p`. Each area-fraction MAT file contains `sig_Dr_Wet` in wetting, drying order. Station coordinates can be supplied as CSV or Excel tables with longitude and latitude columns, or as numeric MATLAB structures with longitude and latitude fields.

## Analysis commands

The `processing` directory includes forcing, atmospheric evaporative demand, SPEI, regional summaries, trend estimation and figure-table preparation.

```sh
python processing/hydroclimate_trends.py --base-dir /path/to/precipitation --output-dir hydroclimate_trends
python processing/analyze_precipitation_substitution.py --input-csv /path/to/global_series.csv --output-dir analysis_output
```

`hydroclimate_trends.py` calculates 1981–2024 precipitation trends from GPCC and MSWEP annual totals and GPCP monthly precipitation rates. Use `--help` for settings. The precipitation map command reads the resulting trend grids and the prepared area fractions and forcing mask.

`analyze_precipitation_substitution.py` reads prepared annual series with columns `year,SSI,GHCD,ERA5,JRA3Q,MIX`. `MIX` is SPEI calculated with GHCD precipitation and ERA5 atmospheric evaporative demand. The command calculates Sen slopes, compares absolute SSI–SPEI trend gaps, and generates paired bootstrap intervals.

## License

Project code is licensed under MIT. The KTAUB implementation is licensed under BSD 3-Clause as described in `KTAUB_LICENSE.txt`.

