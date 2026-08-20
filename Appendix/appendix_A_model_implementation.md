# Appendix A — Model implementation

Supplementary material to *[paper reference]*. This appendix records the weekly
forcing of the model and the parameter values behind it, so that the runs can be
reproduced without reading the source code. The model itself is in
[`../model/kegeyli_model.py`](../model/kegeyli_model.py).

## A.1 Weekly water inputs and boundary conditions

The model year begins on 1 November and runs for 52 weeks; each week is assigned
to the calendar month containing its midpoint. The columns below give the depths
**arriving at the water table**, not the depths applied at the surface: the
applied leaching norm of 5,175 m³/ha (517.6 mm) is reduced by the factor
`LEACHING_SCALE = 0.5` and the applied irrigation norm of 3,490 m³/ha by
`IRRIGATION_SCALE = 0.7`, the remainder being retained in the profile or
evaporated after infiltration.

Collector stage is given as metres above the base elevation of 67 m a.s.l.

| Week | Month | Week type | Leaching | Irrigation 1 | Irrigation 2 | Precipitation | ET₀, mm/d | Collector stage, m |
|---:|:--|:--|---:|---:|---:|---:|---:|---:|
| 1 | Nov | rain | 0.0 | 0.0 | 0.0 | 0.68 | 1.62 | 2.40 |
| 2 | Nov | rain | 0.0 | 0.0 | 0.0 | 0.68 | 1.62 | 2.40 |
| 3 | Nov | rain | 0.0 | 0.0 | 0.0 | 0.68 | 1.62 | 2.40 |
| 4 | Nov | rain | 0.0 | 0.0 | 0.0 | 0.68 | 1.62 | 2.40 |
| 5 | Dec | rain | 0.0 | 0.0 | 0.0 | 2.54 | 0.71 | 2.33 |
| 6 | Dec | rain | 0.0 | 0.0 | 0.0 | 2.54 | 0.71 | 2.33 |
| 7 | Dec | rain | 0.0 | 0.0 | 0.0 | 2.54 | 0.71 | 2.33 |
| 8 | Dec | rain | 0.0 | 0.0 | 0.0 | 2.54 | 0.71 | 2.33 |
| 9 | Dec | rain | 0.0 | 0.0 | 0.0 | 2.54 | 0.71 | 2.33 |
| 10 | Jan | rain | 0.0 | 0.0 | 0.0 | 0.17 | 0.59 | 2.20 |
| 11 | Jan | rain | 0.0 | 0.0 | 0.0 | 0.17 | 0.59 | 2.20 |
| 12 | Jan | rain | 0.0 | 0.0 | 0.0 | 0.17 | 0.59 | 2.20 |
| 13 | Jan | rain | 0.0 | 0.0 | 0.0 | 0.17 | 0.59 | 2.20 |
| 14 | Feb | leaching | 69.3 | 0.0 | 0.0 | 2.33 | 1.14 | 2.90 |
| 15 | Feb | leaching | 33.9 | 0.0 | 0.0 | 2.33 | 1.14 | 2.90 |
| 16 | Feb | leaching | 16.9 | 0.0 | 0.0 | 2.33 | 1.14 | 2.90 |
| 17 | Feb | leaching | 9.2 | 0.0 | 0.0 | 2.33 | 1.14 | 2.90 |
| 18 | Mar | leaching | 64.7 | 0.0 | 0.0 | 1.02 | 3.28 | 3.03 |
| 19 | Mar | leaching | 31.6 | 0.0 | 0.0 | 1.02 | 3.28 | 3.03 |
| 20 | Mar | leaching | 15.8 | 0.0 | 0.0 | 1.02 | 3.28 | 3.03 |
| 21 | Mar | leaching | 8.6 | 0.0 | 0.0 | 1.02 | 3.28 | 3.03 |
| 22 | Mar | leaching | 8.6 | 0.0 | 0.0 | 1.02 | 3.28 | 3.03 |
| 23 | Apr | rain | 0.0 | 0.0 | 0.0 | 1.98 | 5.09 | 3.30 |
| 24 | Apr | rain | 0.0 | 0.0 | 0.0 | 1.98 | 5.09 | 3.30 |
| 25 | Apr | rain | 0.0 | 0.0 | 0.0 | 1.98 | 5.09 | 3.30 |
| 26 | Apr | rain | 0.0 | 0.0 | 0.0 | 1.98 | 5.09 | 3.30 |
| 27 | May | rain | 0.0 | 0.0 | 0.0 | 3.02 | 6.32 | 3.20 |
| 28 | May | rain | 0.0 | 0.0 | 0.0 | 3.02 | 6.32 | 3.20 |
| 29 | May | rain | 0.0 | 0.0 | 0.0 | 3.02 | 6.32 | 3.20 |
| 30 | May | rain | 0.0 | 0.0 | 0.0 | 3.02 | 6.32 | 3.20 |
| 31 | Jun | irrigation_1 | 0.0 | 36.2 | 0.0 | 0.72 | 8.18 | 2.88 |
| 32 | Jun | irrigation_1 | 0.0 | 36.2 | 0.0 | 0.72 | 8.18 | 2.88 |
| 33 | Jun | irrigation_1 | 0.0 | 36.2 | 0.0 | 0.72 | 8.18 | 2.88 |
| 34 | Jun | irrigation_1 | 0.0 | 36.2 | 0.0 | 0.72 | 8.18 | 2.88 |
| 35 | Jun | irrigation_1 | 0.0 | 36.2 | 0.0 | 0.72 | 8.18 | 2.88 |
| 36 | Jul | rain | 0.0 | 0.0 | 0.0 | 0.57 | 7.83 | 2.80 |
| 37 | Jul | rain | 0.0 | 0.0 | 0.0 | 0.57 | 7.83 | 2.80 |
| 38 | Jul | rain | 0.0 | 0.0 | 0.0 | 0.57 | 7.83 | 2.80 |
| 39 | Jul | rain | 0.0 | 0.0 | 0.0 | 0.57 | 7.83 | 2.80 |
| 40 | Aug | irrigation_2 | 0.0 | 0.0 | 7.0 | 1.12 | 6.57 | 3.00 |
| 41 | Aug | irrigation_2 | 0.0 | 0.0 | 7.0 | 1.12 | 6.57 | 3.00 |
| 42 | Aug | irrigation_2 | 0.0 | 0.0 | 7.0 | 1.12 | 6.57 | 3.00 |
| 43 | Aug | irrigation_2 | 0.0 | 0.0 | 7.0 | 1.12 | 6.57 | 3.00 |
| 44 | Sep | irrigation_2 | 0.0 | 0.0 | 7.0 | 0.82 | 4.17 | 2.75 |
| 45 | Sep | irrigation_2 | 0.0 | 0.0 | 7.0 | 0.82 | 4.17 | 2.75 |
| 46 | Sep | irrigation_2 | 0.0 | 0.0 | 7.0 | 0.82 | 4.17 | 2.75 |
| 47 | Sep | irrigation_2 | 0.0 | 0.0 | 7.0 | 0.82 | 4.17 | 2.75 |
| 48 | Sep | irrigation_2 | 0.0 | 0.0 | 7.0 | 0.82 | 4.17 | 2.75 |
| 49 | Oct | rain | 0.0 | 0.0 | 0.0 | 2.65 | 2.18 | 2.46 |
| 50 | Oct | rain | 0.0 | 0.0 | 0.0 | 2.65 | 2.18 | 2.46 |
| 51 | Oct | rain | 0.0 | 0.0 | 0.0 | 2.65 | 2.18 | 2.46 |
| 52 | Oct | rain | 0.0 | 0.0 | 0.0 | 2.65 | 2.18 | 2.46 |
| **Total** | | | **259** | **181** | **63** | **75.6** | | |

Seasonal totals reaching the water table: leaching 259 mm, first irrigation
181 mm, second irrigation 63 mm, precipitation 76 mm.

## A.2 Parameter values

### Fixed from measurement or from the literature

| Parameter | Symbol in code | Value | Source |
|---|---|---|---|
| Thickness of the exchanging zone, reference state | `MIXING_THICKNESS` | 5.5 m | screened interval of the observation wells |
| Total porosity | `POROSITY` | 0.28 | regional hydrogeological data |
| Storage coefficient | `STORAGE_COEFF` | 0.078 | water-balance reconstruction |
| Critical depth for capillary rise | `H_CRITICAL_M` | 3.0 m | regional practice |
| Exponent of the capillary relation | `AVERYANOV_N` | 1.88 | constrained against the AquaCrop estimate (Appendix B) |
| Soil salt content of the leached layer | `SOIL_SALT_KG_M3` | 8.6 kg/m³ | district soil salinity classes |
| Thickness of the leached layer | `LEACHED_DEPTH_M` | 0.69 m | reclamation leaching curve |
| Leaching removal efficiency | `SALT_REMOVAL_EFF` | 0.60 | district measurement |
| Distance to the nearest drain | `FIXED_DRAIN_DISTANCE` | 122 m | digitized drainage network |

### Reconstructed in the water-balance step

| Parameter | Symbol in code | Value |
|---|---|---|
| Leaching water reaching the water table | `LEACHING_SCALE` | 0.5 |
| Irrigation water reaching the water table | `IRRIGATION_SCALE` | 0.7 |
| Seasonal collector stage | `COLLECTOR_LEVELS_BY_MONTH` | see the table above |
| Effective conductivity of the drain connection | `DRAIN_EFFICIENCY` | 0.066 |

The collector curve carries one free value per month against one observed depth
per month, so the agreement it produces is a reconstruction and not an
independent test. This is discussed in Section 3.3 of the paper.

### Fitted

| Parameter | Symbol in code | Value |
|---|---|---|
| Downward loss of dissolved salt | `DOWNWARD_SALT_LOSS_PER_WEEK` | 0.0114 per week, 45% of the dissolved load per year |

This is the only parameter adjusted against the observed mineralization.

## A.3 Reproducing the published runs

```python
import kegeyli_model as M

runner  = M.IPhreeqcRunner(M.DATABASE)
results = M.run_simulation(runner, M.TOTAL_CELLS)
M.print_calibration_report(results)
```

Expected output for the calibrated present state: mean mineralization 2617 mg/L,
April 3338 mg/L, October 2261 mg/L, mean depth to water 2.38 m, water-balance
RMSE 0.21 m.

The seven configurations of the paper are held in `M.SCENARIOS`; see the model
README for the loop.
