# Appendix B — Supporting analyses

Supplementary material to *[paper reference]*. This appendix records three
analyses referred to in the main text but not reported there in full: the
evidence used to set and evaluate the capillary relation, the spatial structure
of the observed quantities, and the full statistics of the network transfer test.

---

## B.1 The capillary relation

The capillary flux is computed from the relation of Averyanov (Equation 2 of the
paper),

```
E = E0 (1 - H / H_cr)^n        for H < H_cr,   zero otherwise
```

where H is the depth to water and E0 the monthly FAO-56 reference
evapotranspiration. Neither parameter was fitted to the groundwater observations
of this study.

**Critical depth.** `H_cr = 3.0 m` is the value in regional use, adopted in
reclamation practice for these soils and used elsewhere in this paper as the
depth below which capillary transport to the root zone ceases.

**Exponent.** `n = 1.88` was set so that the relation reproduces the capillary
contribution measured in the field in the lower Amu Darya. Forkutsa et al. (2009)
report that shallow groundwater supplied up to 399 mm of the water consumed by
cotton at water tables of about 1.0–1.5 m in Khorezm. Requiring the relation to
return that flux at a water table of 1.5 m gives n = 1.87; the value adopted is
1.88.

### Evaluation against independent field data

Both published field studies available for the region measured capillary uptake
at water tables of 1.2–1.5 m, considerably shallower than the 2.38 m of the
calibration well. The relation was therefore evaluated at the depths at which
those measurements were made.

| Quantity | This relation | Measured | Source |
|---|---|---|---|
| Capillary flux at a water table of 1.2–1.5 m, mm/yr | 395–557 | up to 399 | Forkutsa et al. (2009) |
| Salt delivered to the root zone at 1.2–1.5 m, t ha⁻¹ yr⁻¹ | 12–15 | 3.5–14 | Ibrakhimov et al. (2007) |

The flux is reproduced by construction, since it was used to set the exponent.
The salt delivery is an independent check, because it involves the mineralization
of the groundwater as well as the flux, and there the relation lies at the upper
end of the measured range. It therefore does not underestimate capillary
transport.

### Sensitivity of the result to depth

The relation is strongly non-linear in depth, which is why the values at the
calibration well are much smaller than those measured at the shallower Khorezm
sites:

| Depth to water, m | Flux, mm/yr | Salt delivered at 2.62 g/L, t ha⁻¹ yr⁻¹ |
|---|---|---|
| 1.20 | 557 | 14.6 |
| 1.50 | 395 | 10.4 |
| 2.00 | 185 | 4.8 |
| 2.38 | 126 | 3.8 |
| 2.60 | 33 | 0.9 |

This non-linearity is the reason depth to water dominates mineralization as a
control on salt delivery (Section 5.1 of the paper).

### Uncertainty

Absolute salt delivery should be read as uncertain by about a factor of two,
which is the range over which alternative parameter pairs remain compatible with
the field measurements above. Section 4.5 of the paper shows that the ranking of
wells and the relative differences between configurations are not sensitive to
this.

---

## B.2 Spatial structure of the observed quantities

The model uses a single district-wide soil salt input, so simulated
mineralization is nearly uniform across the network by construction (2.58–3.59
g/L) while the observed three-year means span 1.34–8.09 g/L. The residual between
them is not explained by drain distance, drainage response time or depth to
water (Section 4.1). Its spatial structure is examined here.

### Empirical variograms

Semivariance by lag, computed on the 155 wells over the three-year means:

| Lag, km | Pairs | Salt delivery | Depth to water | Mineralization |
|---|---|---|---|---|
| 1–2 | 135 | 15.35 | 0.318 | 1.64 |
| 2–3 | 219 | 12.46 | 0.391 | 1.20 |
| 3–4 | 268 | 16.81 | 0.455 | 1.68 |
| 4–6 | 706 | 17.80 | 0.547 | 1.79 |
| 6–8 | 904 | 24.94 | 0.654 | 1.67 |
| 8–12 | 2,080 | 27.09 | 0.701 | 1.61 |
| 12–18 | 3,079 | 31.33 | 0.690 | 1.69 |
| **Sill (total variance)** | | **23.19** | **0.657** | **1.575** |

Depth to water shows a well-defined structure, rising from about half its sill at
the shortest resolvable lag to the sill at roughly 8 km. Mineralization is at its
sill already in the 1–2 km class: the network does not resolve its spatial
structure at all, and at the scale of the observations it behaves as pure nugget.
Salt delivery is intermediate, inheriting a range of the order of 8–12 km from
depth to water while carrying about two thirds of its variance at separations
below the shortest lag.

The mean distance between neighbouring wells is 1.4 km. The consequence is that
mineralization varies over distances the network cannot see, whereas depth to
water does not. This is why the mapped results are shown as symbols and
nearest-neighbour cells rather than as an interpolated surface, and why the
classification of Section 5.3 is described as applying to the neighbourhood of a
well rather than to an individual field.

### Moran's I

| Lag | Salt delivery | Depth to water | Mineralization |
|---|---|---|---|
| 1.5 km | 0.60 | 0.61 | 0.36 |
| 3 km | 0.44 | 0.48 | 0.18 |
| 5 km | 0.31 | 0.36 | 0.06 |
| 10 km | 0.12 | 0.15 | 0.01 |

Consistent with the variograms: autocorrelation in mineralization decays to
nothing within a few kilometres, while depth to water retains structure to about
10 km.

---

## B.3 Network transfer test — full statistics

The calibrated parameter set was applied to each of the 155 wells with only the
mapped distance to the nearest field drain changed; no parameter was adjusted to
the individual records. Ground surface elevations are not available for the
network, so the comparison is made on the seasonal shape of the depth-to-water
record with the mean of each series removed.

| Statistic | Median | Interquartile range |
|---|---|---|
| RMSE of the demeaned seasonal series, m | 0.36 | 0.30–0.45 |
| Correlation with the observed cycle | 0.76 | 0.60–0.82 |
| Observed seasonal amplitude, m | 1.63 | 1.31–1.96 |
| Simulated seasonal amplitude, m | 1.29 | 1.10–1.64 |
| Drainage response time implied by the mapped drain distance, days | 12.5 | — |

The correlation exceeds 0.5 at 81% of wells and 0.7 at 63%, and is negative at
6%. At the calibration well itself the RMSE is 0.20 m and the correlation 0.89.

Two systematic features follow from applying a single collector regime across the
district. Simulated seasonal amplitude falls short of the observed by 0.34 m at
the median. And agreement deteriorates where that regime is least representative:
with the distance to the drain and the corresponding response time (Spearman
ρ = +0.38 for both) and with the observed amplitude (ρ = +0.40); the dependence
on mean depth is weaker (ρ = +0.22).

The median residual in mineralization is −74 mg/L, so the model sits close to the
district median. As set out in B.2, that residual varies over distances shorter
than the network resolves.

### Per-well results

Per-well statistics are in [`network_transfer_test.csv`](network_transfer_test.csv):
well number, mapped drain distance, RMSE and correlation of the demeaned seasonal
series, observed and simulated amplitude, observed mean depth, simulated and
observed mineralization, and the implied response time.
