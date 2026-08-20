# Appendix B — Supporting analyses

Supplementary material to *[paper reference]*. This appendix records three
analyses referred to in the main text but not reported there in full: the
external estimate used to constrain the capillary parameters, the screening of
spatial predictors of groundwater mineralization, and the full statistics of the
network transfer test.

---

## B.1 The AquaCrop estimate used to constrain the capillary relation

The critical depth `H_cr` and the exponent `n` of the capillary relation
(Equation 2 of the paper) are not identifiable from the groundwater observations
alone, because several pairs reproduce the same flux over the observed range of
water-table depths. They were therefore constrained against an estimate obtained
independently of the groundwater model.

*[To be completed by the authors: the AquaCrop configuration used — crop, soil
profile, irrigation schedule, groundwater boundary condition, and the period
simulated — together with the water and salt balance from which the figure of
2.89 t ha⁻¹ yr⁻¹ of secondary salinization was derived.]*

Reproducing that salt delivery on the observed depth-to-water record of the
calibration well requires an annual capillary flux of about 100 mm. Several
parameter pairs satisfy this constraint. The pair adopted, `H_cr = 3.0 m` and
`n = 1.88`, was chosen because 3.0 m is the critical depth used elsewhere in the
paper and in regional practice.

The consequence of that choice is limited. Across the alternative pairs
consistent with the same constraint, the absolute capillary salt delivery
computed for the monitoring network changes by up to about a fifth, whereas the
ranking of wells is almost unaffected: rank correlations between the resulting
well orderings lie between 0.994 and 0.999. The mapped pattern and the relative
differences between scenarios are therefore insensitive to the choice, while the
absolute values are not.

---

## B.2 Screening of spatial predictors of groundwater mineralization

The model uses a single district-wide soil salt input, so simulated
mineralization is nearly uniform across the network by construction (2.58–3.59
g/L) while the observed three-year means span 1.34–8.09 g/L. Section 4.1 of the
paper states that the residual variation is not explained by any geometric or
land-cover predictor available. The screening behind that statement is recorded
here.

*[To be completed by the authors: the predictors tested — distance to the nearest
drain, distance to the nearest canal, position within the inter-drain span, land
cover fractions within a buffer, distance to bare land, elevation — with the
correlation or regression statistics obtained for each.]*

Two results from that screening are used in the main text. Neither the distance
to the nearest drain, the drainage response time, nor the depth to water explains
the residual between simulated and observed mineralization. And the spatial
structure of mineralization is short-ranged: its variogram range is of the order
of 2 km, comparable with or shorter than the mean spacing between neighbouring
wells of 1.4 km, whereas depth to water is autocorrelated over some 8 km. This is
the basis for the statement that the network resolves the neighbourhood of a well
rather than an individual field.

### Spatial structure of the mapped quantity

Moran's I for the modelled capillary salt delivery, computed on the observed
records of the 155 wells:

| Lag | Capillary salt delivery | Depth to water | Mineralization |
|---|---|---|---|
| 1.5 km | 0.60 | 0.61 | 0.36 |
| 3 km | 0.44 | 0.48 | 0.18 |
| 5 km | 0.31 | 0.36 | 0.06 |
| 10 km | 0.12 | 0.15 | 0.01 |

The mapped quantity is spatially structured, and inherits that structure from
depth to water rather than from mineralization. About half of its variance
nevertheless lies at separations shorter than the mean well spacing, which is why
the results are shown as symbols and nearest-neighbour cells rather than as an
interpolated surface.

---

## B.3 Network transfer test — full statistics

The calibrated parameter set was applied to each of the 155 wells with only the
mapped distance to the nearest field drain changed; no parameter was adjusted to
the individual records. Ground surface elevations are not available for the
network, so the comparison is made on the seasonal shape of the depth-to-water
record with the mean of each series removed.

| Statistic | Median | Interquartile range |
|---|---|---|
| Root mean square error of the seasonal shape, m | 0.36 | 0.30–0.45 |
| Correlation with the observed cycle | 0.76 | 0.60–0.82 |
| Observed seasonal amplitude, m | 1.63 | 1.31–1.96 |
| Simulated seasonal amplitude, m | 1.29 | 1.10–1.64 |
| Drainage response time implied by the mapped drain distance, days | 12.5 | — |

The correlation exceeds 0.5 at 81% of wells and 0.7 at 63%, and is negative at
6%. At the calibration well itself the error is 0.20 m and the correlation 0.89.

Two systematic features follow from applying a single collector regime across the
district. Simulated seasonal amplitude falls short of the observed by 0.34 m at
the median. And agreement deteriorates where that regime is least representative:
with the distance to the drain and the corresponding response time (Spearman
ρ = +0.38 for both) and with the observed amplitude (ρ = +0.40); the dependence
on mean depth is weaker (ρ = +0.22).

The median residual in mineralization is −74 mg/L, so the model sits close to the
district median. As set out in B.2, that residual is not explained by the
predictors tested.

### Per-well results

Per-well statistics are in [`network_transfer_test.csv`](network_transfer_test.csv):
well number, mapped drain distance, error and correlation of the seasonal shape,
observed and simulated amplitude, observed mean depth, simulated and observed
mineralization, and the implied response time.
