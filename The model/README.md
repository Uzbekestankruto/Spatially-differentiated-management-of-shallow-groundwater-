# Kegeyli mixing-cell model

`kegeyli_model.py` simulates the water and salt balance of the shallow aquifer
beneath an irrigated field in the Kegeyli district, Republic of Karakalpakstan,
Uzbekistan. It supports the analysis published in *[paper reference]*.

## What the model represents

The aquifer beneath one field is treated as a **single well-mixed volume** — a
mixing cell — rather than as a spatially resolved flow domain. The cell is the
uppermost part of the saturated zone, the interval that the district monitoring
wells actually sample. Nothing is transported *within* the cell; the model
tracks only what enters it, what leaves it, and how the composition of the
water inside changes as a result.

Five exchanges are represented:

1. **Surface input.** Leaching and irrigation water, and precipitation, arrive
   at the water table carrying salt dissolved from the soil profile.
2. **Capillary rise.** Water is drawn back up towards the surface whenever the
   water table is shallow enough, taking dissolved salt with it into the root
   zone.
3. **Return of that salt.** The salt deposited in the soil is washed back down
   by the next leaching campaign, so the exchange between aquifer and root zone
   is closed rather than one-way.
4. **Drainage.** Water is exchanged with the collector–drainage network in
   either direction, depending on whether the water table stands above or below
   the collector.
5. **Downward loss.** Part of the dissolved load sinks out of the exchanging
   zone into the stagnant deeper aquifer.

Geochemistry — mixing, speciation and cation exchange — is computed by IPhreeqc
at every step. The water and salt bookkeeping is done in Python.

The model advances in weekly steps, 52 weeks per year, and is run for 30 years.
**Those 30 years are numerical spin-up to a stable periodic solution, not a
reconstruction of the past.** The soil salt load is prescribed anew for each
annual cycle and the soil is not represented as a finite reservoir, so the run
says nothing about the long-term evolution of the soil salt inventory.

## Requirements and use

```bash
pip install IPhreeqcPy
```

A PHREEQC database is also needed; the published runs used `wateq_CH2O.dat`.
Set `DATABASE` and `OUTPUT_DIR` at the top of the script, then:

```python
import kegeyli_model as M

runner  = M.IPhreeqcRunner(M.DATABASE)
results = M.run_simulation(runner, M.TOTAL_CELLS)
M.print_calibration_report(results)
```

This reproduces the calibrated present state: mean mineralization 2617 mg/L,
April 3338, October 2261, mean depth to water 2.38 m, water-balance RMSE 0.21 m.

To run the seven configurations of the paper:

```python
for name, lu, dm, ssd, salt, ev, surf, coll in M.SCENARIOS:
    M.apply_scenario(lu, dm, ssd, salt, ev, surf, coll)
    results = M.run_simulation(runner, M.TOTAL_CELLS)
```

## The code, block by block

### Configuration constants (top of file)

Everything a user would normally change sits here: paths, the geometry of the
cell, the calibrated parameters, and the scenario multipliers that are reset by
`apply_scenario`.

`MIXING_THICKNESS` (5.5 m) is the thickness of the cell **in the reference
state**, taken from the screened interval of the observation wells. It is not
fitted. The state variable of the model is the mass of water in the cell, which
rises during leaching and irrigation, so the volume over which mixing is
computed is not constant.

`DOWNWARD_SALT_LOSS_PER_WEEK` (0.0114) is the **only fitted parameter** in the
model. Everything else is measured, taken from the literature, or derived from a
balance calculation.

### Unit conversion — `mm_to_mass_fraction`, `calculate_gwl`

The bridge between hydrology and geochemistry. PHREEQC works in kilograms of
water; the agricultural inputs are in millimetres of applied depth.
`mm_to_mass_fraction` converts one to the other through the cell area, thickness
and porosity.

`calculate_gwl` converts the mass of water in the cell to a water-table
elevation. It uses the **storage coefficient**, not the total porosity, because
only the drainable fraction of the pore water moves the water table, whereas the
whole water body takes part in mixing. This is why a 15% change in stored water
corresponds to a water-table movement of nearly 3 m.

### Calendar — `_build_week_month`, `WEEK_MONTH`, `get_cell_type`

The model year starts on 1 November and is divided into 52 weeks. Each week is
assigned to the calendar month containing its midpoint, so that monthly climate
data and monthly observations can be used without interpolation.

`get_cell_type` labels each week by what happens agriculturally: leaching in
February and March, irrigation in June and in August–September, and precipitation
alone in the remaining weeks. Irrigation is represented as two **periods** of
several weeks, not as discrete irrigation events.

### Climate and capillary rise — `ET0_MM_PER_DAY`, `capillary_flux_mm`

`ET0_MM_PER_DAY` holds FAO-56 Penman–Monteith reference evapotranspiration by
month, computed from district station data. This is the potential rate.

`capillary_flux_mm` attenuates it with depth to water using the relation of
Averyanov,

```
E = E0 (1 - H / H_cr)^n        for H < H_cr,   zero otherwise
```

with `H_CRITICAL_M = 3.0` and `AVERYANOV_N = 1.88`. The strong non-linearity of
this expression is what makes water-table depth, rather than mineralization, the
dominant control on how much salt reaches the root zone. The parameters were
fixed against an independent AquaCrop estimate, not calibrated against the
monitoring data.

### Soil salt balance — `leaching_percolate_mm_per_season`, `derive_leach_conc_scale`

Determines how concentrated the water arriving at the water table is. The salt
present in the leached layer is the product of its mean salt content
(`SOIL_SALT_KG_M3`, from the reported distribution of soil salinity classes in
the district) and its thickness (`LEACHED_DEPTH_M`, from the reclamation leaching
curve of Reeve and Hoffman). A fixed fraction of it dissolves
(`SALT_REMOVAL_EFF`), and the percolate concentration is that mass divided by the
volume of water passing through.

`LEACH_CONC_SCALE` is the resulting multiplier applied to the reference
composition of the leaching percolate. It is **derived, not fitted**, which is
what leaves the downward loss as the single adjustable parameter.

### Water input — `get_water_input`

Returns the fraction of the cell water mass arriving at the water table in a
given week. Precipitation is added in every week; campaign volumes are added on
top in the appropriate weeks.

Two multipliers act here and both remove water together with its dissolved salt,
so scaling the volume scales the salt load proportionally:

- `SURFACE_DIVERSION_MULT` — the **surface function** of the collector–drainage
  network, that is, how much of the applied leaching water the on-farm furrows
  and depressions intercept before it infiltrates;
- `SSD_INTERCEPT` — interception by shallow subsurface drainage.

### Drainage — `calculate_drainage_boussinesq`

Lateral exchange with the collector network, following the transient solution of
the linearized Boussinesq equation for drainage between parallel collectors
(Glover–Dumm):

```
f = 1 - exp(-dt / tau),      tau = Sy L^2 / (pi^2 Kd D)
```

The driving head is the difference between the simulated water table and the
seasonal water level in the collector. The formulation is **bidirectional**:
water leaves the cell while the water table stands above the collector and
returns when it falls below, which is the only lateral inflow in the model.
`DRAIN_MULT` scales the conductance of the drain connection and `COLLECTOR_SCALE`
the collector stage; together they represent the **subsurface function** of the
network.

### PHREEQC input builders — `build_leaching`, `build_irrigation_*`, `build_rain`

Each builder writes the PHREEQC input for one weekly step: a `SOLUTION` block for
the resident groundwater, another for the incoming water, a `MIX` block in the
proportions set by `get_water_input`, an `EXCHANGE` block for cation exchange, and
a `REACTION` block that removes the moles of H₂O lost to capillary rise.

`build_leaching` additionally releases the salt held in `SOIL_SALT_STORE`,
described below, into the percolate.

`BUILDERS` maps the week type produced by `get_cell_type` to the appropriate
builder.

### The capillary loop — `SOIL_SALT_STORE`, `update_solution_from_array`

This closes the exchange between aquifer and root zone, and is the reason the
model treats capillary rise differently from evaporation from an open surface.

Because capillary water moves upward as **liquid** and evaporates at or near the
soil surface, it carries its dissolved salt with it. `update_solution_from_array`
therefore removes from the cell not only the water lost to capillary rise but the
corresponding fraction of the solutes, and adds that salt to `SOIL_SALT_STORE`
(kg/m²). At the start of the next leaching campaign the store is released and
shared over the weeks of the campaign, raising the concentration of the percolate.

No salt is created or destroyed by the loop and it introduces no additional
parameter. Its net effect on the annual mean mineralization is negligible; what
it changes is the seasonal timing, and it makes the salt delivered to the root
zone a flux of the model rather than a quantity computed afterwards.

The same function applies `DOWNWARD_SALT_LOSS_PER_WEEK`, which removes a fixed
fraction of the remaining dissolved load — representing density-driven sinking of
concentrated water and leakage into the underlying unit together, since the
available data cannot separate them.

### Runner — `IPhreeqcRunner`, `run_simulation`

`IPhreeqcRunner` loads the database once and executes each step's input.
`run_simulation` walks the weeks, calls the appropriate builder, passes the
result to `update_solution_from_array`, caps the water table at the ground
surface and carries the surplus forward to the next step (water standing on a
flooded check), and collects the output.

### Diagnostics — `print_calibration_report` and the metric functions

`water_balance_rmse` compares the simulated depth to water with the observed
monthly means at the calibration well. `seasonal_tds` and `seasonal_amplitude`
report the April, July and October mineralization and their ratio, for comparison
with the monitoring network. `blending_potential` computes the share of the
post-flowering irrigation that groundwater could supply within the salinity limit
adopted for the local soils.

### Sensitivity sweeps — `run_salt_loss_sweep`, `run_mixing_sweep`, `run_soil_salinity_sweep`

Helpers for the sensitivity analysis. `run_salt_loss_sweep` is what was used to
calibrate the downward loss. `run_mixing_sweep` and `run_soil_salinity_sweep` test
the consequences of the adopted mixing thickness and soil salt content; neither is
a calibration, since both values come from outside the monitoring data.

Note that the water-balance RMSE is identical in every row of the mixing sweep:
the water balance does not depend on the mixing thickness, only the concentrations
do.

## Scenarios

`SCENARIOS` holds the seven configurations of the paper. Each is the calibrated
model with one or two multipliers displaced to a plausible bound; all other
parameters are held at their calibrated values.

| Configuration | Surface function | Collector stage | Drain conductance | SSD | Evaporation |
|---|---|---|---|---|---|
| Present state | 1.00 | 1.00 | 1.00 | — | 1.00 |
| Both functions maintained | 0.667 | 0.80 | 1.50 | — | 1.00 |
| Furrows maintained, drains silted | 0.667 | 1.20 | 0.25 | — | 1.00 |
| Furrows silted, drains maintained | 1.50 | 0.80 | 1.50 | — | 1.00 |
| Both functions degraded | 1.50 | 1.20 | 0.25 | — | 1.00 |
| Shallow subsurface drainage | 1.00 | 1.00 | 1.00 | 0.70 | 1.00 |
| Plastic film on the field | 1.00 | 1.00 | 1.00 | — | 0.10 |

Because the scenarios are counterfactual, the goodness-of-fit statistics of the
calibration do not apply to them; their value lies in the differences between
them.

## What the model does not do

- It does not simulate the unsaturated zone. The volume, timing and composition
  of the water arriving at the water table are prescribed from balance
  calculations, so preferential flow, root water uptake and reactions within the
  profile are absent.
- It does not resolve transport within the cell, nor vertical stratification of
  composition, which enters only as the one-way downward loss.
- It contains no crop, so the consequences of a deeper water table for crop water
  supply are outside its scope.
- It uses a single district-wide soil salt input, so simulated mineralization is
  nearly uniform across wells by construction.
