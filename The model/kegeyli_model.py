"""
Lumped mixing-cell model of the shallow aquifer beneath an irrigated field,
Kegeyli district, Republic of Karakalpakstan, Uzbekistan.

The cell represents the uppermost part of the saturated zone, the interval the
monitoring wells sample. It is treated as fully mixed at each weekly step. Water
and dissolved salt enter from leaching, irrigation and precipitation; water and
salt leave upward by capillary rise to the root zone and return with the next
leaching campaign; water is exchanged with the collector-drainage network in
either direction; and part of the remaining dissolved load is lost downward to
the stagnant part of the aquifer. Geochemistry is computed with IPhreeqc.

The state variable is the mass of water in the cell, so the volume over which
mixing is computed is not constant: MIXING_THICKNESS is the thickness of the
reference state. Changes in mass are converted to changes in water-table
elevation through the storage coefficient, not the total porosity, since only
the drainable fraction moves the water table.

Calibration (see the paper, Section 3.3): the salt input is fixed from
independent data and the downward salt loss is the single fitted parameter.

The 30 repeated years are numerical spin-up to a stable periodic solution, not
a reconstruction of the past. The district-average leaching salt load is
prescribed anew for every annual cycle; the soil is not represented as a finite
reservoir, so the run says nothing about the long-term evolution of the soil
salt inventory.

Requires IPhreeqcPy and a PHREEQC database (wateq_CH2O.dat).
"""


import os
import math

# ============================================================
# PATHS
# ============================================================
DATABASE = r"/mnt/user-data/uploads/wateq_CH2O.dat"
OUTPUT_DIR = r"/home/claude/out"

# ============================================================
# >>> CALIBRATION KNOBS  (this is the ONLY block I will change) <<<
# ============================================================
FIXED_DRAIN_DISTANCE = 122.0   # distance well 45 -> nearest drain (m)
LEACHING_SCALE   = 0.5         # multiplies leaching water input (lower = less leaching water)
IRRIGATION_SCALE = 0.7         # multiplies irrigation water input
EVAP_SCALE       = 1.0         # multiplies evaporation
COLLECTOR_SCALE  = 1.0         # multiplies seasonal collector levels (lower = less backwater = more drainage)
DRAIN_EFFICIENCY = 0.066         # 0..1, effective K reduction (higher = stronger drainage)
SURFACE_ELEVATION = 72.3       # model ground surface used to convert GWL -> depth-to-water

# ------------------------------------------------------------
# SALT INPUT FROM THE SOIL PROFILE  (replaces the old assumption that leaching
# percolate has the composition of collector-drainage water)
#
# The leaching percolate no longer carries an assumed TDS. Its concentration is
# derived from a salt balance over the leached layer:
#
#   salt present per m2      M   = SOIL_SALT_KG_M3 * LEACHED_DEPTH_M      [kg/m2]
#   salt dissolved by leaching    = M * SALT_REMOVAL_EFF                  [kg/m2]
#   water passing through the profile per season = V                      [m3/m2]
#   percolate concentration  C   = M * SALT_REMOVAL_EFF / V               [kg/m3 = g/L]
#
# SOIL_SALT_KG_M3 is a single representative mean for the whole district, held
# constant and NOT made a function of water-table depth. Two effects oppose each
# other there - a deeper water table leaves a thicker unsaturated column (more
# salt), while a shallower water table drives stronger capillary rise (more
# intense secondary salinization) - and neither is constrained by the available
# data. Holding it constant makes the per-well residual interpretable: where the
# model overpredicts mineralization the local soil is less saline than the
# district mean, and where it underpredicts it is more saline.
# ------------------------------------------------------------
SOIL_SALT_KG_M3   = 8.6        # representative mean salt density of the leached layer [kg/m3]
LEACHED_DEPTH_M   = 0.69       # thickness of the layer flushed by leaching in checks [m]
SALT_REMOVAL_EFF  = 0.60       # fraction of the salt present that dissolves into leaching water

# Reference composition whose TDS is rescaled to the value derived above.
# Ion proportions are kept from the measured collector-drainage water.
REFERENCE_LEACH_TDS_MGL = 3563.0
LEACH_CONC_SCALE = 1.0         # recomputed from the salt balance below; do not set by hand

# ============================================================
# SCENARIO PARAMETERS (set per scenario by the runner below)
# ============================================================
LAND_USE = "field"                 # only irrigated fields are represented
DRAIN_MULT = 1.0                   # multiplier on the calibrated drainage response
SSD_INTERCEPT = 0.0                # fraction of leaching percolate intercepted by shallow drains

# Surface diversion of leaching water by the on-farm depressions/furrows that also
# serve as summer water-delivery channels. This is the SURFACE function of the
# drainage network, separate from the subsurface function scaled by DRAIN_MULT.
# It multiplies the leaching water (and hence its salt load) that actually reaches
# the water table:
#   1.000 = calibrated present state (whatever diversion exists is embedded in
#           LEACHING_SCALE = 0.5 and must not be double-counted)
#   0.667 = better-maintained surface network; matches the input ratio implied by
#           the Tier-1 "60% efficiency" scenario relative to its "40%" reference
#   1.500 = silted / degraded surface network; matches the Tier-1 "10%" scenario
SURFACE_DIVERSION_MULT = 1.0


# Groundwater blending (report section 5.7): irrigation events 2 and 3
SURFACE_WATER_TDS_GL = 1.3         # Akmalov (2024)
TARGET_BLEND_TDS_GL  = 2.5         # Chembarisov et al. (2022), meadow-takyr soils
POST_FLOWERING_MM    = 91          # the model's second irrigation period; 1 mm = 10 m3/ha
# ============================================================


YEARS = 30
WEEKS_PER_YEAR = 52
TOTAL_CELLS = YEARS * WEEKS_PER_YEAR

BR_TO_HCO3 = 61.0171 / 79.904
HCO3_TO_BR = 79.904 / 61.0171

# ============================================================
# AQUIFER PARAMETERS
# ============================================================
CELL_LENGTH = 7.0
CELL_WIDTH = 1.0
CELL_AREA = CELL_LENGTH * CELL_WIDTH
SATURATED_THICKNESS = 22          # full saturated thickness of the aquifer [m]
POROSITY = 0.28

# ------------------------------------------------------------
# ACTIVE MIXING ZONE
# The aquifer is vertically stratified in composition and the observation wells
# sample its upper part, so the cell represents the shallow zone that takes part
# in seasonal exchange, not the full saturated thickness. This thickness sets the
# CHEMICAL buffer volume only. The water-table response is computed from the
# volume of water added or removed divided by the storage coefficient, so it is
# INDEPENDENT of MIXING_THICKNESS: changing it changes concentrations, never the
# water balance.
# ------------------------------------------------------------
MIXING_THICKNESS = 5.5            # active mixing zone [m]

# ------------------------------------------------------------
# DOWNWARD SALT LOSS FROM THE ACTIVE ZONE
# Groundwater in this aquifer is vertically stratified: 0.6-10 g/L near the
# surface against 15-40 g/L at depth (Akmalov, 2024), and archival cross-sections
# show near-vertical contacts between zones of contrasting mineralization. Water
# that has been concentrated by evaporation and by dissolution of soil salt is
# denser than the ambient groundwater and sinks out of the actively exchanging
# upper zone into the stagnant lower part of the same Quaternary aquifer.
# A lumped cell cannot resolve free convection, so the process enters as a
# first-order loss: a fixed fraction of the dissolved load leaves the cell each
# week, at the cell's own composition and without removing water. Calibrate this
# rate against the observed mean mineralization on a run long enough to reach
# equilibrium (see YEARS below).
# ------------------------------------------------------------
DOWNWARD_SALT_LOSS_PER_WEEK = 0.0114  # calibrated value (45% of the dissolved load per year)

# Storage coefficient used to convert a change in stored water to a change in
# water-table elevation. 0.078 is the effective value implied by the existing
# calibration; the literature specific yield for these sediments is 0.15
# (Johnson, 1967) and is retained in the drainage term below.
STORAGE_COEFF = 0.078

WATER_IN_CELL = CELL_AREA * MIXING_THICKNESS * POROSITY * 1000   # litres in the mixing zone
BASE_MASS = 1.0
BASE_GWL = 67.0
# metres of water-table rise per unit of relative water mass
GWL_FACTOR = MIXING_THICKNESS * POROSITY / STORAGE_COEFF


def mm_to_mass_fraction(mm):
    """Depth of water [mm] arriving at the water table -> fraction of the mixing-zone volume."""
    return mm / 1000.0 / (MIXING_THICKNESS * POROSITY)

MAX_GWL = SURFACE_ELEVATION
MAX_MASS = BASE_MASS + (MAX_GWL - BASE_GWL) / GWL_FACTOR
MIN_MASS = 0.5

K_HORIZONTAL_BASE = 13.0
K_HORIZONTAL = K_HORIZONTAL_BASE * DRAIN_EFFICIENCY

DRAIN_SPACING = 1000
MAX_DRAIN_DISTANCE = DRAIN_SPACING / 2
CELL_TRAVEL_DISTANCE = 7.0
EFFECTIVE_THICKNESS = 22.0
TIME_STEP = 7.0
SPECIFIC_YIELD = 0.15




def calculate_drain_distance(cell_index):
    if FIXED_DRAIN_DISTANCE is not None:
        return float(FIXED_DRAIN_DISTANCE)
    distance_from_start = cell_index * CELL_TRAVEL_DISTANCE
    distance_in_field = distance_from_start % DRAIN_SPACING
    if distance_in_field == 0:
        distance_in_field = DRAIN_SPACING
    distance_to_nearest = min(distance_in_field, DRAIN_SPACING - distance_in_field)
    distance_to_nearest = max(distance_to_nearest, CELL_TRAVEL_DISTANCE)
    distance_to_nearest = min(distance_to_nearest, MAX_DRAIN_DISTANCE)
    return distance_to_nearest


# ---------------------------------------------------------------------------
# Evaporation from the water table
# ---------------------------------------------------------------------------
# Potential rate: FAO-56 Penman-Monteith reference evapotranspiration computed
# from district station data (mm/day, by calendar month).
ET0_MM_PER_DAY = {1: 0.59, 2: 1.14, 3: 3.28, 4: 5.09, 5: 6.32, 6: 8.18,
                  7: 7.83, 8: 6.57, 9: 4.17, 10: 2.18, 11: 1.62, 12: 0.71}

# Attenuation with depth to water, after Averyanov:  E = E0 (1 - H/H_cr)^n
H_CRITICAL_M = 3.0
AVERYANOV_N  = 1.88

# Monthly precipitation, mm (district station).
PRECIP_MM_PER_MONTH = {1: 0.7, 2: 9.3, 3: 5.1, 4: 7.9, 5: 12.1, 6: 3.6,
                       7: 2.3, 8: 4.5, 9: 4.1, 10: 10.6, 11: 2.7, 12: 12.7}

# The model year starts on 1 November; each week is assigned to the calendar
# month containing its midpoint.
_MONTH_ORDER = [11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
_MONTH_LEN = {11: 30, 12: 31, 1: 31, 2: 28, 3: 31, 4: 30,
              5: 31, 6: 30, 7: 31, 8: 31, 9: 30, 10: 31}

def _build_week_month():
    bounds, day = [], 0
    for m in _MONTH_ORDER:
        bounds.append((m, day + 1, day + _MONTH_LEN[m]))
        day += _MONTH_LEN[m]
    out = {}
    for wk in range(1, WEEKS_PER_YEAR + 1):
        mid = (wk - 1) * 7 + 4
        for m, a, b in bounds:
            if a <= mid <= b:
                out[wk] = m
                break
    return out

WEEK_MONTH = _build_week_month()
WEEKS_IN_MONTH = {m: sum(1 for w, mm in WEEK_MONTH.items() if mm == m)
                  for m in _MONTH_ORDER}
EVAPORATION_MM_PER_DAY = {wk: ET0_MM_PER_DAY[WEEK_MONTH[wk]]
                          for wk in range(1, WEEKS_PER_YEAR + 1)}
PRECIP_MM_PER_WEEK = {wk: PRECIP_MM_PER_MONTH[WEEK_MONTH[wk]] / WEEKS_IN_MONTH[WEEK_MONTH[wk]]
                      for wk in range(1, WEEKS_PER_YEAR + 1)}

_obs_month_dtw = {  # month -> observed DTW (m)
    "Nov": 2.87, "Dec": 2.91, "Jan": 3.13, "Feb": 2.45, "Mar": 1.50,
    "Apr": 1.97, "May": 2.04, "Jun": 2.38, "Jul": 2.51, "Aug": 2.28,
    "Sep": 2.55, "Oct": 2.87,
}
_MONTH_NAME = {11: "Nov", 12: "Dec", 1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr",
               5: "May", 6: "Jun", 7: "Jul", 8: "Aug", 9: "Sep", 10: "Oct"}
_week_month = {wk: _MONTH_NAME[WEEK_MONTH[wk]] for wk in range(1, WEEKS_PER_YEAR + 1)}

# Composition of rainfall reaching the water table (mg/L): dilute.
RAIN_COMPOSITION = (0.6, 1.5, 1.2, 0.3, 3.0, 1.5, 0.5, 0.1, 6.5)

# Observed depth-to-water at well 45 (m), mapped to the weeks of the model year.
# Source: 2023-2025 monitoring, 3x/month, averaged by month, mapped to model weeks
# (weeks are assigned to the calendar month containing their midpoint)
OBSERVED_DTW_WELL45 = {wk: _obs_month_dtw[_week_month[wk]]
                       for wk in range(1, WEEKS_PER_YEAR + 1)}


def capillary_flux_mm(week_in_year, mass_h2o):
    """Depth-attenuated flux from the water table for one weekly step (mm)."""
    depth = SURFACE_ELEVATION - calculate_gwl(mass_h2o)
    if depth >= H_CRITICAL_M:
        return 0.0
    reduction = (1.0 - max(depth, 0.0) / H_CRITICAL_M) ** AVERYANOV_N
    return EVAPORATION_MM_PER_DAY[week_in_year] * EVAP_SCALE * reduction * TIME_STEP

def calculate_evaporation_moles(week_in_year, mass_h2o=BASE_MASS):
    """Moles of H2O leaving the cell by capillary rise in one step."""
    return mm_to_mass_fraction(capillary_flux_mm(week_in_year, mass_h2o)) * 1000.0 / 18.015

COLLECTOR_LEVELS_BY_MONTH = {11: 2.40, 12: 2.33, 1: 2.20, 2: 2.90, 3: 3.03,
                             4: 3.30, 5: 3.20, 6: 2.88, 7: 2.80, 8: 3.00,
                             9: 2.75, 10: 2.46}
COLLECTOR_LEVELS_SEASONAL = {wk: COLLECTOR_LEVELS_BY_MONTH[WEEK_MONTH[wk]]
                             for wk in range(1, WEEKS_PER_YEAR + 1)}


def get_collector_level(week_in_year):
    return COLLECTOR_LEVELS_SEASONAL.get(week_in_year, 0.5) * COLLECTOR_SCALE


# Water arriving at the water table, in mm per week (was: fractions of the cell
# volume). Totals: leaching 517 mm/season (5,174 m3/ha), irrigation 349 mm/season
# (3,493 m3/ha) -- identical to the calibrated version.
# Applied depths (mm) arriving at the water table, spread over the weeks of each
# campaign. Leaching runs in February-March, irrigation in June and in
# August-September: two periods of several weeks, not two discrete events.
_LEACH_SHAPE = [138.6, 67.8, 33.9, 18.5]

def _weeks_of(*months):
    return [w for w in range(1, WEEKS_PER_YEAR + 1) if WEEK_MONTH[w] in months]

def _spread_shape(weeks, shape):
    vals = [shape[min(i, len(shape) - 1)] for i in range(len(weeks))]
    scale = sum(shape) / sum(vals)
    return {w: v * scale for w, v in zip(weeks, vals)}

# The shape is applied within each month of the campaign, so the seasonal total
# is twice the sum of the shape.
LEACHING_WATER_INPUT = {}
for _m in (2, 3):
    LEACHING_WATER_INPUT.update(_spread_shape(_weeks_of(_m), _LEACH_SHAPE))
IRRIGATION_1_WATER_INPUT = {w: 258.72 / len(_weeks_of(6)) for w in _weeks_of(6)}
IRRIGATION_2_WATER_INPUT = {w: 90.58 / len(_weeks_of(8, 9)) for w in _weeks_of(8, 9)}


def leaching_percolate_mm_per_season():
    return sum(LEACHING_WATER_INPUT.values()) * LEACHING_SCALE


def derive_leach_conc_scale(verbose=True):
    global LEACH_CONC_SCALE
    salt_present = SOIL_SALT_KG_M3 * LEACHED_DEPTH_M           # kg/m2
    salt_dissolved = salt_present * SALT_REMOVAL_EFF           # kg/m2
    water_m = leaching_percolate_mm_per_season() / 1000.0      # m3/m2
    conc_gL = salt_dissolved / water_m                         # kg/m3 = g/L
    LEACH_CONC_SCALE = conc_gL * 1000.0 / REFERENCE_LEACH_TDS_MGL
    if verbose:
        print("Leaching percolate derived from the soil salt balance:")
        print(f"  soil salt {SOIL_SALT_KG_M3:.1f} kg/m3 over {LEACHED_DEPTH_M:.2f} m "
              f"-> {salt_present:.2f} kg/m2 present")
        print(f"  removal efficiency {SALT_REMOVAL_EFF:.0%} "
              f"-> {salt_dissolved:.2f} kg/m2 dissolved")
        print(f"  percolate volume {water_m*1000:.0f} mm/season "
              f"({water_m*10000:.0f} m3/ha)")
        print(f"  percolate TDS = {conc_gL:.2f} g/L "
              f"(LEACH_CONC_SCALE = {LEACH_CONC_SCALE:.3f})\n")
    return LEACH_CONC_SCALE


derive_leach_conc_scale(verbose=False)


def get_water_input(cell_type, week_in_year):
    """Fraction of the cell water mass arriving at the water table this step.

    Precipitation is added in every week; the campaign volumes are added on top.
    SURFACE_DIVERSION_MULT is the surface interception by on-farm furrows and
    SSD_INTERCEPT the shallow subsurface interception. Both remove water together
    with its dissolved salt, so scaling the volume scales the salt proportionally.
    """
    mm = PRECIP_MM_PER_WEEK[week_in_year]
    if cell_type == "leaching":
        mm += (LEACHING_WATER_INPUT.get(week_in_year, 0.0) * LEACHING_SCALE
               * SURFACE_DIVERSION_MULT * (1.0 - SSD_INTERCEPT))
    elif cell_type == "irrigation_1":
        mm += IRRIGATION_1_WATER_INPUT.get(week_in_year, 0.0) * IRRIGATION_SCALE
    elif cell_type == "irrigation_2":
        mm += IRRIGATION_2_WATER_INPUT.get(week_in_year, 0.0) * IRRIGATION_SCALE
    return mm_to_mass_fraction(mm)



initial_solution = {
    "Mg": 17, "Na": 90, "Ca": 73, "K": 5.00, "S(6)": 166, "Cl": 124,
    "Br": 137 * HCO3_TO_BR, "Form": 1.5, "pH": 7.8,
    "mass_H2O": 1.0, "units": "mg/kgw",
}

CONTAMINANT_LEACHING_BR = 797 * HCO3_TO_BR
CONTAMINANT_IRRIGATION_BR = 260 * HCO3_TO_BR


def get_year_and_week(cell_index):
    year = (cell_index - 1) // WEEKS_PER_YEAR + 1
    week_in_year = (cell_index - 1) % WEEKS_PER_YEAR + 1
    return year, week_in_year


def get_cell_type(cell_index):
    """Weeks are typed by the calendar month their midpoint falls in."""
    _, wk = get_year_and_week(cell_index)
    m = WEEK_MONTH[wk]
    if m in (2, 3):
        return "leaching"
    if m == 6:
        return "irrigation_1"
    if m in (8, 9):
        return "irrigation_2"
    return "rain"



def calculate_drainage_boussinesq(mass_h2o, cell_type, week_in_year, cell_index):
    drain_distance = calculate_drain_distance(cell_index)
    excess_mass = mass_h2o - BASE_MASS
    head_above_base = (excess_mass / BASE_MASS) * GWL_FACTOR
    h_collector = get_collector_level(week_in_year)
    effective_head = head_above_base - h_collector

    if abs(effective_head) < 0.001:
        return 0.0, drain_distance, 0.0, 0.0

    k_h = K_HORIZONTAL_BASE * DRAIN_EFFICIENCY * DRAIN_MULT
    tau = SPECIFIC_YIELD * drain_distance ** 2 / (math.pi ** 2 * k_h * EFFECTIVE_THICKNESS)
    tau = max(tau, 0.001)
    fraction = 1.0 - math.exp(-TIME_STEP / tau)
    head_change = effective_head * fraction
    mass_change = head_change / GWL_FACTOR * BASE_MASS

    if mass_change > 0:
        mass_change = min(mass_change, max(excess_mass, 0))
    else:
        max_inflow_head = h_collector - head_above_base
        max_inflow_mass = max_inflow_head / GWL_FACTOR * BASE_MASS
        mass_change = max(mass_change, -max_inflow_mass)

    return mass_change, drain_distance, tau, fraction


# ----- PHREEQC input builders -----
_SKIP = ("units", "mass_H2O", "mass_H2O_phreeqc", "drainage",
         "drain_distance", "tau_days", "fraction_drained", "overflow")


def _solution_block(num, solution, is_initial=False):
    L = [f"SOLUTION {num}", "    units mg/kgw", "    temp 11.4"]
    if is_initial:
        L += ["    Mg 17", "    Na 90", "    Ca 73", "    K 5.00",
              "    S(6) 166 as SO4", "    Cl 124",
              f"    Br {137 * HCO3_TO_BR:.4f}", "    Form 1.5", "    pH 7.8",
              f"    -water {BASE_MASS}"]
    else:
        for species, value in solution.items():
            if species in _SKIP:
                continue
            if species == "S(6)":
                L.append(f"    {species} {value:.4f} as SO4")
            elif species == "Br":
                L.append(f"    Br {value:.4f}")
            elif species == "Form":
                L.append(f"    {species} {value:.4f}")
            elif species == "pH":
                L.append(f"    pH {value:.2f}")
            else:
                L.append(f"    {species} {value:.4f}")
        mass_h2o = solution.get("mass_H2O", BASE_MASS)
        L.append(f"    -water {mass_h2o:.6f}")
    L.append("END")
    return "\n".join(L) + "\n\n"


def _overflow_solution_block(num, solution):
    L = [f"SOLUTION {num}", "    units mg/kgw", "    temp 11.4"]
    for species, value in solution.items():
        if species in _SKIP:
            continue
        if species == "S(6)":
            L.append(f"    {species} {value:.4f} as SO4")
        elif species == "Br":
            L.append(f"    Br {value:.4f}")
        elif species == "Form":
            L.append(f"    {species} {value:.4f}")
        elif species == "pH":
            L.append(f"    pH {value:.2f}")
        else:
            L.append(f"    {species} {value:.4f}")
    L.append("    -water 1.0")
    L.append("END")
    return "\n".join(L) + "\n\n"


def _contaminant_block(num, mg, na, ca, k, so4, cl, br, form, ph):
    return (f"SOLUTION {num}\n  units mg/kgw\n  temp 11.4\n"
            f"    Mg {mg}\n    Na {na}\n    Ca {ca}\n    K {k}\n"
            f"    S(6) {so4} as SO4\n    Cl {cl}\n"
            f"    Br {br:.4f}\n    Form {form}\n  pH {ph}\nEND\n\n")


def _selected_output_block():
    L = ["SELECTED_OUTPUT",
         "  -high_precision true",
         "  -reset true",
         "  -water true",
         "  -user_punch true",
         "  -selected_output true",
         "",
         "USER_PUNCH",
         "  -headings Mg_mgkgw Na_mgkgw Ca_mgkgw K_mgkgw SO4_mgkgw Cl_mgkgw "
         "Br_mgkgw HCO3_tracer_mgkgw TDS_mgkgw pH CH2O_mgkgw",
         "  -start",
         '  10 PUNCH TOT("Mg") * 24.305 * 1000',
         '  20 PUNCH TOT("Na") * 22.98977 * 1000',
         '  30 PUNCH TOT("Ca") * 40.078 * 1000',
         '  40 PUNCH TOT("K") * 39.0983 * 1000',
         '  50 PUNCH TOT("S(6)") * 96.06 * 1000',
         '  60 PUNCH TOT("Cl") * 35.45 * 1000',
         '  70 PUNCH TOT("Br") * 79.904 * 1000',
         f'  80 PUNCH TOT("Br") * 79.904 * 1000 * {BR_TO_HCO3:.6f}',
         f'  90 PUNCH (TOT("Mg")*24.305 + TOT("Na")*22.98977 + TOT("Ca")*40.078 '
         f'+ TOT("K")*39.0983 + TOT("S(6)")*96.06 + TOT("Cl")*35.45 '
         f'+ TOT("Br")*79.904*{BR_TO_HCO3:.6f}) * 1000',
         '  100 PUNCH -LA("H+")',
         '  110 PUNCH MOL("Form") * 30.026 * 1000',
         "  -end",
         "END",
         ""]
    return "\n".join(L) + "\n"


def _mix_block(base_components, overflow_fraction):
    L = ["MIX 2"]
    for sol_num, frac in base_components:
        L.append(f"    {sol_num} {frac}")
    if overflow_fraction > 0:
        L.append(f"    900 {overflow_fraction:.6f}")
    return "\n".join(L) + "\n"


# Salt held in the soil between the capillary rise that delivered it and the
# leaching campaign that returns it (kg/m2). A one-element list so that the
# builders can mutate it.
SOIL_SALT_STORE = [0.0]
_RETURN_PER_WEEK = [0.0]
_LEACH_WEEKS = sorted(w for w in range(1, WEEKS_PER_YEAR + 1) if WEEK_MONTH[w] in (2, 3))


def build_leaching(index, solution, previous_solution, overflow_fraction, overflow_solution):
    year, week_in_year = get_year_and_week(index)
    contamination = get_water_input("leaching", week_in_year)
    # At the start of the campaign the salt accumulated in the soil over the past
    # year is released and shared equally over the weeks of the campaign.
    if week_in_year == _LEACH_WEEKS[0]:
        _RETURN_PER_WEEK[0] = SOIL_SALT_STORE[0] / len(_LEACH_WEEKS)
        SOIL_SALT_STORE[0] = 0.0
    s = _solution_block(0, previous_solution, is_initial=(index == 1))
    s += _solution_block(1, solution)
    vol_m3 = (LEACHING_WATER_INPUT.get(week_in_year, 0.0) * LEACHING_SCALE
              * SURFACE_DIVERSION_MULT * (1.0 - SSD_INTERCEPT) / 1000.0)
    extra_g_per_l = (_RETURN_PER_WEEK[0] / vol_m3) if vol_m3 > 0 else 0.0
    k = LEACH_CONC_SCALE + extra_g_per_l * 1000.0 / REFERENCE_LEACH_TDS_MGL
    s += _contaminant_block(800, 99*k, 526*k, 425*k, 28*k, 966*k, 722*k,
                            CONTAMINANT_LEACHING_BR*k, 7.6, 7.8)
    if overflow_fraction > 0 and overflow_solution is not None:
        s += _overflow_solution_block(900, overflow_solution)
    s += _selected_output_block()
    s += _mix_block([(1, 1.0), (800, contamination)], overflow_fraction)
    s += "SAVE SOLUTION 3\nEND\n\n"
    s += "EXCHANGE 1\n    X 0.01\nSAVE EXCHANGE 1\nEND\n"
    return s


def _build_with_evap(index, solution, previous_solution, overflow_fraction, overflow_solution,
                     contaminant_args):
    year, week_in_year = get_year_and_week(index)
    reaction_moles = calculate_evaporation_moles(week_in_year, solution["mass_H2O"])
    water_input = get_water_input(get_cell_type(index), week_in_year)
    s = _solution_block(0, previous_solution)
    s += _solution_block(1, solution)
    s += _contaminant_block(700, *contaminant_args)
    if overflow_fraction > 0 and overflow_solution is not None:
        s += _overflow_solution_block(900, overflow_solution)
    s += _selected_output_block()
    s += _mix_block([(1, 1.0), (700, water_input)], overflow_fraction)
    s += "SAVE SOLUTION 3\nEND\n\n"
    s += "USE SOLUTION 3\n"
    s += f"REACTION 1\n    H2O -{reaction_moles:.6f}\n"
    s += "SAVE SOLUTION 4\nEND\n\n"
    s += "EXCHANGE 1\n    X 0.01\nSAVE EXCHANGE 1\nEND\n"
    return s


def build_irrigation_1(index, solution, previous_solution, overflow_fraction, overflow_solution):
    return _build_with_evap(index, solution, previous_solution, overflow_fraction, overflow_solution,
                            (39, 164, 105, 3, 452, 232, CONTAMINANT_IRRIGATION_BR, 7.6, 7.8))


def build_irrigation_2(index, solution, previous_solution, overflow_fraction, overflow_solution):
    return _build_with_evap(index, solution, previous_solution, overflow_fraction, overflow_solution,
                            (39, 164, 105, 3, 452, 232, CONTAMINANT_IRRIGATION_BR, 7.6, 7.8))


def _salt_dissolution_lines(solution):
    """Dissolution of accumulated soil salt once the water table contacts the
    salinized upper profile. Adds salt only, no water."""
    if SALT_MULT <= 0.0:
        return ""
    gwl = calculate_gwl(solution.get("mass_H2O", BASE_MASS))
    if gwl < WASTELAND_SALT_GWL_M:
        return ""
    g = WASTELAND_SALT_G_PER_WEEK * SALT_MULT
    # Dissolve as elements (Na, Cl, S) rather than named phases, so the custom
    # database does not need NaCl / Na2SO4 phase definitions.
    #   NaCl    -> 1 Na + 1 Cl   (M = 58.44)
    #   Na2SO4  -> 2 Na + 1 S    (M = 142.04)
    mol_nacl = g * WASTELAND_SALT_NACL_FRAC / 58.44
    mol_na2so4 = g * (1.0 - WASTELAND_SALT_NACL_FRAC) / 142.04
    # In a REACTION block reactants are given as chemical FORMULAS; the formula
    # fixes stoichiometry and valence unambiguously (so sulfur enters as
    # sulfate, not an ambiguous S that triggers a redox solve), and no phase
    # definition in the database is needed.
    out = ""
    if mol_nacl > 0:
        out += f"    NaCl {mol_nacl:.8f}\n"
    if mol_na2so4 > 0:
        out += f"    Na2SO4 {mol_na2so4:.8f}\n"
    return out


def _build_no_input_evap(index, solution, previous_solution, overflow_fraction, overflow_solution):
    year, week_in_year = get_year_and_week(index)
    reaction_moles = calculate_evaporation_moles(week_in_year, solution["mass_H2O"])
    s = _solution_block(0, previous_solution)
    s += _solution_block(1, solution)
    if overflow_fraction > 0 and overflow_solution is not None:
        s += _overflow_solution_block(900, overflow_solution)
    s += _selected_output_block()
    s += _mix_block([(1, 1.0)], overflow_fraction)
    s += "SAVE SOLUTION 3\nEND\n\n"
    s += "USE SOLUTION 3\n"
    s += f"REACTION 1\n    H2O -{reaction_moles:.6f}\n"
    s += _salt_dissolution_lines(solution)
    s += "SAVE SOLUTION 4\nEND\n\n"
    s += "EXCHANGE 1\n    X 0.01\nSAVE EXCHANGE 1\nEND\n"
    return s





def build_rain(index, solution, previous_solution, overflow_fraction, overflow_solution):
    """Weeks outside the campaigns: precipitation only, plus capillary rise."""
    return _build_with_evap(index, solution, previous_solution,
                            overflow_fraction, overflow_solution, RAIN_COMPOSITION)


BUILDERS = {
    "rain": build_rain,
    "leaching": build_leaching,
    "irrigation_1": build_irrigation_1,
    "irrigation_2": build_irrigation_2,
}


# ----- IPhreeqcPy runner -----
class IPhreeqcRunner:
    def __init__(self, database_path):
        import IPhreeqcPy
        self.phreeqc = IPhreeqcPy.IPhreeqc()
        self.phreeqc.LoadDatabase(database_path)

    def _error_string(self):
        try:
            return self.phreeqc.GetErrorString()
        except TypeError:
            return self.phreeqc.GetErrorString
        except Exception:
            return ""

    def run(self, input_string):
        self.phreeqc.RunString(input_string)
        arr = self.phreeqc.GetSelectedOutputArray()
        if not arr or len(arr) < 2:
            raise RuntimeError("PHREEQC returned no selected output.\n"
                               + str(self._error_string()) + "\n--- input ---\n" + input_string)
        return arr


def array_to_last_dict(arr):
    headers = arr[0]
    last = arr[-1]
    d = {}
    for name, val in zip(headers, last):
        if name not in d:
            d[name] = val
    return d


def calculate_gwl(mass_h2o):
    return BASE_GWL + (mass_h2o - BASE_MASS) * GWL_FACTOR


def update_solution_from_array(arr, cell_type, week_in_year, cell_index):
    row = array_to_last_dict(arr)
    if "mass_H2O" not in row:
        raise KeyError("No 'mass_H2O' column. Got: " + str(list(row.keys())))
    mass_h2o_phreeqc = row["mass_H2O"]
    drainage, drain_distance, tau, fraction = calculate_drainage_boussinesq(
        mass_h2o_phreeqc, cell_type, week_in_year, cell_index)
    mass_h2o_after_drain = max(mass_h2o_phreeqc - drainage, MIN_MASS)
    overflow = 0.0
    if mass_h2o_after_drain > MAX_MASS:
        overflow = mass_h2o_after_drain - MAX_MASS
        mass_h2o_after_drain = MAX_MASS
    # Density-driven loss of dissolved salt to the stagnant lower zone.
    # Applied to the solutes only; the water mass is untouched, so the water
    # balance and the water-table position are unaffected.
    keep = 1.0 - max(0.0, min(DOWNWARD_SALT_LOSS_PER_WEEK, 0.5))
    # Capillary rise carries water AND its dissolved salt into the profile, so
    # the fraction of the cell water drawn upward is removed from the solutes as
    # well. The salt is stored and returned to the aquifer by the next leaching
    # campaign, closing the exchange between aquifer and root zone.
    up_mm = capillary_flux_mm(week_in_year, mass_h2o_after_drain)
    up_frac = min(0.5, mm_to_mass_fraction(up_mm) / max(mass_h2o_after_drain, 0.1))
    tds = (row["Mg_mgkgw"] + row["Na_mgkgw"] + row["Ca_mgkgw"] + row["K_mgkgw"]
           + row["SO4_mgkgw"] + row["Cl_mgkgw"] + row["Br_mgkgw"] * BR_TO_HCO3)
    SOIL_SALT_STORE[0] += (up_frac * tds / 1000.0 * mass_h2o_after_drain
                           * WATER_IN_CELL / CELL_AREA / 1000.0)   # kg/m2
    keep *= (1.0 - up_frac)
    return {
        "Mg": row["Mg_mgkgw"] * keep, "Na": row["Na_mgkgw"] * keep,
        "Ca": row["Ca_mgkgw"] * keep, "K": row["K_mgkgw"] * keep,
        "S(6)": row["SO4_mgkgw"] * keep, "Cl": row["Cl_mgkgw"] * keep,
        "Br": row["Br_mgkgw"] * keep, "Form": row["CH2O_mgkgw"], "pH": row["pH"],
        "mass_H2O": mass_h2o_after_drain, "mass_H2O_phreeqc": mass_h2o_phreeqc,
        "drainage": drainage, "drain_distance": drain_distance,
        "tau_days": tau, "fraction_drained": fraction,
        "overflow": overflow, "units": "mg/kgw",
    }


def run_simulation(runner, total_cells=TOTAL_CELLS):
    SOIL_SALT_STORE[0] = 0.0
    _RETURN_PER_WEEK[0] = 0.0
    current_solution = initial_solution.copy()
    all_results = []
    pending_overflow = 0.0
    pending_overflow_solution = None
    for i in range(1, total_cells + 1):
        cell_type = get_cell_type(i)
        year, week_in_year = get_year_and_week(i)
        previous_solution = current_solution.copy()
        overflow_fraction = pending_overflow
        overflow_solution = pending_overflow_solution
        input_string = BUILDERS[cell_type](
            i, current_solution, previous_solution, overflow_fraction, overflow_solution)
        arr = runner.run(input_string)
        current_solution = update_solution_from_array(arr, cell_type, week_in_year, i)
        new_overflow = current_solution.get("overflow", 0.0)
        if new_overflow > 0:
            pending_overflow = new_overflow
            pending_overflow_solution = {k: v for k, v in current_solution.items() if k not in _SKIP}
        else:
            pending_overflow = 0.0
            pending_overflow_solution = None
        gwl = calculate_gwl(current_solution["mass_H2O"])
        tds = (current_solution["Mg"] + current_solution["Na"] + current_solution["Ca"]
               + current_solution["K"] + current_solution["S(6)"] + current_solution["Cl"]
               + current_solution["Br"] * BR_TO_HCO3)
        all_results.append({
            "Cell_Number": i, "Year": year, "Week_in_Year": week_in_year,
            "Cell_Type": cell_type, "GWL_m": gwl, "TDS_mgkgw": tds,
            "mass_H2O_kg": current_solution["mass_H2O"],
            "Drain_Distance_m": current_solution.get("drain_distance"),
            "Tau_days": current_solution.get("tau_days"),
        })
    return all_results


def print_calibration_report(results):
    """Compare LAST simulated year's depth-to-water to observed well-45."""
    last = [r for r in results if r["Year"] == YEARS]
    print("\n" + "=" * 64)
    print(f"CALIBRATION REPORT — well 45 | knobs: leach={LEACHING_SCALE} "
          f"irr={IRRIGATION_SCALE} evap={EVAP_SCALE} coll={COLLECTOR_SCALE} "
          f"drain_eff={DRAIN_EFFICIENCY} L={FIXED_DRAIN_DISTANCE}")
    print("=" * 64)
    print(f"{'wk':>3} {'season':>18} {'GWL_m':>7} {'DTW_mod':>8} {'DTW_obs':>8} {'diff':>7} {'TDS':>7}")
    sq = 0.0
    n = 0
    for r in last:
        wk = r["Week_in_Year"]
        dtw_mod = SURFACE_ELEVATION - r["GWL_m"]
        dtw_obs = OBSERVED_DTW_WELL45[wk]
        diff = dtw_mod - dtw_obs
        sq += diff * diff
        n += 1
        print(f"{wk:>3} {r['Cell_Type']:>18} {r['GWL_m']:>7.2f} {dtw_mod:>8.2f} "
              f"{dtw_obs:>8.2f} {diff:>7.2f} {r['TDS_mgkgw']:>7.0f}")
    rmse = math.sqrt(sq / n)
    dtw_mod_all = [SURFACE_ELEVATION - r["GWL_m"] for r in last]
    print("-" * 64)
    print(f"RMSE depth-to-water (model vs obs): {rmse:.2f} m")
    print(f"model DTW: min {min(dtw_mod_all):.2f}  max {max(dtw_mod_all):.2f}  "
          f"mean {sum(dtw_mod_all)/len(dtw_mod_all):.2f} m")
    print(f"obs   DTW: min 1.50  max 3.13  mean ~2.56 m   (well 45)")
    print(f"equilibrium TDS (last-year mean): "
          f"{sum(r['TDS_mgkgw'] for r in last)/len(last):.0f} mg/l "
          f"(target ~2500-3000)")
    print("=" * 64)
    print(">>> Copy this whole block back to me. <<<")


# ============================================================
# SCENARIO DEFINITIONS
# ============================================================
# Each scenario overrides a few globals; everything else stays at the
# calibrated values. "drain_mult" scales the calibrated drainage response.
SCENARIOS = [
    # name                                land_use  drain_mult  ssd   salt  evap   surf   coll
    ("Present state",                     "field",       1.00,  0.00,  0.0,  1.0,  1.000, 1.00),
    ("Both functions maintained",         "field",       1.50,  0.00,  0.0,  1.0,  0.667, 0.80),
    ("Furrows maintained, drains silted", "field",       0.25,  0.00,  0.0,  1.0,  0.667, 1.20),
    ("Furrows silted, drains maintained", "field",       1.50,  0.00,  0.0,  1.0,  1.500, 0.80),
    ("Both functions degraded",           "field",       0.25,  0.00,  0.0,  1.0,  1.500, 1.20),
    ("Shallow subsurface drainage",       "field",       1.00,  0.70,  0.0,  1.0,  1.000, 1.00),
    ("Plastic film on the field",         "field",       1.00,  0.00,  0.0,  0.1,  1.000, 1.00),

]

CRITICAL_GWL_M = 69.0   # report metric: share of time with GWL at or above this level


def apply_scenario(land_use, drain_mult, ssd, salt_mult, evap, surf=1.0, coll=1.0):
    """Set the scenario multipliers. `coll` scales the seasonal collector stage,
    which together with `drain_mult` represents the subsurface function of the
    network; `surf` represents its surface function."""
    g = globals()
    g["LAND_USE"] = land_use
    g["DRAIN_MULT"] = drain_mult
    g["SSD_INTERCEPT"] = ssd
    g["SALT_MULT"] = salt_mult
    g["EVAP_SCALE"] = evap
    g["SURFACE_DIVERSION_MULT"] = surf
    g["COLLECTOR_SCALE"] = coll


def scenario_metrics(results):
    """Metrics in the format of report Table 12."""
    last = [r for r in results if r["Year"] == YEARS]
    dtw = [SURFACE_ELEVATION - r["GWL_m"] for r in last]
    tds_final = results[-1]["TDS_mgkgw"]
    tds_mean = sum(r["TDS_mgkgw"] for r in last) / len(last)
    crit = 100.0 * sum(1 for r in results if r["GWL_m"] >= CRITICAL_GWL_M) / len(results)
    return dict(tds_final_mgL=tds_final, tds_year5_mean_mgL=tds_mean,
                min_depth_to_surface_m=min(dtw), mean_depth_m=sum(dtw) / len(dtw),
                time_critical_gwl_pct=crit)


def blending_potential(tds_gL):
    """Report section 5.7: share of irrigation events 2-3 that groundwater can
    cover while keeping the blend at or below the target TDS."""
    if tds_gL <= TARGET_BLEND_TDS_GL:
        share = 1.0
    elif tds_gL <= SURFACE_WATER_TDS_GL:
        share = 1.0
    else:
        share = (TARGET_BLEND_TDS_GL - SURFACE_WATER_TDS_GL) / (tds_gL - SURFACE_WATER_TDS_GL)
        share = max(0.0, min(1.0, share))
    saved_m3_ha = share * POST_FLOWERING_MM * 10.0
    return share, saved_m3_ha


def water_balance_rmse(results):
    """Depth-to-water RMSE of the final year against the well-45 observations."""
    last = [r for r in results if r["Year"] == YEARS]
    sq = 0.0
    for r in last:
        sq += (SURFACE_ELEVATION - r["GWL_m"] - OBSERVED_DTW_WELL45[r["Week_in_Year"]]) ** 2
    return (sq / len(last)) ** 0.5


def season_tds_mgL(results, month):
    """Mean simulated TDS over the weeks of the given month, final year.

    The monitoring network measures mineralization in April, July and October.
    April is the first measurement after the Feb-March leaching campaign and
    therefore carries the leaching signal most directly: it reflects the salt
    that had accumulated in the profile over the preceding year. July and
    October show how mineralization then evolves within the current season.
    """
    weeks = [wk for wk in range(1, WEEKS_PER_YEAR + 1) if _week_month[wk] == month]
    vals = [r["TDS_mgkgw"] for r in results
            if r["Year"] == YEARS and r["Week_in_Year"] in weeks]
    return sum(vals) / len(vals) if vals else float("nan")


def april_tds_mgL(results):
    return season_tds_mgL(results, "Apr")


def seasonal_tds(results):
    return {m: season_tds_mgL(results, m) for m in ("Apr", "Jul", "Oct")}


# Observed seasonal mineralization at well 45, 2025 (g/L)
OBS_SEASONAL_GL = {"Apr": 3.56, "Jul": 3.70, "Oct": 1.66}


def seasonal_amplitude(results):
    """Ratio of the highest to the lowest simulated seasonal value.

    Reported for comparison with the monitoring network, whose median ratio is
    about 1.4. MIXING_THICKNESS is set from the screened interval of the wells
    and is not calibrated against this ratio; the sweep below is a sensitivity
    check on that choice.
    """
    v = seasonal_tds(results)
    lo = min(v.values())
    return (max(v.values()) / lo) if lo > 0 else float("nan")


# ============================================================
# OPTIONAL: soil salinity sweep against the observed April mineralization
# Set SWEEP_SOIL_SALINITY = True to calibrate SOIL_SALT_KG_M3, then set the
# chosen value at the top of this file and run the scenarios normally.
# ============================================================
# ============================================================
# OPTIONAL: mixing-thickness sweep against the observed seasonal amplitude.
# The water balance is independent of MIXING_THICKNESS, so the RMSE must stay
# constant across the sweep - that is the check that the refactor is sound.
# ============================================================
# ============================================================
# OPTIONAL: calibrate the downward salt loss against the observed mean
# mineralization. Run this with YEARS large enough for equilibrium (>= 30).
# ============================================================
SWEEP_SALT_LOSS = True
SWEEP_LOSS_VALUES = [float(x) for x in os.environ.get("SWEEP","0.010,0.020").split(",")]
OBS_MEAN_TDS_MGL = 2620.0


def run_salt_loss_sweep(runner):
    global DOWNWARD_SALT_LOSS_PER_WEEK
    keep = DOWNWARD_SALT_LOSS_PER_WEEK
    apply_scenario("field", 1.0, 0.0, 0.0, 1.0, 1.0)
    print(f"Downward salt-loss sweep — {YEARS}-year runs, target "
          f"{OBS_MEAN_TDS_MGL:.0f} mg/L")
    print(f"  {'loss/wk':>8} {'loss/yr':>8} {'annual':>8} {'resid':>8} "
          f"{'Apr':>7} {'Oct':>7} {'RMSE m':>8}")
    for k in SWEEP_LOSS_VALUES:
        DOWNWARD_SALT_LOSS_PER_WEEK = k
        res = run_simulation(runner, TOTAL_CELLS)
        m = scenario_metrics(res)
        sea = seasonal_tds(res)
        ann = m["tds_year5_mean_mgL"]
        print(f"  {k:8.3f} {1-(1-k)**WEEKS_PER_YEAR:8.2f} {ann:8.0f} "
              f"{ann-OBS_MEAN_TDS_MGL:+8.0f} {sea['Apr']:7.0f} {sea['Oct']:7.0f} "
              f"{water_balance_rmse(res):8.2f}")
    DOWNWARD_SALT_LOSS_PER_WEEK = keep
    print("\n  RMSE must be identical in every row: the term removes salt, not water.")
    print("  Pick the rate matching the target, set it at the top, then")
    print("  recalibrate SOIL_SALT_KG_M3 if needed.\n")


SWEEP_MIXING_THICKNESS = False
SWEEP_MIXING_VALUES_M = [3.0, 4.0, 6.0, 9.0, 14.0, 22.0]


def run_mixing_sweep(runner):
    global MIXING_THICKNESS, WATER_IN_CELL, GWL_FACTOR
    keep = MIXING_THICKNESS
    apply_scenario("field", 1.0, 0.0, 0.0, 1.0, 1.0)
    obs = max(OBS_SEASONAL_GL.values()) / min(OBS_SEASONAL_GL.values())
    print("Mixing-thickness sweep (baseline scenario)")
    print(f"  observed seasonal amplitude at well 45 = {obs:.2f}")
    print(f"  {'mix m':>7} {'Apr':>7} {'Jul':>7} {'Oct':>7} {'ampl':>7} "
          f"{'annual':>8} {'RMSE m':>8}")
    for mt in SWEEP_MIXING_VALUES_M:
        MIXING_THICKNESS = mt
        WATER_IN_CELL = CELL_AREA * MIXING_THICKNESS * POROSITY * 1000
        GWL_FACTOR = MIXING_THICKNESS * POROSITY / STORAGE_COEFF
        derive_leach_conc_scale(verbose=False)
        res = run_simulation(runner, TOTAL_CELLS)
        sea = seasonal_tds(res)
        m = scenario_metrics(res)
        print(f"  {mt:7.1f} {sea['Apr']:7.0f} {sea['Jul']:7.0f} {sea['Oct']:7.0f} "
              f"{seasonal_amplitude(res):7.2f} {m['tds_year5_mean_mgL']:8.0f} "
              f"{water_balance_rmse(res):8.2f}")
    MIXING_THICKNESS = keep
    WATER_IN_CELL = CELL_AREA * MIXING_THICKNESS * POROSITY * 1000
    GWL_FACTOR = MIXING_THICKNESS * POROSITY / STORAGE_COEFF
    derive_leach_conc_scale(verbose=False)
    print("\n  RMSE must be identical in every row: the water balance does not")
    print("  depend on the mixing thickness. Pick the row whose amplitude matches,")
    print("  set MIXING_THICKNESS at the top, then recalibrate SOIL_SALT_KG_M3.\n")


SWEEP_SOIL_SALINITY = False
SWEEP_VALUES_KG_M3 = [10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]

# Observed April mineralization at well 45 (g/L): 2.66 (2023), 3.38 (2024), 3.56 (2025).
# The model reproduces one stationary year, so the three-year mean is the target.
# The spread between years is +-14% of the mean and carries directly into the
# inferred soil salinity, which should be quoted with that uncertainty.
OBS_TDS_APRIL_MGL = 3200.0


def run_soil_salinity_sweep(runner):
    global SOIL_SALT_KG_M3
    keep = SOIL_SALT_KG_M3
    apply_scenario("field", 1.0, 0.0, 0.0, 1.0, 1.0)
    print("Soil salinity sweep (baseline scenario, April TDS of the final year)")
    print(f"  percolate volume {leaching_percolate_mm_per_season():.0f} mm/season")
    if OBS_TDS_APRIL_MGL is not None:
        print(f"  observed April TDS = {OBS_TDS_APRIL_MGL:.0f} mg/L")
    print(f"  {'soil kg/m3':>10} {'percolate g/L':>14} {'Apr':>8} {'Jul':>8} {'Oct':>8} "
          f"{'final':>8} {'Apr resid':>10} {'RMSE m':>8}")
    rows = []
    for s in SWEEP_VALUES_KG_M3:
        SOIL_SALT_KG_M3 = s
        conc = derive_leach_conc_scale(verbose=False) * REFERENCE_LEACH_TDS_MGL / 1000.0
        res = run_simulation(runner, TOTAL_CELLS)
        sea = seasonal_tds(res)
        fin = scenario_metrics(res)["tds_final_mgL"]
        rmse = water_balance_rmse(res)
        resid = ("" if OBS_TDS_APRIL_MGL is None
                 else f"{sea['Apr'] - OBS_TDS_APRIL_MGL:+.0f}")
        print(f"  {s:10.1f} {conc:14.2f} {sea['Apr']:8.0f} {sea['Jul']:8.0f} "
              f"{sea['Oct']:8.0f} {fin:8.0f} {resid:>10} {rmse:8.2f}")
        rows.append((s, conc, sea["Apr"], sea["Jul"], sea["Oct"], fin))
    SOIL_SALT_KG_M3 = keep
    derive_leach_conc_scale(verbose=False)
    print("\n  The water balance (RMSE) must not change across the sweep: the salt")
    print("  balance alters concentrations only, never the volumes.\n")
    return rows


# ============================================================
# RUN
# ============================================================
