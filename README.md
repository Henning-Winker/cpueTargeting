# cpueTargeting

`cpueTargeting` provides simulation tools for testing how CPUE standardisation
copes with changes in fishing targeting. It contains three operating models:

1. `simulate_cluster()`: the multi-habitat cluster example of Winker et al. (2014).
2. `simulate_bathymetric()`: the four-species Run 22 model used for the
   WKBMSYSPICT CPUE guidelines. Species follow bathymetry; skippers target
   species through *where* they fish.
3. `simulate_tuna()`: a longline tuna model that adds *how deep* they fish.
   Hook depth by hooks between floats (HBF) follows the habitat-based model of
   Bigelow et al. (2002).

## Installation

```r
remotes::install_github("henning-winker/cpueTargeting")
```

The simulators only need `ggplot2`; `mgcv` is needed for
`standardize_tuna_indices()`.

## Quick start

```r
library(cpueTargeting)

# Run 22 bathymetric model
run22 <- simulate_bathymetric(seed = 2)
run22$Plots$CPUE          # nominal CPUE versus true biomass

# Tuna model with vertical targeting
tuna <- simulate_tuna(seed = 3)
tuna$Plots$Env            # SST and thermocline anomaly
tuna$Plots$SpDens         # species densities
tuna$Plots$HBF            # gear use over time
tuna$Plots$CPUE           # nominal CPUE versus true biomass
plot_tuna_vertical_distributions()

# Standardise with and without HBF
comparison <- standardize_tuna_indices(tuna, spatial = c(FALSE, TRUE), spatial_res = 5)
comparison$diagnostics
plot_indices(comparison, models = c("Base+S", "DPC+S", "Cluster+S", "HBF+S", "HBS+S"))

# One row per set (e.g. for your own GAMs)
sets <- as_set_data(tuna$Data)
```

## The tuna operating model

The ocean is a 50 x 50 grid in a conceptual Southern Hemisphere: north (top) is
warm and south is cool.

| Process | Driven by | Main settings |
|---|---|---|
| Horizontal distribution | SST field: north-south gradient + random field | `sst_opt`, `sst_tol`, `sst_gradient`, `sst_noise` |
| Movement | Exponential distance kernel, weighted by habitat | `dispersal` |
| Vertical distribution | Normal depth profile, shifted by a thermocline anomaly (fixed field + optional annual shift and annual field) | `depth_mu`, `depth_sd`, `thermocline_sd`, `thermocline_annual_sd`, `thermocline_annual_field_sd` |
| Hook depth | Catenary by HBF; mainline angle varies between sets; hooks vary within sets; current shoaling (Bigelow et al. 2002, eq. 2-5) | `HBF_levels`, `phi_mean`, `phi_sd`, `current_speed` |
| Catchability of a set | Vessel q x overlap between hook and fish depth distributions | `q`, `vessel_cv`, `tweedie_power`, `tweedie_phi` |
| Spatial targeting | Effort goes to cells with high expected catch value | `price`, `preference` |
| Vertical targeting | Skippers choose the HBF with the highest expected value | `HBF_choice_temp` |
| Population | Schaefer model in each cell, then movement | `B0`, `r` |

Default species:

| Species | Water | Depth | `sst_opt` | `depth_mu` (m) |
|---|---|---|---|---|
| YFT (yellowfin) | warm | shallow | 25 | 90 |
| BET (bigeye) | warm | deep | 24 | 200 |
| BFT (southern bluefin) | cool | deep | 20 | 180 |
| ALB (albacore) | cool | shallow | 21 | 110 |

The pairs overlap in space and the species overlap in depth, so most sets catch
all four. Prices favour YFT in years 1-10, BET in years 11-20 and BFT in years
21-30, so the fleet moves between warm and cool water and between shallow and
deep gear.

`Data` also contains the latent truth that analysts would not see:
`hook_depth` (mean hook depth of the set), `vertical_avail` (realised vertical
overlap) and `hbs_avail` (overlap predicted from HBF, as used in a
habitat-based standardisation). Set `thermocline_sd = 0` to make vertical
availability depend on HBF only, so that an HBF factor can fully correct for
it.

By default the thermocline pattern is fixed over time. Setting
`thermocline_annual_sd` (a basin-wide AR(1) shift) and/or
`thermocline_annual_field_sd` (a new spatial field each year) makes it move
between years, so the catchability of a given HBF changes from year to year.
The vignette *Vertical targeting and a moving thermocline* walks through
constant, fixed and moving scenarios step by step.

Plots that illustrate each effect:

| Function | Shows |
|---|---|
| `plot_tuna_vertical_distributions()` | species depth profiles and hook profiles by HBF |
| `plot_thermocline_effect()` | how a shallow or deep thermocline moves species relative to the hooks |
| `plot_thermocline(sim)` | annual thermocline series and maps |
| `plot_availability(sim)` | vertical availability by HBF over time |
| `sim$Plots$HBF`, `sim$Plots$CPUE` | gear use over time; nominal CPUE versus truth |
| `plot_indices()`, `plot_diagnostics()` | standardised indices and their accuracy, for one or several scenarios |

To change a scenario, edit the settings list:

```r
s <- tuna_settings()
s$depth_mu <- c(80, 250, 220, 100)   # stronger vertical separation
s$HBF_choice_temp <- 0.05            # skippers choose gear more strictly
sim <- simulate_tuna(s, seed = 1)
```

## Speed

A default 30-year run takes a few seconds. The original code was slow because
it built each species' 2500 x 2500 movement matrix cell by cell, iterated the
movement 1000 times for the initial distribution, and grew data frames with
`rbind()` inside the cell loop. The package uses the closed-form stable
distribution, never forms the movement matrix explicitly, and simulates all
sets of a year at once.

## Original code

The original Run 22 scripts (`OMcpueSIM22_original.R`,
`Generate_scenario_data22_original.R`) and the 2014 cluster simulator are kept
unchanged in `inst/extdata` for reference. See `NEWS.md` for differences.

## References

Bigelow, K.A., Hampton, J. and Miyabe, N. (2002) Application of a habitat-based
model to estimate effective longline fishing effort and relative abundance of
Pacific bigeye tuna (*Thunnus obesus*). Fisheries Oceanography 11, 143-155.

Thorson, J.T., Fonner, R., Haltuch, M.A., Ono, K. and Winker, H. (2017)
Accounting for spatiotemporal variation and fisher targeting when estimating
abundance from multispecies fishery data. Canadian Journal of Fisheries and
Aquatic Sciences 74, 1794-1807.

Winker, H., Kerwath, S.E. and Attwood, C.G. (2013) Comparison of two approaches
to standardize catch-per-unit-effort for targeting behaviour in a
multispecies hand-line fishery. Fisheries Research 139, 118-131.

Winker, H., Kerwath, S.E. and Attwood, C.G. (2014) Proof of concept for a novel
procedure to standardize multispecies catch and effort data. Fisheries
Research 155, 149-159.
