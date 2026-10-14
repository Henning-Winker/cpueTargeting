# cpueTargeting

`cpueTargeting` provides simulation tools for testing how CPUE standardisation
copes with changes in fishing targeting. It contains three operating models:

1. `simulate_cluster()`: the multi-habitat cluster example of Winker et al. (2014).
2. `simulate_bathymetric()`: the four-species model used for the
   WKBMSYSPICT CPUE guidelines. It re-implements the spatial multispecies
   simulator of Thorson et al. (2017, Appendix A; `Generate_scenario_data` in
   [spatial_DFA](https://github.com/James-Thorson/spatial_DFA)), which
   generalises the simulator of Ono et al. (2015). Species follow bathymetry;
   skippers target species through *where* they fish.
3. `simulate_tuna()`: a longline tuna model that adds *how deep* they fish.
   Hook depth depends on hooks between floats (HBF), and habitat-based
   standardisation follows Bigelow et al. (2002).

## A realistic simulated logbook

`simulate_logbook()` simulates a mixed demersal-pelagic fishery with ten
species and returns a messy trip-level logbook together with the truth.
Demersal species live along a depth gradient and on patchy reef or sediment
substrate (shallow/deep reef, shallow/deep sediment, as in Winker et al.
2014); pelagic species follow SST and migrate seasonally. Each trip chooses a
tactic (`Reef` handline, `Sediment` longline or `Pelagic` troll) and a place
by expected value; prices change between periods, a closure bans sediment
gear from deep water, specialist (80% loyal to their home tactics) and mixed vessels enter and leave. The
logbook reports trips (some combining two sets), coarse and sometimes missing
positions, rounded catches, under-reported bycatch and a species under two
names, but not the tactic. See the vignette *A realistic logbook: reef,
sediment and pelagic tactics*.

```r
sim <- simulate_logbook(seed = 1)
sim$Plots$density; sim$Plots$tactics; sim$Plots$cpue
head(sim$logbook)
plot_index(fit, true = logbook_truth(sim, "Pagellus bogaraveo"))
```

## Toolbox for real data

The same methods can be applied step by step to real logbook or observer
data. The vignette *Analysing real catch and effort data: a step-by-step
guide* (`vignettes/real-data-guide.Rmd`) walks through each step with
guidelines for the decisions involved.

| Step | Functions |
|---|---|
| 1. Format: long or wide logbooks to one row per operation, synonyms, processing log | `cpue_data()`, `subset()`, `cpue_matrix()` |
| 2. Fleet: operations by year and vessel, fleet turnover, core vessels | `plot_fleet()`, `select_vessels()` |
| 3. Species: cumulative species curves and the species set | `species_curve()`, `plot_species_curve()`, `select_species()` |
| 4. Nominal CPUE and catch composition | `nominal_cpue()`, `plot_nominal()`, `plot_composition()`, `plot_cpue_map()` |
| 5. Habitat-based subsetting (Stephens and MacCall 2004) | `stephens_maccall()`, `plot_stephens_maccall()` |
| 6. Targeting covariates: composition PCA, Optimal Coordinates test, clustering | `add_targeting()`, `composition_pca()`, `compute_oc()`, `cluster_tactics()`, `plot_oc()`, `plot_biplot()`, `tactic_table()`, `plot_targeting()` |
| 7. Models (mgcv): Base, Covariates, DPC, Cluster, spatial (`+S`), random vessel effects; Tweedie or delta families | `fit_cpue()` |
| 8. Diagnostics: model table, influence (Bentley et al. 2012), step plot, quantile residuals, effects, vessel effects | `cpue_diagnostics()`, `plot_index()`, `plot_step()`, `plot_influence()`, `plot_residuals()`, `plot_effects()`, `plot_vessel_effects()` |
| 9. Index for the assessment, with log SE and CV | `cpue_index()`, `export_index()` |

```r
x <- cpue_data(logbook, year = "year", species = "species", catch = "kg",
               effort = "hours", id = "trip", vessel = "vessel",
               covariates = "month",
               synonyms = c("Trichiurus lepturus" = "Lepidopus caudatus"))
plot_species_curve(x, threshold = 0.95)
x <- select_species(x, threshold = 0.95, keep = "Pagellus bogaraveo")
x <- add_targeting(x)          # PCs (OC test) and n_pc + 1 clusters
plot_oc(x); plot_biplot(x); tactic_table(x)
plot_targeting(x)$time
x$log                          # every selection step, for reporting

fit <- fit_cpue(x, species = "Pagellus bogaraveo",
                models = c("Base", "Covariates", "Cluster", "DPC"),
                terms = 's(month, bs = "cc", k = 8)', vessel = TRUE,
                knots = list(month = c(0.5, 12.5)))
cpue_diagnostics(fit)
plot_step(fit); plot_influence(fit, "DPC"); plot_residuals(fit, "DPC")$qq
export_index(fit, "cpue_index.csv", model = "DPC")
```

## Documentation and worked examples

| Page | Content |
|---|---|
| [Background document](https://henning-winker.github.io/cpueTargeting/appendix.html) | Full description of both simulators, how the tuna model extends the spatial simulator, and the standardisation models (Base, spatial, Cluster, DPC, HBF) |
| [Bathymetric example](https://henning-winker.github.io/cpueTargeting/bathymetric-example.html) | Targeting along a continuous depth gradient with the guideline settings: maps, CPUE distributions, indices, what the PCs and clusters capture, replicates |
| [Tuna longline: flat vs contrasting trends](https://henning-winker.github.io/cpueTargeting/tuna-flat-vs-contrast.html) | Spatial and vertical (HBF) targeting with identical targeting but flat or diverging true trends: when the catch composition captures targeting and when it also absorbs abundance |

All pages are also linked from the [project page](https://henning-winker.github.io/cpueTargeting/).

## Installation

```r
remotes::install_github("henning-winker/cpueTargeting")
```

The simulators only need `ggplot2`; `mgcv` is needed for
`standardize_cpue()`.

## Quick start

```r
library(cpueTargeting)

# Bathymetric model (settings of the CPUE guidelines)
bat <- simulate_bathymetric(seed = 123)
bat$Plots$SpDens          # species densities
bat$Plots$CPUE            # nominal CPUE versus true biomass
fit_bat <- standardize_cpue(bat, models = c("Base", "Base+S", "Cluster+S", "DPC", "DPC+S"),
                            n_pc = 2, n_clusters = 3, spatial_k = 8)
plot_indices(fit_bat, type = "smooth")

# Tuna model with vertical targeting
tuna <- simulate_tuna(seed = 123)
tuna$Plots$Env            # SST and thermocline anomaly
tuna$Plots$SpDens         # species densities
tuna$Plots$HBF            # gear use over time
tuna$Plots$CPUE           # nominal CPUE versus true biomass
plot_tuna_vertical_distributions()

# Standardise: Base, space, clusters, PCs and HBF
comparison <- standardize_cpue(tuna, spatial_res = 5, n_pc = 2, n_clusters = 4,
  models = c("Base", "Base+S", "Cluster+S", "DPC+S", "HBF+S"), engine = "bam")
comparison$diagnostics
plot_indices(comparison, type = "smooth")
plot_pcs(comparison)$time         # what the PCs capture over time
plot_clusters(comparison)$time    # cluster shares over time

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
| Hook depth | Typical (median) depth per HBF level (cf. Bigelow et al. 2002, Fig. 2); lognormal set-level error; lognormal spread of hooks within a set | `HBF_levels`, `HBF_depth`, `set_cv`, `hook_cv` |
| Catchability of a set | Vessel q x overlap between hook and fish depth distributions; fishing catchability `q` sets how hard stocks are fished, index catchability `q_index` sets the zeros and noise of the recorded CPUE | `q`, `q_index`, `vessel_cv`, `tweedie_power`, `tweedie_phi` |
| Spatial targeting | Effort goes to cells with high expected catch value | `price`, `preference` |
| Vertical targeting | Skippers choose the HBF with the highest expected value | `HBF_choice_temp` |
| Population | Schaefer model in each cell, then movement; stocks start at a fraction of carrying capacity | `B0`, `r`, `B_start` |

Default species:

| Species | Water | Depth | `sst_opt` | `depth_mu` (SD) in m |
|---|---|---|---|---|
| YFT (yellowfin) | warm | shallow | 25 | 90 (60) |
| BET (bigeye) | warm | deep | 24 | 260 (50) |
| BFT (southern bluefin) | cool | deep | 20 | 240 (55) |
| ALB (albacore) | cool | shallow | 21 | 110 (60) |

The pairs overlap in space and the species overlap in depth, so most sets catch
all four. Prices favour ALB in years 1-10 (shallow, cool water), BET in years
11-20 (deep, warm water) and both BET and BFT in years 21-30, loosely following
the historical shift from albacore to bigeye targeting. Effort is cut by 25% in
years 21-30. Stocks start partly fished (`B_start`) and the true trends differ:
ALB is fished down, BET declines, BFT recovers and declines again once
targeted, YFT dips and recovers. `tuna_settings("mild")` gives weaker
targeting.

The tuna example uses different prices, which also switch the fleet between
water masses (ALB and BFT in the cool south, then BET and YFT in the warm north,
then BET and BFT), and compares stocks with contrasting true trends (the
defaults) with lightly fished, unfished stocks (`q = 0.2`, `B_start = 1`) under
identical targeting. Both share the same observation model (`q_index = 0.2`);
see `inst/scripts/talk_setup.R`.

`Data` also contains the latent truth that analysts would not see:
`hook_depth` (median hook depth of the set), `vertical_avail` (realised vertical
overlap) and `hbs_avail` (overlap predicted by a habitat model from HBF and the average
thermocline, as used in a habitat-based standardisation). Set `thermocline_sd = 0` to make vertical
availability depend on HBF only, so that an HBF factor can fully correct for
it.

By default the thermocline pattern is fixed over time. Setting
`thermocline_annual_sd` (a basin-wide AR(1) shift) and/or
`thermocline_annual_field_sd` (a new spatial field each year) makes it move
between years, so the catchability of a given HBF changes from year to year;
about 15 m and 10 m give a realistic ENSO-like range. The vignette *Targeting
in a longline tuna fishery* walks through the settings and a fixed and a moving
thermocline scenario step by step; the vignette *Tuna longline targeting: flat
versus contrasting true trends* shows the main results.

Plots that illustrate each effect:

| Function | Shows |
|---|---|
| `plot_tuna_vertical_distributions()` | species depth profiles and hook profiles by HBF |
| `plot_thermocline_effect()` | how a shallow or deep thermocline moves species relative to the hooks |
| `plot_thermocline(sim)` | annual thermocline series and maps |
| `plot_availability(sim)` | vertical availability by HBF over time |
| `sim$Plots$HBF`, `sim$Plots$CPUE` | gear use over time; nominal CPUE versus truth |
| `plot_sections(sim)` | species distributions in 3D: latitude x depth sections and maps by depth layer |
| `plot_pcs(fit)` | what the catch-composition PCs capture: loadings, maps, scores against hook depth, scores over time |
| `plot_clusters(fit)` | what the clusters capture: species composition, shares over time, map, hook depth, clusters in PC space |
| `plot_indices()`, `plot_diagnostics()` | standardised indices and their accuracy, for one or several scenarios |

`bind_comparisons()` and `select_models()` combine and subset standardisation
results, e.g. to attach a slow control fit without refitting it.
`run_replicates()` repeats simulation and standardisation over many seeds,
saving after every fit so that long runs can be resumed; `plot_replicates()`
shows the accuracy across seeds.

## Simulation scripts and pre-run vignettes

Longer runs live in `inst/scripts/` (run from the package folder):

| Script | Runs |
|---|---|
| `run_bathymetric.R` | bathymetric example with the guideline settings: single run (seed 123) and 20 replicates |
| `run_tuna_flat.R` | tuna, flat true trends: single run and 20 replicates, fixed and moving thermocline |
| `run_tuna_contrast.R` | tuna, contrasting true trends: as above |
| `render_docs.R` | renders the vignettes to `docs/` from the pre-run results |

Shared settings (prices, models, seeds, number of parallel workers) are in
`talk_setup.R`. Each script first saves the single run (seed 123) used by the
vignettes and then runs the replicates; `options(cpueTargeting.reps = FALSE)`
runs only the single runs. Results are saved to `rdata/` in the repository,
which is ignored by git, so they are not pushed (`options(cpueTargeting.rdata
= "...")` saves them elsewhere). Finished seeds are skipped, so runs can be
resumed. The vignettes *Targeting along a continuous gradient: the
bathymetric example* and *Tuna longline targeting: flat versus contrasting
true trends* load these results (or compute them with the same code when
`options(cpueTargeting.run_vignettes = TRUE)`); they are also on the [project page](https://henning-winker.github.io/cpueTargeting/).

To change a scenario, edit the settings list:

```r
s <- tuna_settings()
s$depth_sd <- c(80, 80, 80, 80)      # more vertical mixing of species
s$HBF_choice_temp <- 0.05            # skippers choose gear more strictly
sim <- simulate_tuna(s, seed = 1)

# Extreme targeting (price contrast 2.0 vs 0.5, strict skippers)
stress <- simulate_tuna(tuna_settings("stress_test"), seed = 1)
```

## Speed

A default 30-year run takes a few seconds. The original code was slow because
it built each species' 2500 x 2500 movement matrix cell by cell, iterated the
movement 1000 times for the initial distribution, and grew data frames with
`rbind()` inside the cell loop. The package uses the closed-form stable
distribution, never forms the movement matrix explicitly, and simulates all
sets of a year at once.

## Original code

The original scripts used for the guidelines ("Run 22") (`OMcpueSIM22_original.R`,
`Generate_scenario_data22_original.R`) and the 2014 cluster simulator are kept
unchanged in `inst/extdata` for reference. The differences between the
original code and `simulate_bathymetric()` are listed in Appendix A.8 of the
[background document](https://henning-winker.github.io/cpueTargeting/appendix.html)
and in `NEWS.md`.

## References

He, X., Bigelow, K.A. and Boggs, C.H. (1997) Cluster analysis of longline sets
and fishing strategies within the Hawaii-based fishery. Fisheries Research 31,
147-158.

Hoyle, S.D., Langley, A.D. and Campbell, R.A. (2014) Recommended approaches for
standardizing CPUE data from pelagic fisheries. WCPFC-SC10-2014/SA-IP-10.

Ono, K., Punt, A.E. and Hilborn, R. (2015) Think outside the grids: an
objective approach to define spatial strata for catch and effort analysis.
Fisheries Research 170, 89-101.

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
