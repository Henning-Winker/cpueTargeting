# cpueTargeting 0.2.0

## Tuna operating model rebuilt

* `simulate_tuna()` now has four species (YFT, BET, BFT, ALB) in a conceptual
  Southern Hemisphere ocean. Horizontal habitat is a preference for SST;
  vertical habitat is a normal depth profile shifted by a thermocline anomaly
  field, so the effect of HBF varies in space.
* Each HBF level has a typical hook depth (cf. Bigelow et al. 2002, Fig. 2)
  with lognormal set-level error (`set_cv`) and lognormal within-set spread
  (`hook_cv`). A set's catchability is vessel q times the overlap between hook
  and fish depth distributions (effective effort).
* Skippers choose cells and HBF from expected catch value.
* `Data` keeps the latent truth (`hook_depth`, `vertical_avail`) and the
  overlap a habitat model would predict from HBF and the average thermocline
  (`hbs_avail`).
* `standardize_tuna_indices()` follows the guideline workflow (Tweedie GAMs
  with estimated power). The default comparison covers the three practical
  approaches to targeting (Cluster, DPC, HBF) and their combinations
  (`Cluster_HBF`, `DPC_HBF`); habitat-based standardisation (`HBS`) is
  available on request. Optional random vessel effect. It also works for
  `simulate_bathymetric()` output.

## Faster and fewer dependencies

* The bathymetric and tuna models share one set of building blocks
  (`R/om_utils.R`). Movement matrices are never formed explicitly, the stable
  initial distribution uses its closed form, and all sets of a year are
  simulated at once. A 30-year run takes seconds.
* Removed dependencies on `fields`, `reshape2` and `tweedie`: the Gaussian
  random field (same covariance as `RandomFields::RMgauss`) and the Tweedie
  generator are now internal.
* `attach()` of the settings is no longer used. Missing settings are filled
  from the defaults.

## Corrections relative to the original Run 22 code

* Species growth rates `r` were recycled over cells instead of species
  (`r * B` with a cells x species matrix); `r` is now applied by species.
* The initial biomass `B0 * 0.5 * rlnorm(4, ...)` assumed four species; it now
  uses the number of species.
* The initial distribution is the exact stable distribution of the movement
  process (the original iterated 1000 times, which may not converge for short
  movement distances).
* Reported coarse positions are the centres of 5 x 5 blocks aligned with the
  grid (`obs_res`); the original rounding gave 11 uneven blocks.
* A closure set catchability to 1 in open cells; open cells now keep `qq_original`.
* Because random numbers are drawn in a different order, results do not match
  the original runs number by number.

## Moving thermocline

* New settings `thermocline_annual_sd`, `thermocline_annual_rho` and
  `thermocline_annual_field_sd` let the thermocline move between years
  (default 0: fixed). `Data$thermocline` is the anomaly of that year and the
  full field is returned in `Thermocline`.
* `standardize_tuna_indices()` gains user-defined models via `add_models`, and
  `engine = "bam"` for speed. Appending `+S` to any model name
  (e.g. `"DPC+S"`) adds the spatial smoother `s(Xs, Ys)`; `spatial_res` sets
  the reporting block size and `spatial_k` the basis size. The `Spatial`
  model is now `Base+S`.
* New vignette *Targeting in a longline tuna fishery*: fixed and moving
  thermocline scenarios, Cluster vs DPC vs HBF, what the PCs capture, and
  optional spatial-resolution and year-interaction steps.

## Standardisation helpers

* `n_pc` sets the number of principal components used by the DPC models
  (default 3; the guidelines used 2); all PCs and the PCA loadings are
  returned.
* `n_clusters` defaults to `n_pc + 1` (Winker et al. 2014), i.e. 4 clusters.
* `label`, `bind_comparisons()` and `select_models()` to combine and reuse
  fits (e.g. a 1 x 1 spatial control).
* Default `effort` in `tuna_settings()` reduced to 800 sets per year.
* Deep species are deeper and narrower by default (BET 260 m, SD 50; BFT
  240 m, SD 55), so deep gear clearly targets them while most sets still
  catch all four species.
* Default prices follow ALB (years 1-10), BET (11-20), BFT (21-30), with a
  moderate contrast (target species worth 1.5, others 0.75-1.0) and moderate
  skipper response (`preference = 1.25`, `HBF_choice_temp = 0.2`).
* More contrast in the true trends: `q` raised from 0.2 to 2, and new
  per-species setting `B_start` (initial biomass / carrying capacity; default
  0.7, 0.6, 0.5, 0.7). Stocks are fished down to roughly 20-50% of carrying
  capacity, fastest while they are targeted.
* `tuna_settings("stress_test")` keeps the earlier extreme targeting (price
  contrast 2.0 vs 0.5, `preference = 1.5`, `HBF_choice_temp = 0.15`).

## New functions

* `plot_om()`, `plot_indices()`, `plot_diagnostics()`, `plot_thermocline()`,
  `plot_thermocline_effect()`, `plot_availability()`, `plot_pcs()`, `plot_sections()`,
  `bind_comparisons()`, `select_models()`, `run_replicates()`,
  `plot_replicates()` and `as_set_data()`.
