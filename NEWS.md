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
  with estimated power) and adds `HBF`, `DPC_HBF` and habitat-based (`HBS`)
  models; optional random vessel effect. It also works for
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
* New vignette *Vertical targeting and a moving thermocline*: constant, fixed
  and moving scenarios, standardisation, and year-interaction models.

## New functions

* `plot_om()`, `plot_indices()`, `plot_diagnostics()`, `plot_thermocline()`,
  `plot_thermocline_effect()`, `plot_availability()` and `as_set_data()`.
