# cpueTargeting 0.4.0

## Toolbox for real data (steps 1-6)

* `cpue_data()` brings long or wide logbook data into one row per operation
  with species columns, merges species synonyms, removes missing effort, and
  keeps a processing log; `subset()` method and `cpue_matrix()`.
* Exploration: `plot_fleet()` (operations by year, vessel x year tiles, fleet
  turnover, covariate shares), `select_vessels()` (core fleet),
  `nominal_cpue()`, `plot_nominal()`, `plot_composition()`, `plot_cpue_map()`.
* Species selection: `species_curve()`, `plot_species_curve()` and
  `select_species()` (cumulative share of occurrence or weight, optional
  `Other` group).
* Habitat-based subsetting of Stephens and MacCall (2004):
  `stephens_maccall()` and `plot_stephens_maccall()`.
* Targeting covariates: `composition_pca()`, `compute_oc()` (Optimal
  Coordinates test of the CPUE guidelines, robust when all components pass),
  `cluster_tactics()` (CLARA, k-means or Ward), `add_targeting()`,
  `tactic_table()`, `plot_oc()`, `plot_biplot()`, `plot_targeting()`.
* Models for real data (steps 7-9): `fit_cpue()` fits Base, Covariates, DPC
  and Cluster GAMs (optionally `+S` spatial and random vessel effects) with
  `mgcv`, Tweedie (power estimated jointly, no two-step `gamm()`) or
  delta-lognormal / delta-gamma families, and derives indices averaged over a
  reference set with delta-method log SEs. `cpue_index()`, `export_index()`.
* Diagnostics: `cpue_diagnostics()`, `plot_index()`, `plot_step()`,
  `influence_cpue()` and `plot_influence()` (Bentley et al. 2012),
  `plot_residuals()` (randomised quantile residuals for Tweedie),
  `plot_effects()`, `plot_vessel_effects()`.
* New vignette *Analysing real catch and effort data: a step-by-step guide*.
* GitHub Actions: R CMD check on every push; generated documentation is
  committed automatically.
* `cluster` added to Suggests.

# cpueTargeting 0.3.0

* `standardize_tuna_indices()` is renamed `standardize_cpue()`, as it works for
  both simulators. The old name is kept as an alias.
* New `plot_clusters()`: diagnostics for the Cluster models (composition of
  each cluster, cluster shares over time, most common cluster in space, hook
  depth by cluster, clusters in PC space), as in the CPUE guidelines.
* `run_replicates()` gains `simulator`, so it also runs
  `simulate_bathymetric()`; new `read_replicates()` combines results saved
  one file per seed (e.g. from parallel runs).
* New vignettes *Targeting along a continuous gradient: the bathymetric
  example* and *Tuna longline targeting: flat versus contrasting true
  trends*. They load pre-run results; the simulation scripts are in
  `inst/scripts/` (`run_bathymetric.R`, `run_tuna_flat.R`,
  `run_tuna_contrast.R`, with shared settings in `talk_setup.R`), and
  `render_docs.R` renders the vignettes to `docs/`. Results are saved to `rdata/`
  (ignored by git); `options(cpueTargeting.reps = FALSE)` runs only the
  single runs used by the vignettes.
* The 1 x 1 spatial control is no longer part of the standard comparisons:
  positions are reported as 5 x 5 blocks, the realistic, data-limited case.
* Single-run illustrations use seed 123.
* New setting `q_index` (default 0.2) separates the catchability at which
  CPUE is recorded from the fishing catchability `q`. Set-level Tweedie noise
  is drawn at `q_index`, and the catch removed scales the same draws by
  `q / q_index`. Fishing harder (larger `q`) therefore no longer removes the
  zeros and reduces the noise of the CPUE data, and avoids saturating
  `1 - exp(-q A)`, which damped the targeting signal. With `q_index = q` results are
  as before.

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
* Deep species are deeper and narrower by default (BET 260 m, SD 50; BFT
  240 m, SD 55), so deep gear clearly targets them while most sets still
  catch all four species.
* Default scenario: prices favour ALB (years 1-10), BET (11-20), then BET
  and BFT (21-30), with strong targeting (target species worth 4 times the
  others, `preference = 1.5`, `HBF_choice_temp = 0.15`).
* Stocks are fished harder (`q` = 3) and start partly fished: new
  per-species setting `B_start` (initial biomass / carrying capacity;
  default 0.6, 0.7, 0.2, 0.9). With `r` = 0.30, 0.18, 0.25, 0.25 the true
  trends differ: ALB is fished down, BET declines, BFT recovers and declines
  again, YFT dips and recovers.
* `effort` can be one value per year. By default it increases towards 800
  sets per year and is cut by 25% in years 21-30.
* Presets: `tuna_settings("mild")` (moderate targeting) and
  `tuna_settings("stress_test")` (strong targeting of lightly fished,
  unfished stocks, as in the first runs).

## New functions

* `plot_om()`, `plot_indices()`, `plot_diagnostics()`, `plot_thermocline()`,
  `plot_thermocline_effect()`, `plot_availability()`, `plot_pcs()`, `plot_sections()`,
  `bind_comparisons()`, `select_models()`, `run_replicates()`,
  `plot_replicates()` and `as_set_data()`.
