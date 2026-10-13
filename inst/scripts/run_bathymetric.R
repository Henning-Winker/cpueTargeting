# Bathymetric (guideline Run 22) example --------------------------------------------
# Settings and models as in the CPUE guidelines: 2 PCs, 3 clusters, spatial
# smoother s(Xs, Ys, k = 8) on 5 x 5 blocks, Tweedie GAMs; the final DPC model
# is fitted as a GAMM with a random vessel effect.

source("inst/scripts/talk_setup.R")
scenarios <- list(Bathymetric = bathymetric_settings())

# Single run ---------------------------------------------------------------------------
single_file <- file.path(rdata_dir, "bathymetric_single.rds")
if (!file.exists(single_file)) {
  sim <- simulate_bathymetric(seed = single_seed, plots = FALSE)
  fit <- do.call(standardize_cpue, c(list(sim, models = bathy_models, engine = "gam"),
                                     bathy_args))
  fit_re <- do.call(standardize_cpue,
                    c(list(sim, models = "DPC+S", engine = "gam", vessel_re = TRUE,
                           label = " + vessel RE"), bathy_args))
  saveRDS(list(sim = sim, fit = bind_comparisons(fit, fit_re)), single_file)
}

# Replicates (fixed-effect models only, bam for speed) --------------------------------
if (run_reps) {
  out_dir <- file.path(rdata_dir, "replicates_bathymetric")
  run_parallel(scenarios, out_dir, bathy_models, extra = NULL,
               args = c(bathy_args, engine = "bam"), simulator = simulate_bathymetric)
  reps <- read_replicates(out_dir)
  saveRDS(reps, file.path(rdata_dir, "bathymetric_replicates.rds"))
  plot_replicates(reps)
}
