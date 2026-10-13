# Tuna longline: contrasting true trends --------------------------------------------
# Single run (seed 123) for the illustrations and 20 replicates.
# Existing results are reused: copy last night's seed files (seed01.rds, ...)
# into rdata/replicates_contrast/ and they are skipped.

source("inst/scripts/talk_setup.R")
scenarios <- with_thermocline(tuna_contrast())

# Single run ---------------------------------------------------------------------------
single_file <- file.path(rdata_dir, "tuna_contrast_single.rds")
if (!file.exists(single_file)) {
  sims <- lapply(scenarios, simulate_tuna, seed = single_seed, plots = FALSE)
  fits <- lapply(sims, fit_tuna)
  saveRDS(list(sims = sims, fits = fits), single_file)
}

# Replicates ---------------------------------------------------------------------------
if (run_reps) {
  out_dir <- file.path(rdata_dir, "replicates_contrast")
  run_parallel(scenarios, out_dir, tuna_models, tuna_extra, tuna_args)
  reps <- read_replicates(out_dir)
  saveRDS(reps, file.path(rdata_dir, "tuna_contrast_replicates.rds"))
  plot_replicates(reps)
}
