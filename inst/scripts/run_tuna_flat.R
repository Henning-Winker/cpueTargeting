# Tuna longline: flat true trends (same targeting as the contrast scenario) ----------
# Single run (seed 123) for the illustrations and 20 replicates.
# Existing results are reused: copy seed files from an earlier flat run
# into rdata/replicates_flat/ and they are skipped.

source("inst/scripts/talk_setup.R")
scenarios <- with_thermocline(tuna_flat())

# Single run ---------------------------------------------------------------------------
single_file <- file.path(rdata_dir, "tuna_flat_single.rds")
if (!file.exists(single_file)) {
  sims <- lapply(scenarios, simulate_tuna, seed = single_seed, plots = FALSE)
  fits <- lapply(sims, fit_tuna)
  saveRDS(list(sims = sims, fits = fits), single_file)
}

# Replicates ---------------------------------------------------------------------------
if (run_reps) {
  out_dir <- file.path(rdata_dir, "replicates_flat")
  run_parallel(scenarios, out_dir, tuna_models, tuna_extra, tuna_args)
  reps <- read_replicates(out_dir)
  saveRDS(reps, file.path(rdata_dir, "tuna_flat_replicates.rds"))
  plot_replicates(reps)
}
