# Shared settings for the simulation scripts --------------------------------------
#
# Sourced by run_bathymetric.R, run_tuna_flat.R and run_tuna_contrast.R.
# Run the scripts from the package folder (the git repository), e.g.
#   source("inst/scripts/run_tuna_flat.R")
# Results are saved to `rdata_dir` (default ../rdata, i.e. next to the
# repository), so they are not pushed to GitHub. Parallel runs use the
# INSTALLED package, so install the current version first:
#   devtools::install()   (or remotes::install_github(..., ref = "<branch>"))

library(cpueTargeting)

rdata_dir <- getOption("cpueTargeting.rdata", "../rdata")
dir.create(rdata_dir, showWarnings = FALSE, recursive = TRUE)

single_seed <- 123      # seed of the single-run illustrations
seeds       <- 1:20     # replicates
n_workers   <- 4        # parallel workers (~2 GB RAM each)

# Tuna: prices for the talk (columns YFT, BET, BFT, ALB) ----------------------------
talk_prices <- matrix(c(rep(c(0.3, 0.2, 1.0, 2.0), 10),   # years 1-10: ALB (+ BFT)
                        rep(c(1.5, 2.5, 0.7, 0.3), 10),   # years 11-20: BET (+ YFT)
                        rep(c(1.0, 2.0, 1.5, 0.5), 10)),  # years 21-30: BET / BFT
                      ncol = 4, byrow = TRUE)
switch_years <- c(10.5, 20.5)

# Contrasting true trends: package defaults (stocks start partly fished and are
# fished down or recover at different rates) with the talk prices
tuna_contrast <- function() {
  s <- tuna_settings()
  s$price <- talk_prices
  s
}
# Flat true trends: same targeting, light fishing of unfished stocks
tuna_flat <- function() {
  utils::modifyList(tuna_contrast(), list(q = 0.2, B_start = rep(1, 4)))
}
# Fixed and moving thermocline versions of a scenario
with_thermocline <- function(s) {
  list(Fixed = s,
       Moving = utils::modifyList(s, list(thermocline_annual_sd = 15,
                                          thermocline_annual_field_sd = 10)))
}

# Models: the same for all tuna scenarios -------------------------------------------
tuna_models <- c("Base", "Base+S", "Cluster+S", "DPC+S", "HBF+S")
tuna_extra  <- list(list(models = "DPC+S", n_pc = 2, label = " (2 PCs)", engine = "bam"))
tuna_args   <- list(spatial_res = 5, n_pc = 3, engine = "bam")   # 4 clusters (n_pc + 1)

# Bathymetric (guideline Run 22) example: as in the CPUE guidelines -------------------
bathy_models <- c("Base", "Base+S", "Cluster+S", "DPC", "DPC+S")
bathy_args   <- list(spatial_res = 5, spatial_k = 8, n_pc = 2, n_clusters = 3)

# Helpers ----------------------------------------------------------------------------

# Fit the models of a scenario to one simulation (single-run illustration)
fit_tuna <- function(sim) {
  fit <- do.call(standardize_cpue, c(list(sim, models = tuna_models), tuna_args))
  for (args in tuna_extra) {
    fit <- bind_comparisons(fit, do.call(standardize_cpue, c(list(sim), args)))
  }
  fit
}

# Run replicates in parallel, one file per seed; finished seeds are skipped
run_parallel <- function(scenarios, out_dir, models, extra, args,
                         simulator = simulate_tuna) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  run_seed <- function(seed) {
    do.call(cpueTargeting::run_replicates,
            c(list(scenarios, seeds = seed, models = models, extra = extra,
                   simulator = simulator,
                   file = file.path(out_dir, sprintf("seed%02d.rds", seed))),
              args))
    seed
  }
  todo <- seeds[!file.exists(file.path(out_dir, sprintf("seed%02d.rds", seeds)))]
  message(length(todo), " seeds to run in ", out_dir)
  if (!length(todo)) return(invisible(NULL))
  t0 <- Sys.time()
  if (n_workers > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    invisible(future.apply::future_lapply(todo, run_seed, future.seed = TRUE))
  } else {
    for (seed in todo) run_seed(seed)
  }
  print(Sys.time() - t0)
}
