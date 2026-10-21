# run from the repository root: constrained MMSY with a viability floor and a reef closure
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
L <- derive_landscape(); sys <- ms_system(L)
open_all <- rep(1, sys$M); no_reef_sh <- replace(open_all, sys$metiers$name == "Reef_shallow", 0)
runs <- list(
  "Equal split (one global cap)"            = list(levers = "total"),
  "Segment caps, each >= 20%"               = list(min_share = 0.2),
  "Segment caps, each >= 20%, shallow reef closed" = list(min_share = 0.2, open = no_reef_sh),
  "Equal split, shallow reef closed"        = list(levers = "total", open = no_reef_sh))
res <- list()
for (nm in names(runs)) {
  t0 <- Sys.time()
  o <- do.call(ms_optimise, c(list(sys, objective = "yield", gamma = 1.5, phi = 0.5), runs[[nm]]))
  res[[nm]] <- o
  cat(sprintf("%-48s caps %-12s yield %4.0f revenue %5.0f profit %5.0f binding %-28s (%.0fs)\n", nm,
      paste(round(o$caps), collapse = "/"), o$summary$yield, o$summary$revenue, o$summary$profit,
      paste(o$binding, collapse = ", "), as.numeric(Sys.time() - t0, units = "secs")))
  cat("   B/Bmsy:", paste(sprintf("%s %.2f", substr(names(o$BBmsy), 1, 8), o$BBmsy), collapse = " | "), "\n")
}
saveRDS(res, "dev/mmsy/output/viability_closure.rds")
