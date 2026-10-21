# run from the repository root
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
dir.create("dev/mmsy/output", showWarnings = FALSE)
L <- derive_landscape(); sys <- ms_system(L)
lv <- c(0, 10^seq(0.5, 3.5, length.out = 20))
grid <- as.matrix(expand.grid(demersal = lv, pelagic = lv, mixed = lv))
out <- list()
for (b in c("fixed", "dynamic")) {
  t0 <- Sys.time()
  x <- ms_simulate(sys, grid, years = 150, gamma = 1.5, behaviour = b)
  s <- ms_summary(x); BB <- attr(s, "BBmsy"); s <- cbind(as.data.frame(grid), s)
  # segment profits at equilibrium (approx: profit share by segment effort is not tracked; use effort share)
  cat("\n==", b, nrow(grid), "scenarios | time", format(round(Sys.time() - t0, 1)), "\n")
  out[[b]] <- list(s = s, BBmsy = BB, C = attr(ms_summary(x), "catch"), x_Eseg = NULL)
}
saveRDS(out, "dev/mmsy/output/grid_segments.rds")
