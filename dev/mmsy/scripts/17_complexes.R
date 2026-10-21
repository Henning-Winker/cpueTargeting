# run from the repository root: assessments of species complexes (pooled taxa) vs single stocks
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
sys <- readRDS("dev/mmsy/output/tuned_sys.rds"); sq <- attr(sys, "status_quo")
sp <- setNames(sys$species, ms_codes[sys$species])
cx <- list("Reef fish" = unname(sp[c("DEC", "SBR", "COE", "AMB")]),
           "Shelf demersals" = unname(sp[c("PAC", "HKE", "GFB")]),
           "Large pelagics" = unname(sp[c("ALB", "DOL")]))
nsim <- 100; trig <- seq(0.4, 1.2, by = 0.1); ctrig <- c(0.8, 1, 1.5, 2, 3)
run <- function(complexes = NULL, index = "biomass", close = "Dentex dentex") {
  g <- expand.grid(Btrig = trig, Btrig_close = ctrig); k <- nrow(g); n <- k * nsim
  rule <- hcr_effort(sys, Btrig = rep(g$Btrig, each = nsim), Btrig_close = rep(g$Btrig_close, each = nsim),
                     close = close, complexes = complexes, complex_index = index)
  W <- sapply(attr(rule, "units"), function(m) (sys$species %in% m) * 1)       # pooled-biomass truth by unit
  x <- ms_simulate(sys, matrix(sq$E_seg, n, 3, byrow = TRUE), 50, gamma = 1.5, B_init = matrix(sq$B, n, sys$S, byrow = TRUE),
                   E_init = matrix(sq$E_seg, n, 3, byrow = TRUE), capacity = 1000, hcr = rule, crn = nsim, seed = 1,
                   keep = "all", proc_cv = 0.1, rho = 0.3, impl_cv = 0.1)
  yr <- 31:50
  tab <- cbind(g, do.call(rbind, lapply(seq_len(k), function(j) {
    rows <- (j - 1) * nsim + seq_len(nsim)
    Bu <- apply(x$B[yr, rows, , drop = FALSE], 1:2, function(b) b %*% W) / colSums(W * sys$K / 2)   # U x yr x sims
    p_unit <- apply(Bu < 0.5, 1, mean)
    cbind(perf_rows(x, rows = rows, years = yr), worst_P_unit = round(max(p_unit), 3), worst_unit = colnames(W)[which.max(p_unit)])
  })))
  ok <- tab$worst_P_unit <= 0.05                                  # what a manager can check
  b <- tab[ok, ][which.max(tab$yield[ok]), ]
  b
}
res <- list(
  "Single stocks (reference)" = run(),
  "Complexes, pooled biomass; closure on Reef fish" = run(cx, "biomass", close = "Reef fish"),
  "Complexes, pooled CPUE; closure on Reef fish" = run(cx, "cpue", close = "Reef fish"),
  "Complexes, pooled biomass; no closure" = run(cx, "biomass", close = NULL))
out <- do.call(rbind, lapply(names(res), function(nm) with(res[[nm]], data.frame(assessment = nm,
  triggers = if (exists("Btrig_close")) sprintf("%.1f / %s", Btrig, Btrig_close) else sprintf("%.1f", Btrig),
  yield, profit, days, risk_units_seen = worst_P_unit, true_worst_stock_risk = worst_P_below_Bpa,
  worst_stock = ms_codes[worst_stock], P_any_below_Blim))))
print(out, row.names = FALSE)
saveRDS(list(res = res, out = out), "dev/mmsy/output/complexes.rds")
