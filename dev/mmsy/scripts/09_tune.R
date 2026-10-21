# run from the repository root: tune daily costs so that the target fleet composition
# is the open-access equilibrium
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
ab <- c("DEC","PAC","SBR","HKE","COE","GFB","ALB","DOL","AMB","VMA")
L <- derive_landscape(); sys <- ms_system(L)
for (st in list(c(median = 0.9), c(min = 0.4))) {
  t0 <- Sys.time()
  tu <- ms_tune(sys, shares = c(demersal = 0.35, pelagic = 0.35, mixed = 0.30), status = st, gamma = 1.5, verbose = FALSE)
  sq <- attr(tu, "status_quo"); h <- attr(tu, "tuning")
  cat(sprintf("\n== target %s = %.2f | converged %s in %d iterations, max error %.3f (%.0fs)\n", names(st), st,
      attr(tu, "converged"), nrow(h), attr(tu, "max_error"), as.numeric(Sys.time() - t0, units = "secs")))
  cat("daily cost multipliers:", paste(sprintf("%s %.3f", names(tu$cost_seg), tu$cost_seg), collapse = ", "), "\n")
  cat(sprintf("open access: days %s (shares %s) | yield %.0f profit %.0f\n",
      paste(round(sq$E_seg), collapse = "/"), paste(round(sq$E_seg / sum(sq$E_seg), 3), collapse = "/"),
      sq$summary$yield, sq$summary$profit))
  cat("B/Bmsy:", paste(sprintf("%s %.2f", ab, sq$BBmsy), collapse = " "), "\n")
}
