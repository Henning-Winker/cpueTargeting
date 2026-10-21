# run from the repository root: management based on four assessed stocks
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
sys <- readRDS("dev/mmsy/output/tuned_sys.rds"); sq <- attr(sys, "status_quo")
sp <- setNames(sys$species, ms_codes[sys$species])
main4  <- unname(sp[c("HKE", "SBR", "DOL", "VMA")])           # main commercial stocks
choke4 <- unname(sp[c("DEC", "SBR", "HKE", "VMA")])           # including the choke stock
no_reef <- replace(rep(1, sys$M), sys$metiers$name == "Reef_shallow", 0)
base <- list(sys = sys, E_ref = sq$E_seg, B_init = sq$B, gamma = 1.5, proc_cv = 0.1, rho = 0.3,
             impl_cv = 0.1, capacity = 1000, risk = c(worst_P_assessed = 0.05))
grid2 <- list(triggers = seq(0.4, 1.0, by = 0.1), close_triggers = c(1, 1.5, 2, 3, 5), nsim = 100)
sc <- list(
  "All 10 assessed, effort + closure" = c(grid2, list(close = "Dentex dentex")),
  "4 main assessed (no DEC), effort only" = list(triggers = seq(0.5, 1.5, by = 0.05), assessed = main4),
  "4 incl. DEC assessed, effort + closure" = c(grid2, list(close = "Dentex dentex", assessed = choke4)),
  "4 main assessed + permanent reef closure" = list(triggers = seq(0.5, 1.5, by = 0.05), assessed = main4, open = no_reef))
res <- lapply(sc, function(a) do.call(ms_tune_hcr, c(base, a)))
out <- do.call(rbind, lapply(names(res), function(nm) { b <- res[[nm]]$best
  if (is.null(b)) return(data.frame(set_up = nm, triggers = "none feasible"))
  data.frame(set_up = nm, triggers = if (!is.null(b$Btrig_close)) sprintf("effort %.2f, closure %s", b$Btrig, b$Btrig_close) else sprintf("effort %.2f", b$Btrig),
             yield = b$yield, profit = b$profit, days = b$days, risk_assessed = b$worst_P_assessed,
             risk_all = b$worst_P_below_Bpa, worst_stock = ms_codes[b$worst_stock], P_any_below_Blim = b$P_any_below_Blim) }))
print(out, row.names = FALSE)
# per-stock risk at the tuned rule of each set-up
saveRDS(list(res = res, out = out), "dev/mmsy/output/assessed_subset.rds")
