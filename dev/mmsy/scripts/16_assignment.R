# run from the repository root: responsibility by catchability vs nominal assignment of stocks to fleets
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
sys <- readRDS("dev/mmsy/output/tuned_sys.rds"); sq <- attr(sys, "status_quo")
dem <- sys$species[sys$group == "demersal"]; pel <- sys$species[sys$group == "pelagic"]
assign_nom <- list(demersal = dem, pelagic = pel)                         # mixed segment: not managed
assign_all <- list(demersal = dem, pelagic = pel, mixed = sys$species)    # mixed answers for all stocks
base <- list(sys = sys, E_ref = sq$E_seg, B_init = sq$B, gamma = 1.5, proc_cv = 0.1, rho = 0.3,
             impl_cv = 0.1, capacity = 1000, close = "Dentex dentex",
             triggers = seq(0.4, 1.0, by = 0.1), close_triggers = c(1, 1.5, 2, 3, 5), nsim = 100)
sc <- list("By catchability (all segments that catch a stock)" = list(),
           "Nominal: demersal->demersal, pelagic->pelagic; mixed unmanaged" = list(assign = assign_nom),
           "Nominal + mixed answers for all stocks" = list(assign = assign_all))
# ms_tune_hcr passes assign through ... only to ms_simulate, so build the rule grid here
run <- function(a) {
  g <- expand.grid(Btrig = base$triggers, Btrig_close = base$close_triggers); k <- nrow(g); n <- k * base$nsim
  rule <- hcr_effort(sys, Btrig = rep(g$Btrig, each = base$nsim), Btrig_close = rep(g$Btrig_close, each = base$nsim),
                     close = "Dentex dentex", assign = a$assign)
  x <- ms_simulate(sys, matrix(sq$E_seg, n, 3, byrow = TRUE), 50, gamma = 1.5, B_init = matrix(sq$B, n, sys$S, byrow = TRUE),
                   E_init = matrix(sq$E_seg, n, 3, byrow = TRUE), capacity = 1000, hcr = rule, crn = base$nsim, seed = 1,
                   keep = "all", proc_cv = 0.1, rho = 0.3, impl_cv = 0.1)
  tab <- cbind(g, do.call(rbind, lapply(seq_len(k), function(j) perf_rows(x, rows = (j - 1) * base$nsim + seq_len(base$nsim), years = 31:50))))
  ok <- tab$worst_P_below_Bpa <= 0.05
  b <- if (any(ok)) tab[ok, ][which.max(tab$yield[ok]), ] else tab[which.min(tab$worst_P_below_Bpa), ]
  cbind(feasible = any(ok), b)
}
out <- do.call(rbind, lapply(names(sc), function(nm) cbind(responsibility = nm, run(sc[[nm]]))))
print(out[, c("responsibility", "feasible", "Btrig", "Btrig_close", "yield", "profit", "days", "worst_P_below_Bpa", "worst_stock")], row.names = FALSE)
saveRDS(out, "dev/mmsy/output/assignment.rds")
