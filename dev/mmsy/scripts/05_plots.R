# run from the repository root: stock status, F, effort and yield figures
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
library(patchwork)
figdir <- "dev/mmsy/figs"; dir.create(figdir, showWarnings = FALSE)
save_fig <- function(p, f, w = 10, h = 5) ggsave(file.path(figdir, f), p, width = w, height = h, dpi = 150, bg = "#fcfcfb")
L <- derive_landscape(); sys <- ms_system(L)
# tune daily costs: 35/35/30% of active days, weakest stock at 0.4 Bmsy under open access
capacity <- 1000                                     # registered days per segment
sys <- ms_tune(sys, shares = c(demersal = 0.35, pelagic = 0.35, mixed = 0.30),
               status = c(min = 0.4), cap = capacity, gamma = 1.5, verbose = FALSE)

# ---- 1. Equilibrium against total effort (equal split) ----
E <- c(0, 10^seq(1, 4.2, by = 0.04)); caps <- cbind(E, E, E) / 3
runs_eq <- list(
  "Profit-blind (γ = 0)"       = ms_simulate(sys, caps, 150, gamma = 0,   behaviour = "dynamic"),
  "Frozen allocation (γ = 1.5)" = ms_simulate(sys, caps, 150, gamma = 1.5, behaviour = "fixed"),
  "Profit-seeking (γ = 1.5)"    = ms_simulate(sys, caps, 150, gamma = 1.5, behaviour = "dynamic"))
scans <- lapply(runs_eq, function(x) { s <- ms_summary(x); s$E <- E; s })
print(eq_limits(scans))
sub <- "Equilibrium with a third of total effort per segment. Dashed verticals: highest effort with every stock ≥ Bpa."
save_fig(plot_eq_total(scans, E_max = Inf) + labs(title = "Yield, profit and the worst-off stock", subtitle = sub), "eq_total.png", 10, 3.8)
save_fig(plot_eq_species(scans, sys, "BBmsy") + labs(title = "Stock status at equilibrium",
  subtitle = "Lines: Bmsy (solid), Bpa = 0.5 Bmsy (dashed), Blim = 0.25 Bmsy (dotted)"), "eq_status.png", 10, 5)
save_fig(plot_eq_species(scans, sys, "FFmsy") + labs(title = "Fishing mortality at equilibrium",
  subtitle = "F relative to the single-species Fmsy = r/2"), "eq_F.png", 10, 5)
save_fig(plot_eq_species(scans, sys, "yield") + labs(title = "Yield by species at equilibrium",
  subtitle = "Free y-scales"), "eq_yield.png", 10, 5)
save_fig(plot_eq_effort(runs_eq, sys, E) + labs(title = "Where the effort goes",
  subtitle = "Share of total effort by métier at equilibrium"), "eq_effort_metier.png", 10, 3.8)

# ---- 2. Kobe: status quo (capacity-limited) and safeguarded optima (γ = 1.5) ----
oa <- attr(sys, "status_quo")
o_y  <- ms_optimise(sys, "yield",  levers = "total", upper = 3 * capacity, gamma = 1.5)
o_p  <- ms_optimise(sys, "profit", levers = "total", upper = 3 * capacity, gamma = 1.5)
o_s  <- ms_optimise(sys, "yield",  min_share = 0.2, upper = capacity, gamma = 1.5)
lab <- function(nm, E, s) sprintf("%s\neffort %s\nyield %.0f, profit %.0f", nm, paste(round(E), collapse = "/"), s$yield, s$profit)
pts <- rbind(
  data.frame(scenario = lab("Status quo (tuned open access)", oa$E_seg, oa$summary), species = sys$species, BBmsy = oa$BBmsy, FFmsy = oa$FFmsy),
  data.frame(scenario = lab("Safeguarded MMSY, one cap", o_y$E_seg, o_y$summary), species = sys$species, BBmsy = o_y$BBmsy, FFmsy = o_y$FFmsy),
  data.frame(scenario = lab("Safeguarded MMSY, caps ≥ 20% each", o_s$E_seg, o_s$summary), species = sys$species, BBmsy = o_s$BBmsy, FFmsy = o_s$FFmsy),
  data.frame(scenario = lab("MMEY, one cap", o_p$E_seg, o_p$summary), species = sys$species, BBmsy = o_p$BBmsy, FFmsy = o_p$FFmsy))
pts$scenario <- factor(pts$scenario, unique(pts$scenario))
save_fig(plot_kobe_eq(pts, sys) + facet_wrap(~ scenario, nrow = 1) +
  labs(title = "Stock status: status quo against the safeguarded optima",
       subtitle = "Days are demersal/pelagic/mixed at equilibrium; safeguard: every stock ≥ Bpa. Status quo tuned to 35/35/30% of days, weakest stock 0.4 Bmsy"),
  "kobe_caps.png", 12, 4.8)

# ---- 3. Projections with and without a segment HCR ----
# a segment answers for the stocks it can put real pressure on: any metier it can use
# has at least 20% of that stock's highest catchability (not the stocks that dominate
# its catch: a by-catch stock that is a small share of every segment's catch would
# then be nobody's responsibility, and effort leaks to the unconstrained segments)
own <- sapply(seq_len(sys$G), function(g) {
  acc <- if (sys$loyalty[g] > 0) sys$home[g, ] > 0 else rep(TRUE, sys$M)
  apply(sys$Q[acc, , drop = FALSE], 2, max) >= 0.2 * apply(sys$Q, 2, max) })
dimnames(own) <- list(ms_codes[sys$species], names(sys$segments)); print(own)
hcr_hs <- function(Blim = 0.25, Btrig = 0.5) function(st) {
  bb <- st$B / (st$K / 2)
  sapply(seq_len(ncol(own)), function(g) pmin(1, pmax(0, (apply(bb[, own[, g], drop = FALSE], 1, min) - Blim) / (Btrig - Blim)))) }
nsim <- 200; B0 <- matrix(oa$B, nsim, sys$S, byrow = TRUE)
Csq <- matrix(capacity, nsim, 3)                              # registered fleet capacity
E0 <- matrix(oa$E_seg, nsim, 3, byrow = TRUE)
Cmm <- matrix(o_y$caps, nsim, 3, byrow = TRUE)
runs_t <- list(
  "No control (open access)" = ms_simulate(sys, Csq, 50, gamma = 1.5, B_init = B0, E_init = E0, proc_cv = 0.2, rho = 0.3, seed = 1, keep = "all"),
  "Status quo effort + segment HCR" = ms_simulate(sys, matrix(oa$E_seg, nsim, 3, byrow = TRUE), 50, gamma = 1.5, B_init = B0, E_init = E0, hcr = hcr_hs(), proc_cv = 0.2, rho = 0.3, seed = 1, keep = "all"),
  "Fixed cap at safeguarded MMSY" = ms_simulate(sys, Cmm, 50, gamma = 1.5, B_init = B0, E_init = E0, proc_cv = 0.2, rho = 0.3, seed = 1, keep = "all"))
sub2 <- paste0("From the tuned open-access status quo; 200 simulations, process error CV 0.2 (AR1 0.3); median and 90% band.\n",
               "HCR: each segment's status-quo days ramp from 0 at Blim to full at Bpa of the most depleted stock it can catch substantially.")
save_fig(plot_traj_species(runs_t, sys, "BBmsy") + labs(title = "Stock status over time", subtitle = sub2), "traj_status.png", 10, 5.3)
save_fig(plot_traj_species(runs_t, sys, "FFmsy") + labs(title = "Fishing mortality over time", subtitle = sub2), "traj_F.png", 10, 5.3)
save_fig(plot_traj_species(runs_t, sys, "yield") + labs(title = "Catch by species over time", subtitle = sub2), "traj_yield.png", 10, 5.3)
save_fig(plot_traj_fleet(runs_t, sys) + labs(title = "Effort by segment and total yield", subtitle = sub2), "traj_fleet.png", 10, 3.9)
cat("done\n")
