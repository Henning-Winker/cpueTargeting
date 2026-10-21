# run from the repository root: stock status, F, effort and yield figures
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
library(patchwork)
figdir <- "dev/mmsy/figs"; dir.create(figdir, showWarnings = FALSE)
save_fig <- function(p, f, w = 10, h = 5) ggsave(file.path(figdir, f), p, width = w, height = h, dpi = 150, bg = "#fcfcfb")
L <- derive_landscape(); sys <- ms_system(L)

# ---- 1. Equilibrium against total effort (equal split) ----
E <- c(0, 10^seq(0, 3.6, by = 0.04)); caps <- cbind(E, E, E) / 3
runs_eq <- list(
  "Profit-blind (γ = 0)"       = ms_simulate(sys, caps, 150, gamma = 0,   behaviour = "dynamic"),
  "Frozen allocation (γ = 1.5)" = ms_simulate(sys, caps, 150, gamma = 1.5, behaviour = "fixed"),
  "Profit-seeking (γ = 1.5)"    = ms_simulate(sys, caps, 150, gamma = 1.5, behaviour = "dynamic"))
scans <- lapply(runs_eq, function(x) { s <- ms_summary(x); s$E <- E; s })
print(eq_limits(scans))
sub <- "Equilibrium with a third of total effort per segment (effort up to 600). Dashed verticals: highest effort with every stock ≥ Bpa."
save_fig(plot_eq_total(scans) + labs(title = "Yield, profit and the worst-off stock", subtitle = sub), "eq_total.png", 10, 3.8)
save_fig(plot_eq_species(scans, sys, "BBmsy") + labs(title = "Stock status at equilibrium",
  subtitle = "Lines: Bmsy (solid), Bpa = 0.5 Bmsy (dashed), Blim = 0.25 Bmsy (dotted)"), "eq_status.png", 10, 5)
save_fig(plot_eq_species(scans, sys, "FFmsy") + labs(title = "Fishing mortality at equilibrium",
  subtitle = "F relative to the single-species Fmsy = r/2"), "eq_F.png", 10, 5)
save_fig(plot_eq_species(scans, sys, "yield") + labs(title = "Yield by species at equilibrium",
  subtitle = "Free y-scales"), "eq_yield.png", 10, 5)
save_fig(plot_eq_effort(runs_eq, sys, E) + labs(title = "Where the effort goes",
  subtitle = "Share of total effort by métier at equilibrium"), "eq_effort_metier.png", 10, 3.8)

# ---- 2. Kobe at the constrained optimum of three cap structures (γ = 1.5) ----
lv <- c(0, 10^seq(0.5, 3.5, length.out = 20)); grid <- as.matrix(expand.grid(lv, lv, lv))
sg <- ms_summary(ms_simulate(sys, grid, 150, gamma = 1.5))
segmin <- apply(grid, 1, min) / pmax(rowSums(grid), 1e-9)
best <- function(ok) { i <- which(ok); i[which.max(sg$yield[i])] }
pick <- c("Equal split" = best(sg$n_below_Bpa == 0 & grid[,1] == grid[,2] & grid[,2] == grid[,3]),
          "Each segment ≥ 20%" = best(sg$n_below_Bpa == 0 & segmin >= 0.2),
          "Unrestricted caps" = best(sg$n_below_Bpa == 0))
pts <- do.call(rbind, lapply(names(pick), function(nm) { i <- pick[[nm]]
  data.frame(scenario = sprintf("%s\ncaps %s, yield %.0f", nm, paste(round(grid[i, ]), collapse = "/"), sg$yield[i]),
             species = sys$species, BBmsy = attr(sg, "BBmsy")[i, ], FFmsy = attr(sg, "FFmsy")[i, ]) }))
pts$scenario <- factor(pts$scenario, unique(pts$scenario))
save_fig(plot_kobe_eq(pts, sys) + labs(title = "Stock status at the best constrained cap structure",
  subtitle = "Maximum total yield with every stock ≥ Bpa. Caps: demersal/pelagic/mixed effort. Profit-seeking fleet (γ = 1.5)"),
  "kobe_caps.png", 10, 4.4)

# ---- 3. Projections with and without a segment HCR ----
own <- sapply(seq_len(sys$G), function(g) {
  acc <- if (sys$loyalty[g] > 0) sys$home[g, ] else rep(1, sys$M)
  c0 <- colSums(sys$Q * acc) * sys$K; c0 / sum(c0) >= 0.05 })
hcr_hs <- function(Blim = 0.25, Btrig = 0.5) function(st) {
  bb <- st$B / (st$K / 2)
  sapply(seq_len(ncol(own)), function(g) pmin(1, pmax(0, (apply(bb[, own[, g], drop = FALSE], 1, min) - Blim) / (Btrig - Blim)))) }
nsim <- 200; C2 <- matrix(c(19, 58, 19) * 2, nsim, 3, byrow = TRUE)
B0 <- matrix(sys$K * L$species$B_start, nsim, sys$S, byrow = TRUE)
runs_t <- list(
  "Fixed caps (2 × constrained optimum)" = ms_simulate(sys, C2, 50, gamma = 1.5, B_init = B0, proc_cv = 0.2, rho = 0.3, seed = 1, keep = "all"),
  "Same caps + segment hockey-stick HCR" = ms_simulate(sys, C2, 50, gamma = 1.5, B_init = B0, hcr = hcr_hs(), proc_cv = 0.2, rho = 0.3, seed = 1, keep = "all"))
sub2 <- paste0("200 simulations, process error CV 0.2 (AR1 0.3); median and 90% band.\n",
               "HCR: each segment's cap ramps from 0 at Blim to full at Bpa of its own most depleted stock.")
save_fig(plot_traj_species(runs_t, sys, "BBmsy") + labs(title = "Stock status over time", subtitle = sub2), "traj_status.png", 10, 5.3)
save_fig(plot_traj_species(runs_t, sys, "FFmsy") + labs(title = "Fishing mortality over time", subtitle = sub2), "traj_F.png", 10, 5.3)
save_fig(plot_traj_species(runs_t, sys, "yield") + labs(title = "Catch by species over time", subtitle = sub2), "traj_yield.png", 10, 5.3)
save_fig(plot_traj_fleet(runs_t, sys) + labs(title = "Effort by segment and total yield", subtitle = sub2), "traj_fleet.png", 10, 3.9)
cat("done\n")
