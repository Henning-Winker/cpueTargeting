# run from the repository root: feedback control from the tuned open-access status quo
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
dir.create("dev/mmsy/output", showWarnings = FALSE); dir.create("dev/mmsy/figs", showWarnings = FALSE)
f_sys <- "dev/mmsy/output/tuned_sys.rds"
if (file.exists(f_sys)) sys <- readRDS(f_sys) else {
  sys <- ms_tune(ms_system(derive_landscape()), status = c(min = 0.4), cap = 1000, gamma = 1.5, verbose = FALSE)
  saveRDS(sys, f_sys)
}
sq <- attr(sys, "status_quo")
nsim <- 200; yrs <- 50
B0 <- matrix(sq$B, nsim, sys$S, byrow = TRUE); E0 <- matrix(sq$E_seg, nsim, sys$G, byrow = TRUE)
Esq <- matrix(sq$E_seg, nsim, sys$G, byrow = TRUE)                 # status-quo days: reference for the rules
run <- function(caps, hcr = NULL) ms_simulate(sys, caps, yrs, gamma = 1.5, B_init = B0, E_init = E0,
  capacity = matrix(1000, nsim, sys$G), hcr = hcr, proc_cv = 0.2, rho = 0.3, seed = 1, keep = "all")
t0 <- Sys.time()
runs <- list(
  "No control (open access)" = run(matrix(1000, nsim, sys$G)),
  "Effort rule, Btrig 0.8 Bmsy" = run(Esq, hcr_effort(sys, Btrig = 0.8)),
  "Effort rule + Dentex reef closure" = run(Esq, hcr_effort(sys, Btrig = 0.8, close = "Dentex dentex")))
extra <- list(
  "Effort rule, ±20%/yr" = run(Esq, hcr_effort(sys, Btrig = 0.8, max_change = 0.2)),
  "Effort rule + closure, ±20%/yr" = run(Esq, hcr_effort(sys, Btrig = 0.8, close = "Dentex dentex", max_change = 0.2)),
  "Effort rule + closure, no assessment error" = run(Esq, hcr_effort(sys, Btrig = 0.8, close = "Dentex dentex", assess_cv = 0)))
cat("runs:", format(round(Sys.time() - t0, 1)), "\n")
print(attr(hcr_effort(sys, close = "Dentex dentex"), "drive")[, "Dentex dentex"])
perf <- ms_performance(c(runs, extra)); print(perf, row.names = FALSE)
cat("\nStatus quo: yield", round(sq$summary$yield), "profit", round(sq$summary$profit), "\n")
# median B/Bmsy by species in years 31-50
bbm <- sapply(c(runs, extra), function(x) round(apply(sweep(x$B[31:50, , ], 3, sys$K / 2, "/"), 3, median), 2))
rownames(bbm) <- ms_codes[sys$species]; print(t(bbm))
saveRDS(list(perf = perf, bbm = bbm), "dev/mmsy/output/feedback.rds")
sub <- paste0("From the tuned open-access status quo; 200 simulations, process error CV 0.2 (AR1 0.3), assessment error CV 0.2.\n",
              "Effort rule: each segment's status-quo days ramp from 0 at Blim to full at 0.8 Bmsy of its most depleted stock.")
sv <- function(p, f, h = 5.3) ggsave(file.path("dev/mmsy/figs", f), p, width = 10, height = h, dpi = 150, bg = "#fcfcfb")
sv(plot_traj_species(runs, sys, "BBmsy") + labs(title = "Feedback control: stock status", subtitle = sub), "fb_status.png")
sv(plot_traj_species(runs, sys, "FFmsy") + labs(title = "Feedback control: fishing mortality", subtitle = sub), "fb_F.png")
sv(plot_traj_fleet(runs, sys) + labs(title = "Feedback control: days by segment and total yield", subtitle = sub), "fb_fleet.png", 3.9)
