# run from the repository root
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
dir.create("dev/mmsy/output", showWarnings = FALSE)
L <- derive_landscape(); sys <- ms_system(L)
# segment-specific hockey-stick: each segment responds to the most depleted species
# among those that make up >= 5% of its unfished catch (its "own" choke species)
seg_species <- function(sys, thr = 0.05) {
  sapply(seq_len(sys$G), function(g) {
    acc <- if (sys$loyalty[g] > 0) sys$home[g, ] else rep(1, sys$M)
    c0 <- colSums(sys$Q * acc) * sys$K; c0 / sum(c0) >= thr })
}
own <- seg_species(sys); dimnames(own) <- list(sys$species, names(sys$segments)); print(own)
hcr_hs <- function(Blim = 0.25, Btrig = 0.5) function(st) {
  bb <- st$B / (st$K / 2)
  sapply(seq_len(ncol(own)), function(g) {
    d <- apply(bb[, own[, g], drop = FALSE], 1, min)
    pmin(1, pmax(0, (d - Blim) / (Btrig - Blim))) })
}
caps <- matrix(c(19, 58, 19) * 2, 1)            # caps at twice the constrained optimum
nsim <- 200; C <- caps[rep(1, nsim), ]
for (h in list(NULL, hcr_hs())) {
  t0 <- Sys.time()
  x <- ms_simulate(sys, C, years = 50, gamma = 1.5, hcr = h, proc_cv = 0.2, rho = 0.3,
                   seed = 1, keep = "all")
  bb <- sweep(x$B, 3, sys$K / 2, "/")
  pBlim <- apply(bb[41:50, , ] < 0.25, 3, mean)
  cat(if (is.null(h)) "no HCR" else "segment HS", "| 200 sims x 50 yr:", format(round(Sys.time() - t0, 2)),
      "| mean yield yrs 41-50:", round(mean(apply(x$C[41:50, , ], 1:2, sum))),
      "| max P(B<Blim):", round(max(pBlim), 2), names(which.max(pBlim)), "\n")
}
