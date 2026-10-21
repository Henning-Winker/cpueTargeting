# run from the repository root
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
dir.create("dev/mmsy/output", showWarnings = FALSE)
L <- derive_landscape(); sys <- ms_system(L)
Bpa <- 0.5
pick <- function(s, crit, cons = TRUE) { ok <- if (cons) s$n_below_Bpa == 0 else TRUE; s[ok, ][which.max(s[ok, crit]), ] }
Etot <- c(0, 10^seq(0, 4, by = 0.025)); caps1 <- cbind(Etot, Etot, Etot) / 3
res <- list()
for (b in c("fixed", "dynamic")) {
  t0 <- Sys.time(); x <- ms_simulate(sys, caps1, years = 150, gamma = 1.5, behaviour = b)
  s <- ms_summary(x); s$E <- Etot; res[[b]] <- s
  cat("\n==", b, "| time", format(round(Sys.time() - t0, 2)), "\n")
  print(rbind(MMSY = s[which.max(s$yield), ], MMEY = s[which.max(s$profit), ],
              `max yield, all>=Bpa` = pick(s, "yield"), `max profit, all>=Bpa` = pick(s, "profit"))[,
        c("E","yield","revenue","profit","min_BBmsy","worst","n_below_Bpa","n_below_Blim")], digits = 3)
}
saveRDS(res, "dev/mmsy/output/scan_global.rds")
