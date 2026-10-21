# run from the repository root
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
dir.create("dev/mmsy/output", showWarnings = FALSE)
L <- derive_landscape(); sys <- ms_system(L)
# 1. nesting: gamma = 0 fixed vs dynamic
E <- c(5, 20, 60); caps <- cbind(E, E, E)
a <- ms_summary(ms_simulate(sys, caps, 150, gamma = 0, behaviour = "fixed"))
b <- ms_summary(ms_simulate(sys, caps, 150, gamma = 0, behaviour = "dynamic"))
cat("gamma=0 max |diff| yield fixed vs dynamic:", max(abs(a$yield - b$yield)), "\n")
# 2. analytic: fixed allocation -> effective annual F_i; Schaefer equilibrium B/K = 1 - F/r (continuous approx)
x <- ms_simulate(sys, caps, 300, gamma = 0, behaviour = "fixed")
Em <- apply(x$E_metier, c(2,3), mean)            # n x M annual effort by metier
Fi <- t(sapply(1:3, function(k) rowSums(sapply(1:12, function(mo)
        (Em[k,] * sys$month_w[mo]) %*% sys$Q * sys$avail[mo,]))))
Bsim <- apply(x$B, c(2,3), mean) / matrix(sys$K, 3, sys$S, byrow = TRUE)
Bana <- pmax(1 - Fi / matrix(sys$r, 3, sys$S, byrow = TRUE), 0)
cat("B/K simulated vs analytic (E = 20 per segment):\n"); print(round(rbind(sim = Bsim[2,], analytic = Bana[2,]), 3))
# 3. gamma sensitivity, equal split and best segment split (each >= 20%)
lv <- c(0, 10^seq(0.5, 3.5, length.out = 20)); grid <- as.matrix(expand.grid(lv, lv, lv))
best <- function(sysx, g) {
  s <- ms_summary(ms_simulate(sysx, grid, 150, gamma = g, behaviour = "dynamic"))
  segmin <- apply(grid, 1, min) / pmax(rowSums(grid), 1e-9)
  eq <- which(s$n_below_Bpa == 0 & grid[,1] == grid[,2] & grid[,2] == grid[,3]); i <- eq[which.max(s$yield[eq])]
  sg <- which(s$n_below_Bpa == 0 & segmin >= 0.2); j <- sg[which.max(s$yield[sg])]
  sprintf("equal split: yield %4.0f profit %5.0f (E/seg %4.0f) | segments>=20%%: yield %4.0f profit %5.0f caps %s",
          s$yield[i], s$profit[i], grid[i,1], s$yield[j], s$profit[j], paste(round(grid[j,]), collapse = "/"))
}
for (g in c(0, 0.5, 1.5, 3)) cat("gamma", g, ":", best(sys, g), "\n")
# 4. choke species becomes the most valuable: Dentex price 12 -> 30
sys2 <- ms_system(L, price = replace(L$species$price, 1, 30))
cat("\nDentex price 30:\n")
for (g in c(0, 1.5, 3)) cat("gamma", g, ":", best(sys2, g), "\n")
