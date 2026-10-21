# run from the repository root: does lower reef catchability or value relax the Dentex choke?
# usage: Rscript dev/mmsy/scripts/10_reef_scenarios.R <scenario number>
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
ab <- c("DEC","PAC","SBR","HKE","COE","GFB","ALB","DOL","AMB","VMA")
reef_fish <- c("Dentex dentex", "Pagellus bogaraveo", "Seriola dumerili", "Conger conger")
scen <- list(
  "Baseline" = list(q_dec = 1, price_mult = 1),
  "Dentex reef q 0.6" = list(q_dec = 0.6, price_mult = 1),
  "Reef fish 20% cheaper" = list(q_dec = 1, price_mult = 0.8),
  "Both" = list(q_dec = 0.6, price_mult = 0.8))
k <- as.integer(commandArgs(TRUE)[1]); nm <- names(scen)[k]; sc <- scen[[k]]
st <- ms_settings(); st$tactics["Reef", "Dentex dentex"] <- sc$q_dec
L <- derive_landscape(st)
pr <- L$species$price; pr[L$species$species %in% reef_fish] <- pr[L$species$species %in% reef_fish] * sc$price_mult
sys <- ms_system(L, price = pr)
sys <- ms_tune(sys, status = c(min = 0.4), gamma = 1.5, verbose = FALSE)
sq <- attr(sys, "status_quo")
o1 <- ms_optimise(sys, "yield", levers = "total", upper = 3000, gamma = 1.5)
o2 <- ms_optimise(sys, "yield", min_share = 0.2, upper = 1000, gamma = 1.5)
row <- function(lab, E, s, BB) data.frame(scenario = nm, option = lab, days = paste(round(E), collapse = "/"),
  yield = round(s$yield), profit = round(s$profit), weakest = ab[which.min(BB)], min_BBmsy = round(min(BB), 2),
  DEC = round(BB[1], 2), reef_share = NA)
out <- rbind(row("Status quo", sq$E_seg, sq$summary, sq$BBmsy),
             row("MMSY one cap", o1$E_seg, o1$summary, o1$BBmsy),
             row("MMSY caps >= 20%", o2$E_seg, o2$summary, o2$BBmsy))
reef <- function(o) { Em <- colSums(o$E_seg_metier); round(Em["Reef_shallow"] / sum(Em), 2) }
xq <- ms_simulate(sys, matrix(1000, 1, 3), 300, gamma = 1.5, E_init = sq$E_seg, B_init = matrix(sq$B, 1))
Emq <- colMeans(apply(xq$E_metier, c(2, 3), mean)); out$reef_share <- c(round(Emq[1] / sum(Emq), 2), reef(o1), reef(o2))
out$costs <- paste(sprintf("%.2f", sys$cost_seg), collapse = "/")
saveRDS(out, sprintf("dev/mmsy/output/reef_scen_%d.rds", k)); print(out)
