# run from the repository root: export the tuned system as plain tables (input format for the new package)
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
sys <- readRDS("dev/mmsy/output/tuned_sys.rds")
L <- derive_landscape()
sp <- data.frame(species = sys$species, code = c("DEC","PAC","SBR","HKE","COE","GFB","ALB","DOL","AMB","VMA"),
                 group = sys$group, r = sys$r, K = sys$K, price = sys$price,
                 season_amp = L$species$season_amp, season_peak = L$species$season_peak)
met <- data.frame(metier = sys$metiers$name, tactic = sys$metiers$tactic, ground = sys$metiers$ground,
                  cost = round(sys$cost, 6))
q <- data.frame(metier = rownames(sys$Q), signif(sys$Q, 6), check.names = FALSE)
seg <- data.frame(segment = names(sys$segments), share_target = c(0.35, 0.35, 0.30), loyalty = sys$loyalty,
                  efficiency = sys$eff[, 1], cost_multiplier = round(sys$cost_seg, 4), vessel_cv = sys$vessel_cv,
                  home = sapply(names(sys$segments), function(g) paste(colnames(sys$home)[sys$home[g, ] > 0], collapse = ";")))
write.csv(sp, "dev/mmsy/extdata/atlantic_species.csv", row.names = FALSE)
write.csv(met, "dev/mmsy/extdata/atlantic_metiers.csv", row.names = FALSE)
write.csv(q, "dev/mmsy/extdata/atlantic_catchability.csv", row.names = FALSE)
write.csv(seg, "dev/mmsy/extdata/atlantic_segments.csv", row.names = FALSE)
print(seg)
