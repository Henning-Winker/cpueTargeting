# run from the repository root after 02_segment_grid.R
out <- readRDS("dev/mmsy/output/grid_segments.rds")
for (b in names(out)) {
  s <- out[[b]]$s; BB <- out[[b]]$BBmsy
  segmin <- pmin(s$demersal, s$pelagic, s$mixed) / pmax(s$E_total, 1e-9)
  rules <- list(`max yield, all>=Bpa` = s$n_below_Bpa == 0,
                `.. + each segment >=20% effort` = s$n_below_Bpa == 0 & segmin >= 0.2,
                `.. + equal 1/3 split` = s$n_below_Bpa == 0 & s$demersal == s$pelagic & s$pelagic == s$mixed)
  cat("\n==", b, "\n")
  for (nm in names(rules)) { ok <- which(rules[[nm]]); i <- ok[which.max(s$yield[ok])]
    cat(sprintf("%-32s dem %6.0f pel %6.0f mix %6.0f | yield %5.0f profit %5.0f | worst %s %.2f\n", nm,
      s$demersal[i], s$pelagic[i], s$mixed[i], s$yield[i], s$profit[i], s$worst[i], s$min_BBmsy[i])) }
}
