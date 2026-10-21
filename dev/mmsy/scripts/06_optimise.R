# run from the repository root: constrained MMSY (all stocks >= Bpa) by lever set,
# and demersal + pelagic caps with the mixed segment as an uncontrolled free agent
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R")
dir.create("dev/mmsy/output", showWarnings = FALSE)
L <- derive_landscape(); sys <- ms_system(L)
ctl <- list(gamma = 1.5, phi = 0.5)       # profit-seeking fleet with economic exit
opt <- function(...) do.call(ms_optimise, c(list(sys, ...), ctl))

# ---- 1. Constrained MMSY and MMEY by lever set ----
lever_sets <- list("One global cap (equal split)" = list(levers = "total"),
                   "Demersal + pelagic caps, mixed = 0" = list(levers = c("demersal", "pelagic"), fixed = c(mixed = 0)),
                   "Three segment caps" = list(levers = c("demersal", "pelagic", "mixed")))
res1 <- list()
for (nm in names(lever_sets)) for (obj in c("yield", "revenue")) {
  t0 <- Sys.time(); o <- do.call(opt, c(list(objective = obj), lever_sets[[nm]]))
  res1[[paste(nm, obj, sep = " | ")]] <- o
  cat(sprintf("%-38s %-7s caps %-14s effort %-14s yield %4.0f revenue %5.0f profit %5.0f binding %s (%.0fs)\n",
      nm, obj, paste(round(o$caps), collapse = "/"), paste(round(o$E_seg), collapse = "/"),
      o$summary$yield, o$summary$revenue, o$summary$profit, paste(o$binding, collapse = ", "),
      as.numeric(Sys.time() - t0, units = "secs")))
}

# ---- 2. Free agent: mixed effort uncontrolled, optimise demersal + pelagic caps ----
dem_m <- grep("^(Reef|Sediment)", sys$metiers$name)
Mlev <- c(0, 5, 10, 15, 20, 25, 30, 40, 60)
res2 <- do.call(rbind, lapply(Mlev, function(M) {
  o <- opt(objective = "yield", levers = c("demersal", "pelagic"), fixed = c(mixed = M))
  if (!isTRUE(o$feasible)) return(data.frame(mixed_cap = M, feasible = FALSE, demersal = NA, pelagic = NA,
    yield = NA, revenue = NA, mixed_effort = NA, mixed_in_demersal = NA, binding = NA))
  data.frame(mixed_cap = M, feasible = TRUE, demersal = o$caps[["demersal"]], pelagic = o$caps[["pelagic"]],
             yield = o$summary$yield, revenue = o$summary$revenue, mixed_effort = o$E_seg[["mixed"]],
             mixed_in_demersal = sum(o$E_seg_metier["mixed", dem_m]) / max(sum(o$E_seg_metier["mixed", ]), 1e-9),
             binding = paste(o$binding, collapse = ", "))
}))
print(res2, digits = 3)

# ---- 3. Leakage: tighten the demersal cap with the mixed segment free at 20 ----
o20 <- opt(objective = "yield", levers = c("demersal", "pelagic"), fixed = c(mixed = 20))
dcap <- seq(0, 40, by = 1)
caps3 <- cbind(demersal = dcap, pelagic = o20$caps[["pelagic"]], mixed = 20)
x3 <- do.call(ms_simulate, c(list(sys, caps3, years = 150), ctl))
Eg <- apply(x3$E_seg_metier, c(2, 3, 4), mean)            # n x G x M
leak <- data.frame(demersal_cap = dcap,
  demersal_effort_dem = rowSums(Eg[, 1, dem_m]), mixed_effort_dem = rowSums(Eg[, 3, dem_m]),
  mixed_effort_pel = rowSums(Eg[, 3, -dem_m]),
  DEC = attr(ms_summary(x3), "BBmsy")[, "Dentex dentex"], yield = ms_summary(x3)$yield)
print(leak[seq(1, nrow(leak), 5), ], digits = 3)
saveRDS(list(res1 = res1, res2 = res2, leak = leak, o20 = o20), "dev/mmsy/output/optimise.rds")
