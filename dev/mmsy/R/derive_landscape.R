# Uses internal helpers of cpueTargeting (install it, or devtools::load_all() from the repo root)
for (f in c("merge_logbook_settings", "simulate_grf", "habitat_density",
           "dispersal_kernels", "movement_operator", "logbook_settings"))
  assign(f, get(f, asNamespace("cpueTargeting")))
#' Logbook settings balanced for the fleet system: Seriola troll catchability
#' 0.8 -> 0.4 (with 0.8 the reef + troll by-catch collapses it under open access)
ms_settings <- function() {
  st <- logbook_settings()
  st$tactics["Pelagic", "Seriola dumerili"] <- 0.4
  st
}

derive_landscape <- function(settings = ms_settings(), seed = 1, depth_cut = 300, coastal_cut = 150) {
  s <- merge_logbook_settings(settings); set.seed(seed)
  sp <- s$species; ns <- nrow(sp)
  nx <- s$grid[1]; ny <- s$grid[2]; n_cells <- nx*ny
  x <- rep(seq_len(nx), ny); y <- rep(seq_len(ny), each = nx)
  depth <- s$shelf_depth[1] + diff(s$shelf_depth) * ((x - 1)/(nx - 1))^s$shelf_power +
    as.vector(simulate_grf(nx, ny, s$depth_sd, s$depth_scale))
  depth <- pmax(depth, 5)
  reef_field <- as.vector(simulate_grf(nx, ny, 1, s$reef_scale))
  reef <- reef_field >= stats::quantile(reef_field, 1 - s$reef_share)
  sst <- s$sst_range[2] - diff(s$sst_range) * (y - 1)/(ny - 1) +
    as.vector(simulate_grf(nx, ny, s$sst_sd, s$sst_scale))
  habitat <- vapply(seq_len(ns), function(i) {
    h <- habitat_density(depth, "Lognorm", sp$depth_mu[i], sp$depth_cv[i])
    h <- h * ifelse(reef, sp$reef[i], sp$sediment[i])
    if (!is.na(sp$sst_mu[i])) h <- h * stats::dnorm(sst, sp$sst_mu[i], sp$sst_sd[i])
    h / max(h)}, numeric(n_cells))
  kernels <- dispersal_kernels(x, y, sp$dispersal)
  carrying <- vapply(seq_len(ns), function(i) movement_operator(kernels[[i]], habitat[, i])$stable * sp$B0[i], numeric(n_cells))
  list(depth = depth, reef = reef, sst = sst, K = carrying, species = sp, settings = s)
}
