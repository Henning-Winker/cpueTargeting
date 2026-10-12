#' Default settings for the four-species bathymetric scenario (Run 22)
#'
#' Settings of the original `OMcpueSIM22.R` example (scenario 2, targeting
#' changes). Species distributions are driven by bathymetry and distance; the
#' model has no vertical component.
#'
#' Prices change in three 10-year blocks: years 1-10 favour species 1, years
#' 11-20 species 2, and years 21-30 species 4.
#'
#' @return A named list accepted by [simulate_bathymetric()].
#' @export
bathymetric_settings <- function() {
  list(
    n_years = 30,
    # Species (petrale, dover, darkblotched, sablefish in the original)
    n_species = 4,
    species = paste0("Sp", 1:4),
    price_fish = matrix(c(rep(c(1.5, 1, 1, 0.5), 10),
                          rep(c(0.5, 1.5, 1, 1), 10),
                          rep(c(0.5, 1, 0.5, 2), 10)), ncol = 4, byrow = TRUE),
    # Habitat preference and movement
    func_mvt_dist = rep("Exp", 4),        # distance kernel
    func_mvt_depth = rep("Lognorm", 4),   # depth preference: Norm, Exp, Lognorm, Unif
    func_mvt_lat = rep("Unif", 4),        # x-axis range preference
    Fish_dist_par1 = c(0.4, 0.4, 0.4, 1), # mean movement distance (cells)
    Fish_dist_par2 = rep(0, 4),           # unused (kept for compatibility)
    Fish_depth_par1 = c(220, 550, 250, 1000), # preferred depth (m)
    Fish_depth_par2 = c(0.5, 0.5, 0.3, 0.5),  # CV of depth preference
    Fish_range_par1 = rep(0, 4),
    Fish_range_par2 = rep(50, 4),
    # Schaefer dynamics
    B0 = c(320000, 420000, 56000, 540000),
    r = c(0.3, 0.17, 0.12, 0.08),
    no.mvt = FALSE,
    # Fleet
    Nvessels = 30,
    Tot_effort = 500,         # asymptotic number of sets per year
    CV_effort = 0.2,
    qq_original = 0.1,        # mean catchability per set
    CV_vessel = 0.2,          # between-vessel CV in catchability
    do.tweedie = TRUE,        # set-level Tweedie noise in catchability
    xi = 1.2,                 # Tweedie power
    phi = 0.1,                # Tweedie dispersion
    Preference = 1,           # effort concentration exponent
    Changing_preference = FALSE,
    # Space
    Range_X = 1:50,
    Range_Y = 1:50,
    SD_O = 150,               # SD of the bathymetry field (m)
    SpatialScale = 4,         # correlation scale of the bathymetry field (cells)
    obs_res = 5,              # reported positions are centres of obs_res x obs_res blocks
    year.depth.restriction = 1000, # first year of a closure (> n_years: none)
    depth.restriction = c(0, 1)    # closed depth range (m), or "random"
  )
}

#' Simulate the four-species bathymetric Run 22 scenario
#'
#' Four species are distributed over a simulated bathymetry according to depth
#' preferences and short-distance movement. Fishing effort is allocated in
#' proportion to expected catch value, which changes with species prices, so
#' targeting shifts over time.
#'
#' This is a re-implementation of the spatial multispecies simulator of
#' Thorson et al. (2017, Appendix A; function `Generate_scenario_data` in the
#' `spatial_DFA` package, <https://github.com/James-Thorson/spatial_DFA>),
#' which generalises the simulator of Ono et al. (2015), with the "Run 22"
#' settings used for the WKBMSYSPICT CPUE guidelines. A full description is
#' in `system.file("appendix", "cpueTargeting_simulator_appendix.html",
#' package = "cpueTargeting")`.
#'
#' Compared with the original script, the dynamics are identical except for
#' three corrections: species-specific `r` is now applied by species (the
#' original recycled `r` over cells), the initial distribution is the exact
#' stable distribution of the movement process, and reported block centres are
#' aligned with the grid.
#'
#' @param settings A named settings list, usually from [bathymetric_settings()].
#' @param seed Optional integer seed for reproducibility.
#' @param plots Logical; build the diagnostic ggplots in `Plots`.
#'
#' @return A list with
#' * `Data`: one row per set and species (`year`, `SetID`, `X`, `Y`, `Xobs`,
#'   `Yobs`, `Species`, `Sp`, `CPUE`, `depth`, `vessel`);
#' * `Biomass`: annual nominal and true indices, scaled to mean 1;
#' * `Biomass.XY`: initial density by cell and species;
#' * `Env`: the bathymetry;
#' * `Plots`: diagnostic ggplots (if `plots = TRUE`);
#' * `settings`: the settings used.
#' @references Thorson, J.T., Fonner, R., Haltuch, M.A., Ono, K. and Winker, H.
#'   (2017) Accounting for spatiotemporal variation and fisher targeting when
#'   estimating abundance from multispecies fishery data. Canadian Journal of
#'   Fisheries and Aquatic Sciences 74, 1794-1807.
#'
#'   Ono, K., Punt, A.E. and Hilborn, R. (2015) Think outside the grids: an
#'   objective approach to define spatial strata for catch and effort
#'   analysis. Fisheries Research 170, 89-101.
#' @export
#'
#' @examples
#' \dontrun{
#' om <- simulate_bathymetric(seed = 42)
#' head(om$Data)
#' om$Plots$CPUE
#' }
simulate_bathymetric <- function(settings = bathymetric_settings(), seed = NULL,
                                 plots = TRUE) {
  if (!is.list(settings) || is.null(names(settings))) {
    stop("`settings` must be a named list.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  s <- utils::modifyList(bathymetric_settings(), settings)
  ns <- s$n_species
  species <- if (length(s$species) == ns) s$species else paste0("Sp", seq_len(ns))

  # Grid and bathymetry (x varies fastest, as in the original)
  nx <- length(s$Range_X)
  ny <- length(s$Range_Y)
  n_cells <- nx * ny
  x <- rep(s$Range_X, ny) + 0.5
  y <- rep(s$Range_Y, each = nx) + 0.5
  z <- as.vector(simulate_grf(nx, ny, s$SD_O, s$SpatialScale))
  z <- z - min(z)                       # depth in m, positive downwards

  # Habitat suitability and movement
  habitat <- vapply(seq_len(ns), function(i) {
    habitat_density(z, s$func_mvt_depth[i], s$Fish_depth_par1[i],
                    s$Fish_depth_par2[i]) *
      habitat_density(x, s$func_mvt_lat[i], s$Fish_range_par1[i],
                      s$Fish_range_par2[i])
  }, numeric(n_cells))
  kernels <- dispersal_kernels(x, y, s$Fish_dist_par1)
  movement <- lapply(seq_len(ns), function(i) movement_operator(kernels[[i]], habitat[, i]))

  # Initial biomass = carrying capacity in each cell. As in the original,
  # stocks start at 0.5 * B0 (with lognormal noise) distributed at the stable
  # distribution, and cell-level carrying capacity equals this initial state.
  b_init <- s$B0 * 0.5 * stats::rlnorm(ns, 0, 0.2)
  biomass <- vapply(seq_len(ns), function(i) movement[[i]]$stable * b_init[i],
                    numeric(n_cells))
  carrying <- biomass

  # Fleet
  sig_vessel <- sqrt(log(s$CV_vessel^2 + 1))
  q_vessel <- stats::rlnorm(s$Nvessels, log(s$qq_original) - sig_vessel^2 / 2, sig_vessel)
  preference <- s$Preference

  records <- vector("list", s$n_years)
  total_biomass <- matrix(0, s$n_years, ns)
  set_counter <- 0L

  for (yr in seq_len(s$n_years)) {
    total_biomass[yr, ] <- colSums(biomass)
    if (isTRUE(s$Changing_preference)) preference <- 2 - exp(-0.1 * yr)
    n_sets <- annual_effort(yr, s$Tot_effort, s$CV_effort)

    # Expected catch value by cell drives effort allocation
    q_cell <- rep(s$qq_original, n_cells)
    if (yr >= s$year.depth.restriction) {
      if (is.numeric(s$depth.restriction) && length(s$depth.restriction) == 2) {
        q_cell[z >= s$depth.restriction[1] & z <= s$depth.restriction[2]] <- 0
      } else if (identical(s$depth.restriction, "random")) {
        q_cell[sample.int(n_cells, min(20, n_cells))] <- 0
      }
    }
    price <- matrix(s$price_fish[yr, ], n_cells, ns, byrow = TRUE)
    value <- rowSums(biomass * (1 - exp(-q_cell)) * price)
    if (sum(value^preference) <= 0) stop("No positive expected catch in year ", yr,
                                         call. = FALSE)
    effort <- stats::rmultinom(1, n_sets, value^preference / sum(value^preference))[, 1]

    # Sets: cell, vessel and effective catchability
    cell <- rep.int(seq_len(n_cells), effort)
    k <- length(cell)
    vessel <- sample.int(s$Nvessels, k, replace = TRUE)
    q_mean <- rep(q_vessel[vessel], ns)
    q_set <- matrix(if (isTRUE(s$do.tweedie)) rtweedie_cp(q_mean, s$xi, s$phi) else q_mean,
                    k, ns)
    catch <- allocate_catch(biomass, cell, q_set)

    records[[yr]] <- data.frame(year = yr, SetID = set_counter + seq_len(k),
      X = x[cell], Y = y[cell], depth = z[cell], vessel = vessel,
      catch$set, check.names = FALSE)
    names(records[[yr]])[-(1:6)] <- species
    set_counter <- set_counter + k

    # Population dynamics, then end-of-year movement
    biomass <- schaefer_update(biomass, carrying, s$r, catch$cell)
    if (!isTRUE(s$no.mvt)) {
      biomass <- vapply(seq_len(ns), function(i) move_biomass(movement[[i]], biomass[, i]),
                        numeric(n_cells))
    }
  }

  sets <- do.call(rbind, records)
  data <- do.call(rbind, lapply(seq_len(ns), function(i) {
    data.frame(sets[, c("year", "SetID", "X", "Y")],
               Xobs = coarse_coord(sets$X, s$obs_res),
               Yobs = coarse_coord(sets$Y, s$obs_res),
               Species = i, Sp = species[i], CPUE = sets[[species[i]]],
               depth = sets$depth, vessel = sets$vessel)
  }))
  rownames(data) <- NULL

  env <- data.frame(x = x, y = y, Depth = -z)
  out <- list(
    Data = data,
    Biomass = index_table(data, total_biomass, species),
    Biomass.XY = density_table(env, carrying, species),
    Env = env,
    settings = s
  )
  if (plots) out$Plots <- plot_om(out)
  out
}

#' Long table of initial densities by cell and species
#' @noRd
density_table <- function(env, density, species) {
  out <- do.call(rbind, lapply(seq_along(species), function(i) {
    data.frame(env, variable = species[i], Density = density[, i])
  }))
  out$variable <- factor(out$variable, levels = species)
  out
}
