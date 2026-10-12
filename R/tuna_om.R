#' Default settings for the longline tuna scenario
#'
#' Four tuna species in a conceptual Southern Hemisphere ocean (y points north,
#' towards the equator, so the north is warm).
#'
#' * **Horizontal habitat** is a normal preference for sea surface temperature
#'   (`sst_opt`, `sst_tol`). Yellowfin (YFT) and bigeye (BET) prefer warm water;
#'   southern bluefin (BFT) and albacore (ALB) prefer cooler water. Wide
#'   tolerances keep the species mixed rather than separated.
#' * **Vertical habitat** is a normal depth profile (`depth_mu`, `depth_sd`, in
#'   m). Within each pair one species is shallow (YFT, ALB) and one deep (BET,
#'   BFT). A thermocline-depth anomaly shifts all profiles up or down
#'   together in each cell, so the catch of a given gear configuration varies
#'   in space and, optionally, between years. The anomaly in year t and cell i
#'   is `fixed[i] + shift[t] + annual[i, t]`:
#'   - `fixed`: a spatial field that stays the same every year
#'     (`thermocline_sd`; 0 = same vertical profile everywhere);
#'   - `shift`: a basin-wide AR(1) shift (`thermocline_annual_sd`,
#'     `thermocline_annual_rho`), e.g. ENSO-like years;
#'   - `annual`: a new spatial field each year (`thermocline_annual_field_sd`).
#' * **Hook depth**: each HBF (hooks between floats) level has a typical
#'   (median) hook depth `HBF_depth`, roughly as in Fig. 2 of Bigelow et al.
#'   (2002). Each set's hook profile is shifted by a lognormal set-level error
#'   (`set_cv`; e.g. line angle, currents), so HBF is an imperfect proxy for
#'   depth. Hooks within a set are lognormally spread around the set's median
#'   (`hook_cv`).
#' * **Targeting**: species prices change in three 10-year blocks: ALB
#'   (shallow, cool water), then BET (deep, warm water), then BFT (deep, cool
#'   water), loosely following the historical shift from albacore to bigeye
#'   targeting and the later interest in recovering southern bluefin. Skippers
#'   choose where to fish and which HBF to use from the expected catch value.
#'
#' All values are conceptual and meant to be tuned. With the same seed,
#' scenarios that differ only in the thermocline settings share the same SST
#' field, the same fixed thermocline pattern (scaled by `thermocline_sd`) and
#' the same vessels; the fishery itself then evolves differently.
#'
#' Two presets are available:
#' * `"default"`: moderate targeting (prices differ by a factor of 2 between
#'   the target species and the others, `preference = 1.25`,
#'   `HBF_choice_temp = 0.2`); nominal CPUE is biased by roughly 10-20% by
#'   decade.
#' * `"stress_test"`: extreme targeting (prices differ by a factor of 4,
#'   `preference = 1.5`, `HBF_choice_temp = 0.15`); abrupt, fleet-wide target
#'   switches with nominal CPUE biased by up to a factor of two. These are the
#'   settings of the first runs of the vignette.
#'
#' @param preset `"default"` or `"stress_test"`.
#' @return A named list accepted by [simulate_tuna()].
#' @export
#' @examples
#' s <- tuna_settings()
#' extreme <- tuna_settings("stress_test")
#' extreme$price[c(1, 11, 21), ]
tuna_settings <- function(preset = c("default", "stress_test")) {
  preset <- match.arg(preset)
  s <- tuna_settings_default()
  if (preset == "stress_test") {
    s$price <- matrix(c(rep(c(1.0, 0.5, 0.5, 2.0), 10),   # years 1-10: ALB
                        rep(c(1.0, 2.0, 1.0, 0.5), 10),   # years 11-20: BET
                        rep(c(0.5, 1.0, 2.0, 0.5), 10)),  # years 21-30: BFT
                      ncol = 4, byrow = TRUE)
    s$preference <- 1.5
    s$HBF_choice_temp <- 0.15
  }
  s
}

#' Default tuna settings (see tuna_settings())
#' @noRd
tuna_settings_default <- function() {
  list(
    n_years = 30,
    nx = 50, ny = 50,          # grid size (cells)
    obs_res = 5,               # reported positions: centres of 5 x 5 blocks

    # Species ----------------------------------------------------------------
    species  = c("YFT", "BET", "BFT", "ALB"),
    B0       = c(400000, 300000, 150000, 350000), # carrying capacity (t)
    r        = c(0.30, 0.18, 0.12, 0.18),         # intrinsic growth rate
    B_start  = c(0.7, 0.6, 0.5, 0.7), # initial biomass / carrying capacity (1 = unfished)
    sst_opt  = c(25, 24, 20, 21),     # preferred SST (deg C)
    sst_tol  = c(3, 3.5, 3, 3.5),     # SD of SST preference (deg C)
    depth_mu = c(90, 260, 240, 110),  # mean depth (m): YFT, ALB shallow; BET, BFT deep
    depth_sd = c(60, 50, 55, 60),     # SD of depth (m)
    dispersal = rep(1.5, 4),          # mean movement distance per year (cells)

    # Environment ------------------------------------------------------------
    sst_mean = 23,            # mean SST (deg C)
    sst_gradient = 10,        # SST increase from south to north (deg C)
    sst_noise = 1,            # SD of the random SST field (deg C)
    sst_scale = 6,            # correlation scale of the SST field (cells)
    thermocline_sd = 30,      # SD of the fixed thermocline-depth anomaly field (m)
    thermocline_scale = 8,    # correlation scale of the anomaly fields (cells)
    thermocline_annual_sd = 0,       # SD of a basin-wide annual shift (m), AR(1);
                                     #   ~10-15 m is a realistic ENSO-like range
    thermocline_annual_rho = 0.5,    # autocorrelation of the annual shift
    thermocline_annual_field_sd = 0, # SD of an extra year-specific anomaly field (m);
                                     #   e.g. 5-10 m
    max_depth = 500,          # depth range of vertical profiles (m)
    depth_bin = 10,           # vertical resolution (m)

    # Longline gear ------------------------------------------------------------
    HBF_levels = c(5, 10, 15, 20),
    HBF_depth = c(110, 160, 205, 250), # median hook depth for each HBF (m)
    set_cv = 0.2,             # between-set error in depth: how well HBF predicts depth
    hook_cv = 0.35,           # within-set spread of hook depths

    # Fleet -------------------------------------------------------------------
    # columns: YFT, BET, BFT, ALB
    price = matrix(c(rep(c(1.0, 0.75, 0.75, 1.5), 10),   # years 1-10: ALB
                     rep(c(1.0, 1.5, 1.0, 0.75), 10),    # years 11-20: BET
                     rep(c(0.75, 1.0, 1.5, 0.75), 10)),  # years 21-30: BFT
                   ncol = 4, byrow = TRUE),
    n_vessels = 30,
    effort = 800,             # asymptotic number of sets per year (~19,000 in 30 years)
    effort_cv = 0.2,
    q = 2,                    # mean catchability per set at full vertical overlap;
                              #   with effort, sets how strongly stocks are fished down
    vessel_cv = 0.2,          # between-vessel CV in catchability
    tweedie_power = 1.2,      # set-level noise in catchability
    tweedie_phi = 0.1,
    preference = 1.25,        # concentration of effort on high-value cells
    HBF_choice_temp = 0.2,    # lower = skippers pick the best HBF more strictly
    avail_ref = 100           # vertical availability of 1 = fish spread evenly
                              # over a 100-m layer fished by all hooks
  )
}

#' Simulate a longline fishery for four tuna species with vertical targeting
#'
#' Extends the spatial multispecies simulator of Thorson et al. (2017) and
#' Ono et al. (2015), as implemented in [simulate_bathymetric()], with a
#' pelagic horizontal habitat (SST), species vertical distributions, longline
#' hook depth by HBF and vertical targeting, following the habitat-based
#' approach of Bigelow et al. (2002). A full description of both models is in
#' `system.file("appendix", "cpueTargeting_simulator_appendix.html",
#' package = "cpueTargeting")`. The
#' catchability of a set for a species is the vessel's catchability times the
#' overlap between the set's hook-depth distribution and the species' vertical
#' distribution in that cell (Bigelow's effective effort). See
#' [tuna_settings()] for the processes and parameters.
#'
#' @param settings A named list, usually a modified [tuna_settings()].
#' @param seed Optional integer seed for reproducibility.
#' @param plots Logical; build the diagnostic ggplots in `Plots`.
#'
#' @return A list with
#' * `Data`: one row per set and species with `year`, `SetID`, `X`, `Y`,
#'   `Xobs`, `Yobs`, `vessel`, `HBF`, `sst`, `thermocline` (anomaly in that
#'   year, m),
#'   `Species`, `Sp`, `CPUE`, and the latent truth `hook_depth` (mean hook
#'   hook depth of the set), `vertical_avail` (realised vertical overlap) and
#'   `hbs_avail` (overlap predicted by a habitat model from HBF, the expected
#'   hook profile and the average thermocline, as in a habitat-based
#'   standardisation);
#' * `Biomass`: annual nominal and true indices, scaled to mean 1;
#' * `Biomass.XY`: initial density by cell and species;
#' * `Env`: SST and the fixed thermocline anomaly by cell;
#' * `Thermocline`: thermocline anomaly by cell and year;
#' * `Plots`: diagnostic ggplots (if `plots = TRUE`);
#' * `settings`: the settings used.
#'
#' @references Bigelow, K.A., Hampton, J. and Miyabe, N. (2002) Application of
#'   a habitat-based model to estimate effective longline fishing effort and
#'   relative abundance of Pacific bigeye tuna (*Thunnus obesus*). Fisheries
#'   Oceanography 11, 143-155.
#'
#'   Thorson, J.T., Fonner, R., Haltuch, M.A., Ono, K. and Winker, H. (2017)
#'   Accounting for spatiotemporal variation and fisher targeting when
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
#' om <- simulate_tuna(seed = 42)
#' head(om$Data)
#' om$Plots$CPUE
#' }
simulate_tuna <- function(settings = tuna_settings(), seed = NULL, plots = TRUE) {
  if (!is.list(settings) || is.null(names(settings))) {
    stop("`settings` must be a named list.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  s <- utils::modifyList(tuna_settings(), settings)
  check_tuna_settings(s)
  species <- s$species
  ns <- length(species)

  # Grid and environment -------------------------------------------------------
  n_cells <- s$nx * s$ny
  x <- rep(seq_len(s$nx), s$ny) + 0.5
  y <- rep(seq_len(s$ny), each = s$nx) + 0.5
  sst <- s$sst_mean + s$sst_gradient * (y - mean(y)) / s$ny +
    as.vector(simulate_grf(s$nx, s$ny, s$sst_noise, s$sst_scale))
  thermocline <- as.vector(simulate_grf(s$nx, s$ny, s$thermocline_sd, s$thermocline_scale))

  # Horizontal habitat and movement --------------------------------------------
  habitat <- vapply(seq_len(ns), function(i) stats::dnorm(sst, s$sst_opt[i], s$sst_tol[i]),
                    numeric(n_cells))
  kernels <- dispersal_kernels(x, y, s$dispersal)
  movement <- lapply(seq_len(ns), function(i) movement_operator(kernels[[i]], habitat[, i]))

  # Vertical structure ------------------------------------------------------------
  edges <- seq(0, s$max_depth, by = s$depth_bin)
  avail_scale <- s$avail_ref / s$depth_bin
  hbf_profile <- nominal_hook_profiles(s, edges)   # expected hook profile by HBF
  # What an analyst's habitat model would predict for each HBF (used for HBS):
  # expected hook profiles and fish profiles under the fixed (average)
  # thermocline, i.e. without the annual anomalies.
  hbs_avail_cell <- lapply(fish_profiles(thermocline, s, edges),
                           function(p) p %*% t(hbf_profile) * avail_scale)
  annual_variation <- s$thermocline_annual_sd > 0 || s$thermocline_annual_field_sd > 0
  shift <- 0
  thermo_years <- matrix(0, n_cells, s$n_years)

  # Carrying capacity = B0 (with lognormal noise), spread over the stable
  # distribution of the movement process; stocks start at B_start x capacity --
  k_total <- s$B0 * stats::rlnorm(ns, 0, 0.2)
  carrying <- vapply(seq_len(ns), function(i) movement[[i]]$stable * k_total[i],
                     numeric(n_cells))
  biomass <- sweep(carrying, 2, s$B_start, "*")

  # Fleet ------------------------------------------------------------------------
  sig_vessel <- sqrt(log(s$vessel_cv^2 + 1))
  q_vessel <- stats::rlnorm(s$n_vessels, log(s$q) - sig_vessel^2 / 2, sig_vessel)
  n_hbf <- length(s$HBF_levels)

  records <- vector("list", s$n_years)
  total_biomass <- matrix(0, s$n_years, ns)
  set_counter <- 0L

  for (yr in seq_len(s$n_years)) {
    total_biomass[yr, ] <- colSums(biomass)

    # Thermocline anomaly this year and the resulting vertical distributions
    # (stationary AR(1): the first year is drawn from the marginal distribution)
    rho <- if (yr == 1) 0 else s$thermocline_annual_rho
    shift <- rho * shift + sqrt(1 - rho^2) * s$thermocline_annual_sd * stats::rnorm(1)
    annual_field <- as.vector(simulate_grf(s$nx, s$ny, s$thermocline_annual_field_sd,
                                           s$thermocline_scale))
    thermo_years[, yr] <- thermocline + shift + annual_field
    if (yr == 1 || annual_variation) {
      fish_profile <- fish_profiles(thermo_years[, yr], s, edges)
      # expected availability for each HBF as seen by skippers this year
      expected_avail <- lapply(fish_profile, function(p) p %*% t(hbf_profile) * avail_scale)
    }
    n_sets <- annual_effort(yr, s$effort, s$effort_cv)

    # Expected value of each cell x HBF combination
    value_hbf <- matrix(0, n_cells, n_hbf)
    for (i in seq_len(ns)) {
      value_hbf <- value_hbf + s$price[yr, i] * biomass[, i] * (1 - exp(-s$q * expected_avail[[i]]))
    }
    best <- value_hbf[cbind(seq_len(n_cells), max.col(value_hbf, ties.method = "first"))]
    if (sum(best) <= 0) stop("No positive expected catch in year ", yr, call. = FALSE)

    # Where to fish (spatial targeting)
    p_cell <- best^s$preference / sum(best^s$preference)
    cell <- rep.int(seq_len(n_cells), stats::rmultinom(1, n_sets, p_cell)[, 1])
    k <- length(cell)

    # Which HBF to use (vertical targeting)
    w <- (value_hbf[cell, , drop = FALSE] / best[cell])^(1 / s$HBF_choice_temp)
    hbf <- s$HBF_levels[sample_rows(w)]

    # Realised hook depths and vertical overlap of each set
    hooks <- set_hook_profiles(hbf, s, edges)
    avail <- vapply(seq_len(ns), function(i) {
      rowSums(hooks$profile * fish_profile[[i]][cell, , drop = FALSE]) * avail_scale
    }, numeric(k))
    avail <- matrix(avail, k, ns)
    hbs <- vapply(seq_len(ns), function(i) {
      hbs_avail_cell[[i]][cbind(cell, match(hbf, s$HBF_levels))]
    }, numeric(k))
    hbs <- matrix(hbs, k, ns)

    # Catch
    vessel <- sample.int(s$n_vessels, k, replace = TRUE)
    q_set <- matrix(rtweedie_cp(rep(q_vessel[vessel], ns), s$tweedie_power, s$tweedie_phi),
                    k, ns) * avail
    catch <- allocate_catch(biomass, cell, q_set)

    records[[yr]] <- list(
      sets = data.frame(year = yr, SetID = set_counter + seq_len(k),
                        X = x[cell], Y = y[cell], vessel = vessel, HBF = hbf,
                        sst = sst[cell], thermocline = thermo_years[cell, yr],
                        hook_depth = hooks$median_depth),
      catch = catch$set, avail = avail, hbs = hbs)
    set_counter <- set_counter + k

    biomass <- schaefer_update(biomass, carrying, s$r, catch$cell)
    biomass <- vapply(seq_len(ns), function(i) move_biomass(movement[[i]], biomass[, i]),
                      numeric(n_cells))
  }

  sets <- do.call(rbind, lapply(records, `[[`, "sets"))
  catch <- do.call(rbind, lapply(records, `[[`, "catch"))
  avail <- do.call(rbind, lapply(records, `[[`, "avail"))
  hbs <- do.call(rbind, lapply(records, `[[`, "hbs"))
  sets$Xobs <- coarse_coord(sets$X, s$obs_res)
  sets$Yobs <- coarse_coord(sets$Y, s$obs_res)
  set_cols <- c("year", "SetID", "X", "Y", "Xobs", "Yobs", "vessel", "HBF",
                "sst", "thermocline")
  data <- do.call(rbind, lapply(seq_len(ns), function(i) {
    data.frame(sets[, set_cols], Species = i, Sp = species[i], CPUE = catch[, i],
               hook_depth = sets$hook_depth, vertical_avail = avail[, i],
               hbs_avail = hbs[, i])
  }))
  rownames(data) <- NULL

  env <- data.frame(x = x, y = y, SST = sst, thermocline = thermocline)
  out <- list(
    Data = data,
    Biomass = index_table(data, total_biomass, species),
    Biomass.XY = density_table(env, carrying, species),
    Env = env,
    Thermocline = data.frame(x = rep(x, s$n_years), y = rep(y, s$n_years),
                             year = rep(seq_len(s$n_years), each = n_cells),
                             anomaly = as.vector(thermo_years)),
    settings = s
  )
  if (plots) out$Plots <- plot_om(out)
  out
}

#' Vertical distribution of each species in each cell
#'
#' Normal depth profile with mean `depth_mu + anomaly` and SD `depth_sd`,
#' as proportions per depth bin.
#' @return A list with one cells x depth-bins matrix per species.
#' @noRd
fish_profiles <- function(anomaly, s, edges) {
  mids <- edges[-1] - diff(edges) / 2
  lapply(seq_along(s$species), function(i) {
    p <- stats::dnorm(matrix(mids, length(anomaly), length(mids), byrow = TRUE),
                      mean = s$depth_mu[i] + anomaly, sd = s$depth_sd[i])
    p / rowSums(p)
  })
}

#' Draw one column index per row with probabilities proportional to `w`
#' @noRd
sample_rows <- function(w) {
  w[!is.finite(w)] <- 0
  w[rowSums(w) == 0, ] <- 1
  cum <- w
  if (ncol(w) > 1) for (j in 2:ncol(w)) cum[, j] <- cum[, j - 1] + w[, j]
  u <- stats::runif(nrow(w)) * cum[, ncol(w)]
  pmin(rowSums(u > cum) + 1L, ncol(w))
}

#' Convert a coefficient of variation to a lognormal SD (log scale)
#' @noRd
cv_to_sdlog <- function(cv) sqrt(log(1 + cv^2))

#' Hook-depth distribution (proportion of hooks per depth bin) of individual sets
#'
#' A set with a given HBF has median hook depth `HBF_depth * exp(e)` with
#' `e ~ N(0, sdlog(set_cv))` (set-level error), and its hooks are lognormally
#' spread around that median with `sdlog(hook_cv)`.
#'
#' @return A list with `profile` (sets x depth bins) and `median_depth`.
#' @noRd
set_hook_profiles <- function(hbf, s, edges) {
  median_depth <- s$HBF_depth[match(hbf, s$HBF_levels)] *
    exp(cv_to_sdlog(s$set_cv) * stats::rnorm(length(hbf)))
  list(profile = lognormal_bin_mass(median_depth, cv_to_sdlog(s$hook_cv), edges),
       median_depth = median_depth)
}

#' Expected hook-depth distribution of each HBF level
#'
#' Lognormal with median `HBF_depth` and the combined between- and within-set
#' spread: the aggregate distribution an analyst would use (cf. Fig. 2 of
#' Bigelow et al. 2002). Rows are HBF levels, columns depth bins.
#' @noRd
nominal_hook_profiles <- function(s, edges) {
  sdlog <- sqrt(cv_to_sdlog(s$set_cv)^2 + cv_to_sdlog(s$hook_cv)^2)
  lognormal_bin_mass(s$HBF_depth, sdlog, edges)
}

#' Proportion of a lognormal depth distribution in each depth bin
#'
#' One row per median; hooks deeper than the last edge are reassigned
#' proportionally, so every row sums to one.
#' @noRd
lognormal_bin_mass <- function(median, sdlog, edges) {
  cdf <- stats::plnorm(matrix(edges, length(median), length(edges), byrow = TRUE),
                       meanlog = log(median), sdlog = sdlog)
  mass <- cdf[, -1, drop = FALSE] - cdf[, -ncol(cdf), drop = FALSE]
  mass / rowSums(mass)
}

#' Basic checks on tuna settings
#' @noRd
check_tuna_settings <- function(s) {
  ns <- length(s$species)
  per_species <- c("B0", "r", "B_start", "sst_opt", "sst_tol", "depth_mu", "depth_sd", "dispersal")
  bad <- per_species[vapply(per_species, function(v) length(s[[v]]) != ns, logical(1))]
  if (length(bad)) stop("These settings need one value per species: ",
                        paste(bad, collapse = ", "), call. = FALSE)
  if (!is.matrix(s$price) || ncol(s$price) != ns || nrow(s$price) < s$n_years) {
    stop("`price` must be a matrix with one row per year and one column per species.",
         call. = FALSE)
  }
  if (any(s$sst_tol <= 0) || any(s$depth_sd <= 0)) {
    stop("`sst_tol` and `depth_sd` must be positive.", call. = FALSE)
  }
  if (any(s$B_start <= 0)) stop("`B_start` must be positive.", call. = FALSE)
  if (s$HBF_choice_temp <= 0) stop("`HBF_choice_temp` must be positive.", call. = FALSE)
  if (abs(s$thermocline_annual_rho) >= 1) {
    stop("`thermocline_annual_rho` must be between -1 and 1.", call. = FALSE)
  }
  if (any(c(s$thermocline_sd, s$thermocline_annual_sd,
            s$thermocline_annual_field_sd, s$set_cv) < 0) || s$hook_cv <= 0) {
    stop("Thermocline SDs and `set_cv` cannot be negative; `hook_cv` must be positive.",
         call. = FALSE)
  }
  if (length(s$HBF_depth) != length(s$HBF_levels) || any(s$HBF_depth <= 0)) {
    stop("`HBF_depth` needs one positive depth per HBF level.", call. = FALSE)
  }
  if (s$max_depth %% s$depth_bin != 0) {
    stop("`max_depth` must be a multiple of `depth_bin`.", call. = FALSE)
  }
  invisible(TRUE)
}
