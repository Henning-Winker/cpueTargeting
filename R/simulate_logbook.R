#' Default settings for the logbook simulator
#'
#' Ten species in a mixed demersal-pelagic fishery. Demersal species live
#' along a continuous depth gradient and on a patchy substrate (reef patches
#' within sediment), which gives four habitats as in Winker et al. (2014):
#' shallow reef, deep reef, shallow sediment and deep sediment. Pelagic
#' species follow sea surface temperature and distance offshore and migrate
#' seasonally through the area. Species names are illustrative.
#'
#' Main elements (all can be changed):
#' * `species`: one row per species with `group` (demersal/pelagic), depth
#'   preference (`depth_mu`, `depth_cv`, lognormal), substrate affinity
#'   (`reef`, `sediment`, multipliers of density), SST preference (`sst_mu`,
#'   `sst_sd`, normal; `NA` for none), seasonal availability (`season_amp`,
#'   0 = resident, 1 = absent at the low point; `season_peak`, month),
#'   `dispersal`, Schaefer `B0` and `r`, `B_start` (initial biomass / B0) and
#'   base `price`.
#' * `tactics`: catchability of each species by tactic (rows `Reef`,
#'   `Sediment`, `Pelagic`), i.e. gear selectivity.
#' * `price_periods`: years at which prices change and multipliers for
#'   demersal and pelagic species in each period (default: demersal premium,
#'   then a pelagic boom, then neutral prices).
#' * Fleet: `vessel_types` (numbers of vessels by type), `home` (home
#'   tactics of each type), `loyalty` (share of a type's trips restricted to
#'   its home tactics; the other trips choose among all tactics by expected
#'   value), vessel catchability CV, late entries with efficiency creep
#'   (`entry_share`, `creep`), early exits (`exit_share`). The default has
#'   demersal and pelagic specialists with loyalty 0.8 and mixed vessels that
#'   choose purely by expected value. (`allowed`, from earlier versions, is
#'   still accepted and means `home` with loyalty 1.)
#' * `trips`: trips per year (default: an increase, cut by 20% with the
#'   closure); `preference`: how strongly effort concentrates on the most
#'   valuable cells and tactics.
#' * `closure`: year, depth range and tactics of a spatial closure (default:
#'   sediment gear banned deeper than 300 m from year 17).
#' * `logbook`: how trips are reported (see [simulate_logbook()]).
#'
#' @return A named list for [simulate_logbook()].
#' @export
logbook_settings <- function() {
  species <- data.frame(
    species = c("Dentex dentex", "Pagellus erythrinus", "Pagellus bogaraveo",
                "Merluccius merluccius", "Conger conger", "Phycis blennoides",
                "Thunnus alalunga", "Coryphaena hippurus", "Seriola dumerili",
                "Scomber colias"),
    group = rep(c("demersal", "pelagic"), c(6, 4)),
    habitat = c("shallow reef", "shallow sediment", "deep reef", "deep sediment",
                "reef, all depths", "deep, mixed", "offshore", "warm, offshore",
                "coastal, reef-associated", "coastal"),
    depth_mu = c(60, 50, 350, 400, 150, 500, 700, 300, 100, 60),
    depth_cv = c(0.5, 0.5, 0.3, 0.35, 1.0, 0.3, 0.6, 1.0, 0.7, 1.0),
    reef = c(1, 0.05, 1, 0.05, 1, 0.5, 1, 1, 1, 1),
    sediment = c(0.02, 1, 0.02, 1, 0.3, 1, 1, 1, 0.3, 1),
    sst_mu = c(NA, NA, NA, NA, NA, NA, 18, 21, 20, 18),
    sst_sd = c(NA, NA, NA, NA, NA, NA, 1.5, 1.2, 2, 2.5),
    season_amp = c(0, 0, 0, 0, 0, 0, 0.7, 0.9, 0.4, 0.3),
    season_peak = c(1, 1, 1, 1, 1, 1, 9, 9, 8, 5),
    dispersal = c(0.5, 0.5, 0.5, 0.8, 0.5, 0.8, 2, 2, 1, 1),
    B0 = c(500, 900, 800, 2500, 600, 700, 3000, 1500, 600, 4000),
    r = c(0.25, 0.35, 0.18, 0.30, 0.30, 0.30, 0.40, 0.80, 0.35, 0.80),
    B_start = c(0.9, 0.8, 0.7, 0.6, 0.9, 0.9, 0.9, 1, 0.8, 1),
    price = c(12, 8, 15, 5, 3, 2, 6, 7, 9, 2),
    stringsAsFactors = FALSE)
  tactics <- rbind(
    Reef     = c(1.00, 0.20, 1.00, 0.05, 0.80, 0.10, 0.001, 0.001, 0.40, 0.20),
    Sediment = c(0.10, 1.00, 0.10, 1.00, 0.50, 1.00, 0.001, 0.001, 0.02, 0.02),
    Pelagic  = c(0.005, 0.001, 0.001, 0.001, 0.001, 0.001, 1.00, 1.00, 0.80, 0.60))
  colnames(tactics) <- species$species
  list(
    n_years = 25,
    start_year = 2000,
    grid = c(40, 40),
    shelf_depth = c(10, 1000),   # depth at the inshore and offshore edge (m)
    shelf_power = 1.5,           # > 1: narrow shelf, steeper slope offshore
    depth_sd = 40, depth_scale = 4,
    reef_share = 0.25,           # share of cells with reef substrate
    reef_scale = 2.5,            # patch size of reefs (cells)
    sst_range = c(17, 22),       # SST from north to south
    sst_sd = 0.5, sst_scale = 8,
    species = species,
    tactics = tactics,
    q = 0.06,                    # catchability scale per set
    price_periods = list(years = c(1, 9, 17), demersal = c(1.4, 0.7, 1),
                         pelagic = c(0.6, 1.6, 1)),
    price = NULL,                # optional years x species price matrix
    vessel_types = c(demersal = 12, pelagic = 8, mixed = 20),
    home = list(demersal = c("Reef", "Sediment"), pelagic = "Pelagic",
                mixed = c("Reef", "Sediment", "Pelagic")),
    loyalty = c(demersal = 0.8, pelagic = 0.8, mixed = 0),  # share of trips on home tactics
    vessel_cv = 0.25,
    entry_share = 0.3,           # share of vessels entering after year 1
    creep = 0.02,                # efficiency gain per year of later entry
    exit_share = 0.2,
    trips = NULL,                # trips per year (one value per year), default curve
    trips_cv = 0.1,
    effort_season_amp = 0.3,     # more trips in summer
    hours_mean = 8, hours_cv = 0.3,
    preference = 1.5,
    closure = list(year = 17, depth = c(300, Inf), tactics = "Sediment"),
    tweedie_power = 1.3, tweedie_phi = 1,
    logbook = list(
      obs_res = 5,               # positions reported as centres of obs_res x obs_res blocks
      origin = c(lon = 3, lat = 36), cell_size = 0.1,
      missing_position = 0.02,
      two_sets = 0.3,            # share of trips with a second set (reported as one trip)
      switch_tactic = 0.7,       # probability that the second set uses another tactic
                                 # (with the default fleet ~18% of trips mix two tactics)
      kg_trip = 60,              # median total catch per trip (kg)
      round_to = 0.5,
      underreport = list(max_price = 3, below_kg = 2, prob = 0.5),
      synonyms = c("Pagellus bogaraveo" = "Blackspot seabream"),
      synonym_rate = 0.1)
  )
}

#' Simulate a realistic multispecies logbook
#'
#' A spatial operating model of a mixed demersal-pelagic fishery with
#' discrete tactic choice, combining the spatial multispecies simulator
#' ([simulate_bathymetric()]) with the habitat-based tactics of Winker et al.
#' (2014) ([simulate_cluster()]), and a reporting layer that turns sets into
#' a typical trip-level logbook.
#'
#' **Environment.** A shelf-to-slope bathymetry (depth increases offshore,
#' plus a Gaussian random field), a patchy substrate (reef patches from a
#' thresholded random field) and an SST gradient. Demersal densities are the
#' product of a depth preference and a substrate affinity, so the four
#' habitats shallow/deep x reef/sediment arise as continuous density
#' gradients. Pelagic densities follow SST and depth, and their availability
#' in the area changes seasonally (migrants).
#'
#' **Dynamics.** Schaefer production in each cell once a year, then movement
#' (exponential distance kernel weighted by habitat). Catches are removed
#' every month.
#'
#' **Fleet.** Each month, trips choose a tactic (`Reef`, `Sediment`,
#' `Pelagic`) and a cell with probability proportional to
#' `expected value ^ preference`, where expected value is the sum over species
#' of price x catchability of the tactic x available biomass. A share
#' `loyalty` of a specialist's trips is restricted to its home tactics; its
#' other trips, and all trips of mixed vessels, choose among all tactics. Prices change between periods, which
#' shifts tactics; a closure bans a tactic from part of the area; vessels
#' enter (more efficient) and leave.
#'
#' **Logbook.** Some trips make a second set (possibly with another tactic)
#' and report the sum, with one position; positions are reported at a coarse
#' resolution and are sometimes missing; catches are rounded; small catches
#' of low-value species are often not recorded; one species is sometimes
#' recorded under another name. The tactic is not reported.
#'
#' @param settings Settings list, usually a modified [logbook_settings()].
#' @param seed Optional seed.
#' @param plots Logical; build diagnostic plots.
#' @return A list with
#' * `logbook`: one row per trip and recorded species (`trip`, `year`,
#'   `month`, `vessel`, `hours`, `lon`, `lat`, `species`, `kg`), ready for
#'   [cpue_data()];
#' * `trips`: the truth per trip (vessel type, whether the trip stayed on home tactics (`loyal`), tactic of each set, number of
#'   sets, cell, depth, substrate, habitat);
#' * `truth`: true total biomass by year and species (also scaled to mean 1);
#' * `nominal`: nominal CPUE (kg/hour) by year and species from the logbook;
#' * `env`: depth, substrate and SST by cell; `density`: initial density by
#'   cell and species; `prices`: prices by year and species;
#' * `Plots` (if `plots = TRUE`) and `settings`.
#' @seealso [logbook_truth()], [plot_logbook_om()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_logbook(seed = 1)
#' head(sim$logbook)
#' sim$Plots$density
#' sim$Plots$tactics
#' x <- cpue_data(sim$logbook, year = "year", species = "species", catch = "kg",
#'                effort = "hours", id = "trip", vessel = "vessel",
#'                lon = "lon", lat = "lat", covariates = "month",
#'                synonyms = c("Blackspot seabream" = "Pagellus bogaraveo"))
#' }
simulate_logbook <- function(settings = logbook_settings(), seed = NULL, plots = TRUE) {
  s <- merge_logbook_settings(settings)
  if (!is.null(seed)) set.seed(seed)
  sp <- s$species
  ns <- nrow(sp)
  spn <- sp$species
  ny_yr <- s$n_years
  tac_names <- rownames(s$tactics)

  # Environment
  nx <- s$grid[1]
  ny <- s$grid[2]
  n_cells <- nx * ny
  x <- rep(seq_len(nx), ny)
  y <- rep(seq_len(ny), each = nx)
  depth <- s$shelf_depth[1] + diff(s$shelf_depth) * ((x - 1) / (nx - 1))^s$shelf_power +
    as.vector(simulate_grf(nx, ny, s$depth_sd, s$depth_scale))
  depth <- pmax(depth, 5)
  reef_field <- as.vector(simulate_grf(nx, ny, 1, s$reef_scale))
  reef <- reef_field >= stats::quantile(reef_field, 1 - s$reef_share)
  sst <- s$sst_range[2] - diff(s$sst_range) * (y - 1) / (ny - 1) +
    as.vector(simulate_grf(nx, ny, s$sst_sd, s$sst_scale))

  habitat <- vapply(seq_len(ns), function(i) {
    h <- habitat_density(depth, "Lognorm", sp$depth_mu[i], sp$depth_cv[i])
    h <- h * ifelse(reef, sp$reef[i], sp$sediment[i])
    if (!is.na(sp$sst_mu[i])) h <- h * stats::dnorm(sst, sp$sst_mu[i], sp$sst_sd[i])
    h / max(h)
  }, numeric(n_cells))
  kernels <- dispersal_kernels(x, y, sp$dispersal)
  movement <- lapply(seq_len(ns), function(i) movement_operator(kernels[[i]], habitat[, i]))
  rm(kernels)
  carrying <- vapply(seq_len(ns), function(i) movement[[i]]$stable * sp$B0[i], numeric(n_cells))
  biomass <- sweep(carrying, 2, sp$B_start, `*`)

  # Prices
  price <- if (!is.null(s$price)) as.matrix(s$price) else {
    period <- findInterval(seq_len(ny_yr), s$price_periods$years)
    t(vapply(seq_len(ny_yr), function(yr) {
      k <- period[yr]
      sp$price * ifelse(sp$group == "demersal", s$price_periods$demersal[k],
                        s$price_periods$pelagic[k])
    }, numeric(ns)))
  }

  # Fleet
  types <- rep(names(s$vessel_types), s$vessel_types)
  nv <- length(types)
  entry <- rep(1L, nv)
  late <- sample.int(nv, round(s$entry_share * nv))
  late_years <- 2:max(2, ceiling(0.7 * ny_yr))
  entry[late] <- late_years[sample.int(length(late_years), length(late), replace = TRUE)]
  exit <- rep(ny_yr, nv)
  candidates <- setdiff(seq_len(nv), late)
  early <- candidates[sample.int(length(candidates), min(length(candidates),
                                                          round(s$exit_share * nv)))]
  exit_years <- max(2, ceiling(0.3 * ny_yr)):max(2, ny_yr - 1)
  exit[early] <- exit_years[sample.int(length(exit_years), length(early), replace = TRUE)]
  sig_v <- cv_to_sdlog(s$vessel_cv)
  q_vessel <- stats::rlnorm(nv, -sig_v^2 / 2, sig_v) * exp(s$creep * (entry - 1))

  trips_year <- if (!is.null(s$trips)) rep_len(s$trips, ny_yr) else {
    base <- 400 + 800 * stats::plogis((seq_len(ny_yr) - 8) / 2.5)
    cut <- if (!is.null(s$closure)) ifelse(seq_len(ny_yr) >= s$closure$year, 0.8, 1) else 1
    round(base * cut)
  }
  month_w <- 1 + s$effort_season_amp * cos(2 * pi * (1:12 - 7) / 12)
  avail <- t(vapply(1:12, function(m) {
    (1 - sp$season_amp) + sp$season_amp * (1 + cos(2 * pi * (m - sp$season_peak) / 12)) / 2
  }, numeric(ns)))
  q_tac <- s$tactics[, spn, drop = FALSE] * s$q
  sig_h <- cv_to_sdlog(s$hours_cv)
  sig_t <- cv_to_sdlog(s$trips_cv)

  all_sets <- list()
  catch_sets <- list()
  true_b <- matrix(0, ny_yr, ns)
  trip_counter <- 0L
  for (yr in seq_len(ny_yr)) {
    true_b[yr, ] <- colSums(biomass)
    n_trips <- round(trips_year[yr] * stats::rlnorm(1, -sig_t^2 / 2, sig_t))
    trips_m <- stats::rmultinom(1, n_trips, month_w)[, 1]
    closed <- rep(FALSE, n_cells)
    if (!is.null(s$closure) && yr >= s$closure$year) {
      closed <- depth >= s$closure$depth[1] & depth <= s$closure$depth[2]
    }
    active <- which(entry <= yr & exit >= yr)
    act_types <- factor(types[active], levels = names(s$vessel_types))
    for (m in 1:12) {
      if (trips_m[m] == 0 || !length(active)) next
      b_av <- sweep(biomass, 2, avail[m, ], `*`)
      value <- (b_av * matrix(price[yr, ], n_cells, ns, byrow = TRUE)) %*% t(q_tac)
      if (any(closed)) value[closed, intersect(s$closure$tactics, tac_names)] <- 0
      n_type <- stats::rmultinom(1, trips_m[m], as.numeric(table(act_types)))[, 1]
      names(n_type) <- levels(act_types)
      month_sets <- list()
      for (ty in names(n_type)[n_type > 0]) {
        # loyal trips use the home tactics only, the others all tactics
        n_loyal <- stats::rbinom(1, n_type[[ty]], s$loyalty[[ty]])
        groups <- list(home = list(n = n_loyal, tac = intersect(s$home[[ty]], tac_names)),
                       free = list(n = n_type[[ty]] - n_loyal, tac = tac_names))
        for (g in names(groups)) {
          allowed <- groups[[g]]$tac
          if (groups[[g]]$n == 0 || !length(allowed)) next
          w <- value[, allowed, drop = FALSE]^s$preference
          if (sum(w) <= 0) next
          alloc <- stats::rmultinom(1, groups[[g]]$n, as.vector(w) / sum(w))[, 1]
          idx <- rep.int(seq_along(alloc), alloc)
          k <- length(idx)
          cell <- (idx - 1) %% n_cells + 1
          tac <- allowed[(idx - 1) %/% n_cells + 1]
          pool <- active[act_types == ty]
          vessel <- pool[sample.int(length(pool), k, replace = TRUE)]
          trip <- trip_counter + seq_len(k)
          trip_counter <- trip_counter + k
          # second sets, reported with the first as one trip
          two <- stats::runif(k) < s$logbook$two_sets
          tac2 <- tac
          if (length(allowed) > 1) {
            sw <- two & stats::runif(k) < s$logbook$switch_tactic
            tac2[sw] <- vapply(tac[sw], function(t) {
              o <- setdiff(allowed, t)
              o[sample.int(length(o), 1)]
            }, character(1))
          }
          cell2 <- cell
          for (t in allowed) {
            j <- which(two & tac2 == t)
            wt <- w[, t]
            if (length(j) && sum(wt) > 0) cell2[j] <- sample.int(n_cells, length(j), TRUE, prob = wt)
          }
          month_sets[[paste(ty, g)]] <- data.frame(
            trip = c(trip, trip[two]), set = rep(1:2, c(k, sum(two))),
            year = yr, month = m, vessel = c(vessel, vessel[two]),
            vessel_type = ty, loyal = g == "home", tactic = c(tac, tac2[two]),
            cell = c(cell, cell2[two]), stringsAsFactors = FALSE)
        }
      }
      if (!length(month_sets)) next
      ms <- do.call(rbind, month_sets)
      ms$hours <- stats::rlnorm(nrow(ms), log(s$hours_mean) - sig_h^2 / 2, sig_h)
      q_eff <- q_tac[ms$tactic, , drop = FALSE] * q_vessel[ms$vessel] * (ms$hours / s$hours_mean)
      noise <- matrix(rtweedie_cp(rep(1, length(q_eff)), s$tweedie_power, s$tweedie_phi),
                      nrow(ms))
      catch <- allocate_catch(b_av, ms$cell, q_eff * noise)
      biomass <- pmax(biomass - catch$cell, 0)
      all_sets[[length(all_sets) + 1]] <- ms
      catch_sets[[length(catch_sets) + 1]] <- catch$set
    }
    biomass <- schaefer_update(biomass, carrying, sp$r, 0 * biomass)
    biomass <- vapply(seq_len(ns), function(i) move_biomass(movement[[i]], biomass[, i]),
                      numeric(n_cells))
  }
  sets <- do.call(rbind, all_sets)
  set_catch <- do.call(rbind, catch_sets)
  colnames(set_catch) <- spn

  # Trips: sum the sets, keep the first set's position and tactic
  trip_catch <- rowsum(set_catch, sets$trip)
  trip_hours <- as.numeric(rowsum(sets$hours, sets$trip))
  first <- sets[sets$set == 1, ]
  first <- first[match(as.integer(rownames(trip_catch)), first$trip), ]
  second <- sets[sets$set == 2, ]
  shallow <- depth < 150
  trips <- data.frame(
    trip = first$trip, year = first$year + s$start_year - 1, month = first$month,
    vessel = first$vessel, vessel_type = first$vessel_type, loyal = first$loyal,
    tactic = first$tactic,
    tactic2 = second$tactic[match(first$trip, second$trip)],
    n_sets = 1 + first$trip %in% second$trip, hours = trip_hours,
    X = x[first$cell], Y = y[first$cell], depth = depth[first$cell],
    substrate = ifelse(reef[first$cell], "reef", "sediment"),
    habitat = paste(ifelse(shallow[first$cell], "shallow", "deep"),
                    ifelse(reef[first$cell], "reef", "sediment")),
    stringsAsFactors = FALSE)
  rownames(trips) <- NULL

  # Reporting
  lb <- s$logbook
  tot <- rowSums(trip_catch)
  kg_scale <- if (any(tot > 0)) lb$kg_trip / stats::median(tot[tot > 0]) else 1
  kg <- trip_catch * kg_scale
  if (!is.null(lb$underreport)) {
    low <- which(sp$price <= lb$underreport$max_price)
    for (i in low) {
      drop <- kg[, i] > 0 & kg[, i] < lb$underreport$below_kg &
        stats::runif(nrow(kg)) < lb$underreport$prob
      kg[drop, i] <- 0
    }
  }
  if (!is.null(lb$round_to) && lb$round_to > 0) kg <- round(kg / lb$round_to) * lb$round_to
  rownames(kg) <- NULL
  lon <- lb$origin[["lon"]] + coarse_coord(trips$X + 0.5, lb$obs_res) * lb$cell_size
  lat <- lb$origin[["lat"]] + coarse_coord(trips$Y + 0.5, lb$obs_res) * lb$cell_size
  miss <- stats::runif(nrow(trips)) < lb$missing_position
  lon[miss] <- NA
  lat[miss] <- NA
  logbook <- do.call(rbind, lapply(seq_len(ns), function(i) {
    rows <- which(kg[, i] > 0)
    data.frame(trip = trips$trip[rows], year = trips$year[rows], month = trips$month[rows],
               vessel = sprintf("V%02d", trips$vessel[rows]), hours = round(trips$hours[rows], 1),
               lon = lon[rows], lat = lat[rows], species = spn[i], kg = kg[rows, i],
               stringsAsFactors = FALSE)
  }))
  if (length(lb$synonyms)) {
    for (orig in names(lb$synonyms)) {
      hit <- logbook$species == orig & stats::runif(nrow(logbook)) < lb$synonym_rate
      logbook$species[hit] <- lb$synonyms[[orig]]
    }
  }
  logbook <- logbook[order(logbook$trip, logbook$species), ]
  rownames(logbook) <- NULL

  years_cal <- seq_len(ny_yr) + s$start_year - 1
  truth <- do.call(rbind, lapply(seq_len(ns), function(i) {
    data.frame(year = years_cal, species = spn[i], group = sp$group[i],
               biomass = true_b[, i], relative = true_b[, i] / mean(true_b[, i]))
  }))
  truth$species <- factor(truth$species, levels = spn)
  cp <- kg / round(trips$hours, 1)
  nominal <- do.call(rbind, lapply(seq_len(ns), function(i) {
    m <- tapply(cp[, i], trips$year, mean)
    data.frame(year = as.numeric(names(m)), species = spn[i], nominal = as.numeric(m),
               relative = as.numeric(m) / mean(m))
  }))
  nominal$species <- factor(nominal$species, levels = spn)
  env <- data.frame(x = x, y = y, depth = depth,
                    substrate = ifelse(reef, "reef", "sediment"), sst = sst)
  density <- do.call(rbind, lapply(seq_len(ns), function(i) {
    data.frame(x = x, y = y, species = spn[i], density = carrying[, i])
  }))
  density$species <- factor(density$species, levels = spn)
  colnames(price) <- spn
  out <- list(logbook = logbook, trips = trips, truth = truth, nominal = nominal,
              env = env, density = density,
              prices = data.frame(year = years_cal, price, check.names = FALSE),
              settings = s)
  if (plots) out$Plots <- plot_logbook_om(out)
  out
}

#' Merge user settings with the defaults (species and tactics replaced whole)
#' @noRd
merge_logbook_settings <- function(settings) {
  if (!is.list(settings)) stop("`settings` must be a list.", call. = FALSE)
  def <- logbook_settings()
  whole <- c("species", "tactics", "price", "trips", "closure", "vessel_types", "home",
             "loyalty")
  # earlier versions: `allowed` = home tactics without leaving them
  if (!is.null(settings$allowed) && is.null(settings$home)) {
    settings$home <- settings$allowed
    if (is.null(settings$loyalty)) {
      settings$loyalty <- stats::setNames(rep(1, length(settings$allowed)),
                                          names(settings$allowed))
    }
  }
  settings$allowed <- NULL
  keep <- settings[intersect(names(settings), whole)]
  rest <- settings[setdiff(names(settings), whole)]
  s <- utils::modifyList(def, rest)
  for (n in names(keep)) s[n] <- list(keep[[n]])
  sp <- s$species
  need <- c("species", "group", "depth_mu", "depth_cv", "reef", "sediment", "sst_mu",
            "sst_sd", "season_amp", "season_peak", "dispersal", "B0", "r", "B_start", "price")
  if (!all(need %in% names(sp))) {
    stop("`species` needs columns: ", paste(setdiff(need, names(sp)), collapse = ", "),
         call. = FALSE)
  }
  if (!all(sp$species %in% colnames(s$tactics))) {
    stop("`tactics` needs one column per species.", call. = FALSE)
  }
  if (!all(names(s$vessel_types) %in% names(s$home))) {
    stop("`home` needs an entry for every vessel type.", call. = FALSE)
  }
  if (!all(names(s$vessel_types) %in% names(s$loyalty)) ||
      any(s$loyalty < 0 | s$loyalty > 1)) {
    stop("`loyalty` needs a value between 0 and 1 for every vessel type.", call. = FALSE)
  }
  if (s$n_years < 4) stop("`n_years` must be at least 4.", call. = FALSE)
  s
}

#' True relative biomass of one species from simulate_logbook()
#'
#' @param sim Output of [simulate_logbook()].
#' @param species Species name.
#' @return A data frame with `year` and `true`, e.g. for `plot_index(true = )`.
#' @export
logbook_truth <- function(sim, species) {
  tr <- sim$truth[sim$truth$species == species, ]
  if (!nrow(tr)) stop("Species not found: ", species, call. = FALSE)
  data.frame(year = tr$year, true = tr$relative)
}

#' Diagnostic plots of the logbook simulator
#'
#' @param sim Output of [simulate_logbook()].
#' @return A list of ggplots: `env` (depth with reef patches), `sst`,
#'   `density` (initial density by species), `tactics` (true share of trips by
#'   tactic over years), `effort` (trips by vessel type), `prices`, and `cpue`
#'   (nominal CPUE from the logbook against the true biomass).
#' @export
plot_logbook_om <- function(sim) {
  env <- sim$env
  reef_cells <- env[env$substrate == "reef", ]
  tr <- sim$trips
  lev <- rownames(sim$settings$tactics)
  pr <- sim$prices
  pr_long <- data.frame(year = rep(pr$year, ncol(pr) - 1),
                        species = factor(rep(names(pr)[-1], each = nrow(pr)),
                                         levels = names(pr)[-1]),
                        price = unlist(pr[-1], use.names = FALSE))
  pr_long$group <- sim$settings$species$group[match(pr_long$species,
                                                    sim$settings$species$species)]
  both <- rbind(data.frame(year = sim$truth$year, species = sim$truth$species,
                           value = sim$truth$relative, what = "true biomass"),
                data.frame(year = sim$nominal$year, species = sim$nominal$species,
                           value = sim$nominal$relative, what = "nominal CPUE"))
  list(
    env = ggplot(env, aes(x = .data$x, y = .data$y)) +
      geom_raster(aes(fill = .data$depth)) +
      geom_point(data = reef_cells, shape = 15, size = 0.7, colour = "grey15") +
      scale_fill_viridis_c(option = "mako", direction = -1) + coord_equal() +
      labs(x = "x (offshore ->)", y = "y", fill = "Depth (m)",
           subtitle = "dark squares: reef patches") + theme_bw(),
    sst = ggplot(env, aes(x = .data$x, y = .data$y, fill = .data$sst)) +
      geom_raster() + scale_fill_viridis_c(option = "inferno") + coord_equal() +
      labs(x = "x", y = "y", fill = "SST") + theme_bw(),
    density = ggplot(sim$density, aes(x = .data$x, y = .data$y, fill = log1p(.data$density))) +
      geom_raster() + facet_wrap(~species, ncol = 5) +
      scale_fill_viridis_c(option = "turbo") + coord_equal() +
      labs(x = NULL, y = NULL, fill = "log(1 + density)") + theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            strip.text = element_text(size = 7)),
    tactics = ggplot(data.frame(year = tr$year, tactic = factor(tr$tactic, levels = lev)),
                     aes(x = .data$year, fill = .data$tactic)) +
      geom_bar(position = "fill") + scale_fill_manual(values = cluster_colours(length(lev))) +
      labs(x = "Year", y = "Share of trips (true tactic)", fill = "Tactic") + theme_bw(),
    effort = ggplot(data.frame(year = tr$year, type = tr$vessel_type),
                    aes(x = .data$year, fill = .data$type)) +
      geom_bar() + labs(x = "Year", y = "Trips", fill = "Vessel type") + theme_bw(),
    prices = ggplot(pr_long, aes(x = .data$year, y = .data$price, group = .data$species,
                                 colour = .data$group)) +
      geom_step() + labs(x = "Year", y = "Price", colour = NULL) + theme_bw(),
    cpue = ggplot(both, aes(x = .data$year, y = .data$value, colour = .data$what)) +
      geom_line() + geom_hline(yintercept = 1, linetype = 3) +
      facet_wrap(~species, ncol = 5, scales = "free_y") +
      scale_colour_manual(values = c(`true biomass` = "black", `nominal CPUE` = "#d95f0e")) +
      labs(x = "Year", y = "Relative to mean", colour = NULL) + theme_bw() +
      theme(strip.text = element_text(size = 7), legend.position = "bottom")
  )
}
