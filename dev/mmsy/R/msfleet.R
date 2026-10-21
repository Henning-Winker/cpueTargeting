# =================================================================
# msfleet.R: multispecies fleet system (prototype)
#
# Schaefer stock pools + fleet segments + metiers (tactic x ground)
# with a preference (expected-profit) allocation of effort.
# All functions are vectorised over rows of a state matrix (n rows =
# scenarios x simulations), so equilibria for thousands of effort
# scenarios run in one call.
#
# Henning Winker, FAO NFIFM (prototype, dev-mmsy)
# =================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# -----------------------------------------------------------------
# 1. System definition from the logbook landscape
# -----------------------------------------------------------------

#' Build the system: species, metiers, catchability matrix, segments
#'
#' Metiers are tactic x ground. Demersal tactics (Reef, Sediment) fish
#' the reef or sediment cells of the shallow (< depth_cut) or deep
#' ground; the Pelagic tactic fishes all cells of a ground. The
#' catchability of species i in metier m is
#'   q[m, i] = q_scale * tac[tactic(m), i] * sum_c w_mc * K_ic / K_i
#' where w_mc is the distribution of sets over the cells of the ground
#' (proportional to unfished value ^ preference, as in simulate_logbook).
#' Biomass is pooled per species (fast mixing), so q[m, i] is the share
#' of species i a unit of effort in metier m meets.
ms_system <- function(L, depth_cut = 300, preference = 0,
                      segments = c(demersal = 1/3, pelagic = 1/3, mixed = 1/3),
                      home = list(demersal = c("Reef_shallow", "Reef_deep",
                                               "Sediment_shallow", "Sediment_deep"),
                                  pelagic  = c("Pelagic_shallow", "Pelagic_deep"),
                                  mixed    = NULL),
                      loyalty = c(demersal = 0.8, pelagic = 0.8, mixed = 0),
                      price = NULL, cost_share = 0.4, cost_mode = c("relative", "common"),
                      cost_ground = c(shallow = 1, deep = 1.5),
                      eff_home = NULL, eff_away = NULL, cost_seg = NULL, vessel_cv = NULL,
                      season = TRUE) {
  s   <- L$settings
  sp  <- L$species
  ns  <- nrow(sp)
  pref <- preference   # within-ground set concentration (0 = uniform; pooled biomass cannot deplete locally)
  price <- price %||% sp$price
  tac <- s$tactics[, sp$species, drop = FALSE]
  zone <- ifelse(L$depth >= depth_cut, "deep", "shallow")
  share <- sweep(L$K, 2, colSums(L$K), "/")             # cell share of each species' K

  met <- expand.grid(tactic = rownames(tac), ground = c("shallow", "deep"),
                     stringsAsFactors = FALSE)
  met$name <- paste(met$tactic, met$ground, sep = "_")
  M <- nrow(met)
  Q <- matrix(0, M, ns, dimnames = list(met$name, sp$species))
  met$n_cells <- NA
  for (m in seq_len(M)) {
    cells <- zone == met$ground[m]
    if (met$tactic[m] == "Reef")     cells <- cells & L$reef
    if (met$tactic[m] == "Sediment") cells <- cells & !L$reef
    met$n_cells[m] <- sum(cells)
    # unfished value per set in each cell, then set distribution w
    v <- as.vector((L$K[cells, , drop = FALSE] * rep(price, each = sum(cells))) %*% tac[met$tactic[m], ])
    w <- v^pref / sum(v^pref)
    Q[m, ] <- tac[met$tactic[m], ] * colSums(w * share[cells, , drop = FALSE])
  }

  # seasonal availability (months x species), as in simulate_logbook
  avail <- if (season) t(vapply(1:12, function(mo)
    (1 - sp$season_amp) + sp$season_amp * (1 + cos(2 * pi * (mo - sp$season_peak) / 12)) / 2,
    numeric(ns))) else matrix(1, 12, ns)
  month_w <- 1 + s$effort_season_amp * cos(2 * pi * (1:12 - 7) / 12)
  month_w <- month_w / sum(month_w)

  # segment access: home metiers (loyal effort) and all metiers (free effort)
  G <- length(segments)
  home_mat <- matrix(0, G, M, dimnames = list(names(segments), met$name))
  for (g in names(segments)) home_mat[g, home[[g]] %||% met$name] <- 1
  # fishing efficiency by segment and metier (mean of the vessel effects of the segment;
  # in simulate_logbook terms the type multiplier on q_vessel): multiplies catchability, so a day of a
  # specialist on its home metiers removes more fish than a day of a generalist
  # defaults: the type efficiencies of the logbook settings (1 if absent), on and
  # off the home metiers; daily cost proportional to efficiency (smaller, cheaper
  # vessels), so generalists stay viable at the specialists' break-even
  te <- s$type_efficiency %||% setNames(rep(1, length(segments)), names(segments))
  eff_home <- eff_home %||% te; eff_away <- eff_away %||% te; cost_seg <- cost_seg %||% te
  eff <- home_mat * eff_home[names(segments)] + (1 - home_mat) * eff_away[names(segments)]

  # costs per unit effort: cost_share of mean unfished revenue rate, by ground
  # cost_mode "common":   cost_share x mean unfished revenue rate, x ground factor
  #           "relative": cost_share x the metier's own unfished revenue rate, so every
  #                       metier breaks even at the same relative revenue loss
  cost_mode <- match.arg(cost_mode)
  rev0 <- as.vector(Q %*% (price * sp$B0 * colMeans(avail)))
  cost <- if (cost_mode == "common") cost_share * mean(rev0) * cost_ground[met$ground]
          else cost_share * rev0

  structure(list(
    species = sp$species, group = sp$group, r = sp$r, K = sp$B0, price = price,
    metiers = met, Q = Q, avail = avail, month_w = month_w,
    segments = segments, home = home_mat, loyalty = loyalty[names(segments)], eff = eff,
    cost_seg = cost_seg[names(segments)],   # daily cost of a segment relative to the specialists
    vessel_cv = vessel_cv %||% s$vessel_cv %||% 0.25,   # CV of vessel efficiency within a segment
    cost = setNames(cost, met$name), rev0 = setNames(rev0, met$name),
    S = ns, M = M, G = G), class = "ms_system")
}

# -----------------------------------------------------------------
# 2. Effort allocation (preference function)
# -----------------------------------------------------------------

#' Allocate segment effort to metiers
#'
#' share_gm = loyalty_g * P(m | home_g) + (1 - loyalty_g) * P(m | all)
#' P(m | set) proportional to open_m * pmax(profit_gm, 0)^gamma, with the
#' segment's own profit per unit (nominal) effort
#'   profit_gm = eff_gm * value_m - cost_m
#' gamma = 0: profit-blind (equal shares); larger gamma: stronger preference
#' for the most profitable metier. Returns an n x M matrix of nominal effort
#' by metier, with attributes "by_seg" (n x G x M nominal effort) and
#' "effective" (n x M effort weighted by efficiency, which drives F).
#' @param value  n x M expected revenue per unit of effective effort
#' @param cost   n x M cost per unit of nominal effort
#' @param emult  n x G mean efficiency of the active vessels relative to the segment mean
#' @param E_seg  n x G effort by segment
#' @param open   n x M or length-M 0/1 (closures)
allocate_effort <- function(sys, value, E_seg, gamma, open = 1, cost = 0, by_seg = TRUE, emult = NULL,
                            cseg = NULL) {
  n <- nrow(value); M <- sys$M
  cost <- if (is.matrix(cost)) cost else matrix(cost, n, M, byrow = TRUE)
  op <- if (is.matrix(open)) open else matrix(open, n, M, byrow = TRUE)
  E_m <- matrix(0, n, M); E_eff <- E_m
  E_gm <- array(0, c(n, sys$G, M))
  if (is.null(cseg)) cseg <- matrix(sys$cost_seg, n, sys$G, byrow = TRUE)
  all_set <- rep(1, M)
  for (g in seq_len(sys$G)) {
    effg <- matrix(sys$eff[g, ], n, M, byrow = TRUE)
    if (!is.null(emult)) effg <- effg * emult[, g]     # mean efficiency of the active vessels
    util <- if (gamma == 0) matrix(1, n, M) else pmax(value * effg - cost * cseg[, g], 0)^gamma
    util <- util * op
    for (part in c("home", "free")) {
      lam <- if (part == "home") sys$loyalty[g] else 1 - sys$loyalty[g]
      if (lam == 0) next
      acc <- if (part == "home") sys$home[g, ] else all_set
      u <- util * matrix(acc, n, M, byrow = TRUE)
      tot <- rowSums(u)
      p <- u / ifelse(tot > 0, tot, 1)          # rows with no profitable option: no effort
      E_add <- p * (lam * E_seg[, g])
      E_m <- E_m + E_add; E_eff <- E_eff + E_add * effg
      E_gm[, g, ] <- E_gm[, g, ] + E_add
    }
  }
  attr(E_m, "by_seg") <- E_gm; attr(E_m, "effective") <- E_eff
  E_m
}

# -----------------------------------------------------------------
# 3. Dynamics
# -----------------------------------------------------------------

#' Simulate the system forward
#'
#' Annual Schaefer growth on start-of-year biomass; catches removed
#' monthly (Baranov-type, C = B (1 - exp(-F))). Each month every segment
#' allocates its effort among metiers by expected profit computed from
#' the biomass it perceives (`info`: "current" = true current biomass,
#' "lagged" = last year's realised CPUE by metier).
#'
#' @param caps    n x G matrix of annual effort caps by segment
#' @param gamma   preference exponent (0 = profit-blind)
#' @param behaviour "dynamic" (reallocate every month) or "fixed"
#'   (shares frozen at the unfished allocation: the classical shared-q model)
#' @param phi     entry/exit rate of each segment: E_t+1 = min(cap, E_t * exp(phi * margin_g)),
#'   margin_g = segment profit / segment revenue; 0 = effort always at the cap
#' @param open    metier access (length M, or function(year) -> length M)
#' @param hcr     optional function(state) -> n x G multiplier on caps
#' @param hcr_lag years between the data a rule sees and the year its advice applies
#'   (default 2: data of year y - 2, assessed in y - 1, advice for y). Before the
#'   projection the data years are taken to be the initial state.
#' @param B_init  n x S initial biomass (default K)
#' @param E_init  initial effort by segment (vector of length G or n x G; default the caps)
#' @param proc_cv process error CV at the median r (JABBA form: the whole
#'   process equation is multiplied by exp(eta), eta ~ N(-sigma^2/2, sigma^2))
#' @param rho     AR(1) autocorrelation of the process error in time
#' @param proc_r_power the process CV of species i is proc_cv * (r_i / median r)^power:
#'   more productive species are more variable (0 = same CV for all)
#' @param proc_rho_sp correlation of the process error across species (shared environment)
#' @param impl_cv,impl_bias implementation error on the effort caps by segment:
#'   realised cap = cap * exp(bias + N(-s^2/2, s^2)), s from impl_cv
#' @param crn     common random numbers: rows i and i + crn share the same process
#'   (and rule) error draws, e.g. crn = nsim when the rows are nsim simulations
#'   repeated for each of several control settings
#' @param cost_seg daily cost by segment (length G or n x G), overrides sys$cost_seg
#' @param capacity fleet days by segment (registered vessels x days); default the caps.
#'   With vessel heterogeneity (sys$vessel_cv > 0) and phi > 0, the days used move
#'   towards capacity x the share of vessels that break even (lognormal vessel
#'   efficiency, mean 1), and the active vessels' mean efficiency is the mean of
#'   that upper tail. phi is then the annual adjustment rate (0-1).
#' @param proc_cv process error CV on growth; ncv year draws per row
#' @param keep    "final" (equilibrium summaries) or "all" (trajectories)
ms_simulate <- function(sys, caps, years = 100, gamma = 1.5,
                        behaviour = c("dynamic", "fixed"),
                        info = c("current", "lagged"),
                        phi = 0.5, open = 1, hcr = NULL, hcr_lag = 2, B_init = NULL, E_init = NULL, capacity = NULL,
                        cost_seg = NULL,
                        proc_cv = 0, rho = 0, proc_r_power = 1, proc_rho_sp = 0,
                        impl_cv = 0, impl_bias = 0,
                        seed = NULL, crn = NULL, keep = c("final", "all"),
                        avg_years = 10) {
  behaviour <- match.arg(behaviour); info <- match.arg(info); keep <- match.arg(keep)
  if (!is.null(seed)) set.seed(seed)
  caps <- as.matrix(caps); n <- nrow(caps); S <- sys$S; M <- sys$M; G <- sys$G
  K <- matrix(sys$K, n, S, byrow = TRUE); r <- matrix(sys$r, n, S, byrow = TRUE)
  price <- sys$price; Q <- sys$Q; tQ <- t(Q)
  cost_m <- matrix(sys$cost, n, M, byrow = TRUE)
  B <- B_init %||% K
  E_seg <- if (is.null(E_init)) caps else pmin(caps, matrix(E_init, n, G, byrow = is.null(dim(E_init))))
  open_fun <- if (is.function(open)) open else function(y) open
  Bbuf <- replicate(hcr_lag + 1, B, simplify = FALSE)      # start-of-year biomass, most recent first
  capacity <- if (is.null(capacity)) caps else matrix(capacity, n, G, byrow = is.null(dim(capacity)))
  emult <- matrix(1, n, G)
  idx <- if (is.null(crn)) seq_len(n) else (seq_len(n) - 1) %% crn + 1
  nd <- if (is.null(crn)) n else crn
  cseg <- if (is.null(cost_seg)) matrix(sys$cost_seg, n, G, byrow = TRUE) else
    matrix(cost_seg, n, G, byrow = is.null(dim(cost_seg)))
  sv <- sys$vessel_cv; sdv <- if (sv > 0) sqrt(log(1 + sv^2)) else 0; muv <- -sdv^2 / 2
  w_fixed <- NULL
  if (behaviour == "fixed") {
    # shares frozen at the unfished allocation (per unit segment effort)
    w_fixed <- (K * matrix(price, n, S, byrow = TRUE)) %*% tQ   # unfished revenue rate, frozen
  }
  cpue_last <- (B * matrix(price, n, S, byrow = TRUE)) %*% tQ   # value per effort, last year
  eps <- 0
  cv_i <- proc_cv * (sys$r / stats::median(sys$r))^proc_r_power
  sig_i <- matrix(sqrt(log(1 + cv_i^2)), n, S, byrow = TRUE)          # by species
  sig <- proc_cv
  sig_impl <- if (impl_cv > 0) sqrt(log(1 + impl_cv^2)) else 0

  nkeep <- if (keep == "all") years else min(avg_years, years)
  Bk <- array(NA_real_, c(nkeep, n, S)); Ck <- Bk; Fk <- Bk
  Ek <- array(NA_real_, c(nkeep, n, M)); Rk <- matrix(NA_real_, nkeep, n); Pk <- Rk
  Esk <- array(NA_real_, c(nkeep, n, G)); Egk <- array(NA_real_, c(nkeep, n, G, M))

  for (y in seq_len(years)) {
    op <- open_fun(y)
    opm <- matrix(if (is.null(dim(op))) rep_len(op, M) else op, n, M, byrow = is.null(dim(op)))
    E_cap_y <- caps
    if (!is.null(hcr)) {
      # a rule returns multipliers on the caps (n x G), or list(mult, open) where
      # open (n x M, 0-1) is the share of each metier's ground/season left open
      Bbuf <- c(list(B), Bbuf)[seq_len(hcr_lag + 1)]
      ctl <- hcr(list(B = Bbuf[[hcr_lag + 1]], K = K, year = y, sys = sys, idx = idx, nd = nd, lag = hcr_lag))
      if (is.list(ctl)) {
        if (!is.null(ctl$mult)) E_cap_y <- caps * ctl$mult
        if (!is.null(ctl$open)) opm <- opm * ctl$open
      } else E_cap_y <- caps * ctl
    }
    if (sig_impl > 0 || impl_bias != 0) {                    # implementation error
      zi <- matrix(rnorm(nd * G), nd, G)[idx, , drop = FALSE]
      E_cap_y <- E_cap_y * exp(impl_bias + sig_impl * zi - sig_impl^2 / 2)
    }
    E_seg <- pmin(E_seg, E_cap_y); if (phi == 0) E_seg <- E_cap_y
    Bm <- B
    C_y <- matrix(0, n, S); F_y <- C_y; Em_y <- matrix(0, n, M); rev_m <- matrix(0, n, M); Eg_y <- array(0, c(n, G, M)); Eeff_y <- Em_y
    for (mo in 1:12) {
      av <- matrix(sys$avail[mo, ], n, S, byrow = TRUE)
      Bav <- Bm * av
      if (behaviour == "fixed") {
        E_m <- allocate_effort(sys, w_fixed * opm, E_seg * sys$month_w[mo], gamma, opm > 0, cost = cost_m, emult = emult, cseg = cseg)
      } else {
        exp_val <- if (info == "current") (Bav * matrix(price, n, S, byrow = TRUE)) %*% tQ
                   else cpue_last
        # a partly closed metier fishes only its open share: value per day scales with it
        E_m <- allocate_effort(sys, exp_val * (if (info == "current") opm else 1), E_seg * sys$month_w[mo],
                               gamma, opm > 0, cost = cost_m, emult = emult, cseg = cseg)
      }
      E_eff <- attr(E_m, "effective") * opm         # effective effort on the open share
      Fm <- (E_eff %*% Q) * av                     # n x S fishing mortality this month
      Cm <- Bm * (1 - exp(-Fm))
      # catch by metier (proportional split of each species' catch by F share)
      share_m <- array(0, c(n, M))
      Fs <- pmax(Fm, 1e-12)
      for (m in seq_len(M)) {
        Fim <- outer(E_eff[, m], Q[m, ]) * av
        rev_m[, m] <- rev_m[, m] + rowSums((Cm * Fim / Fs) %*% diag(price, S))
      }
      Bm <- Bm - Cm
      C_y <- C_y + Cm; F_y <- F_y + Fm; Em_y <- Em_y + E_m; Eg_y <- Eg_y + attr(E_m, "by_seg"); Eeff_y <- Eeff_y + E_eff
    }
    cpue_last <- ifelse(Eeff_y > 0, rev_m / pmax(Eeff_y, 1e-12), cpue_last)   # per effective effort
    revenue <- as.vector(C_y %*% price)
    costs <- rowSums(matrix(sapply(seq_len(G), function(g)
      rowSums(matrix(Eg_y[, g, ], n, M) * cost_m) * cseg[, g]), n, G))
    # growth on start-of-year biomass with process error
    if (sig > 0) {
      z <- matrix(rnorm(nd * S), nd, S)
      if (proc_rho_sp > 0) z <- sqrt(1 - proc_rho_sp) * z + sqrt(proc_rho_sp) * stats::rnorm(nd)
      z <- z[idx, , drop = FALSE]
      eps <- if (y == 1 || rho == 0) z else rho * eps + sqrt(1 - rho^2) * z
      pe <- exp(sig_i * eps - sig_i^2 / 2)
    } else pe <- 1
    keep_idx <- if (keep == "all") y else y - (years - nkeep)
    if (keep_idx >= 1) {
      Bk[keep_idx, , ] <- B; Ck[keep_idx, , ] <- C_y; Fk[keep_idx, , ] <- F_y; Ek[keep_idx, , ] <- Em_y; Egk[keep_idx, , , ] <- Eg_y
      Rk[keep_idx, ] <- revenue; Pk[keep_idx, ] <- revenue - costs; Esk[keep_idx, , ] <- E_seg
    }
    B <- pmax((B + r * B * (1 - B / K) - C_y) * pe, 1e-6 * K)     # JABBA-type process error
    # entry / exit below the cap
    if (phi > 0) {
      for (g in seq_len(G)) {
        Egm <- matrix(Eg_y[, g, ], n, M); tot <- rowSums(Egm)
        acc <- if (sys$loyalty[g] < 1) rep(TRUE, M) else sys$home[g, ] > 0
        # revenue and cost per day of a vessel of mean efficiency: allocation-weighted
        # when the segment fishes, otherwise its best open metier
        val <- cpue_last * matrix(sys$eff[g, ], n, M, byrow = TRUE) * (opm > 0)   # realised: already net of closure
        cst <- cost_m * cseg[, g]
        net <- val - cst; net[, !acc] <- -Inf
        best <- max.col(net, ties.method = "first")
        rev_day <- ifelse(tot > 0, rowSums(Egm * val) / pmax(tot, 1e-12), val[cbind(seq_len(n), best)])
        cost_day <- ifelse(tot > 0, rowSums(Egm * cst) / pmax(tot, 1e-12), cst[cbind(seq_len(n), best)])
        if (sdv > 0) {
          # heterogeneous vessels: those with efficiency e > e* = cost / revenue break even
          lest <- log(pmax(cost_day, 1e-12) / pmax(rev_day, 1e-12))
          pg <- 1 - stats::pnorm((lest - muv) / sdv)
          emult[, g] <- pmin(stats::pnorm((muv + sdv^2 - lest) / sdv) / pmax(pg, 1e-9), exp(3 * sdv))
          target <- pmin(E_cap_y[, g], capacity[, g] * pg)
          E_seg[, g] <- E_seg[, g] + min(phi, 1) * (target - E_seg[, g])
        } else {
          margin <- ifelse(rev_day > 0, (rev_day - cost_day) / rev_day, -1)
          E_seg[, g] <- pmin(E_cap_y[, g], pmax(E_seg[, g], 1e-4 * E_cap_y[, g]) * exp(phi * margin))
        }
      }
    }
  }
  out <- list(B = Bk, C = Ck, F = Fk, E_metier = Ek, E_seg_metier = Egk, E_seg = Esk, revenue = Rk, profit = Pk,
              caps = caps, sys = sys,
              settings = list(gamma = gamma, behaviour = behaviour, info = info, phi = phi,
                              years = years, proc_cv = proc_cv))
  class(out) <- "ms_sim"
  out
}

#' Summaries over the kept years (equilibrium when keep = "final")
ms_summary <- function(x, Bpa = 0.5, Blim = 0.25) {
  sys <- x$sys
  Bmsy <- sys$K / 2
  Bbar <- apply(x$B, c(2, 3), mean)               # n x S
  Cbar <- apply(x$C, c(2, 3), mean)
  BBmsy <- sweep(Bbar, 2, Bmsy, "/")
  colnames(BBmsy) <- colnames(Cbar) <- sys$species
  data.frame(
    yield = rowSums(Cbar), revenue = colMeans(x$revenue), profit = colMeans(x$profit),
    E_total = rowSums(apply(x$E_seg, c(2, 3), mean)),
    min_BBmsy = apply(BBmsy, 1, min),
    worst = sys$species[apply(BBmsy, 1, which.min)],
    max_FFmsy = apply(sweep(apply(x$F, c(2, 3), mean), 2, sys$r / 2, "/"), 1, max),
    n_below_Bpa = rowSums(BBmsy < Bpa), n_below_Blim = rowSums(BBmsy < Blim),
    check.names = FALSE) -> out
  FFmsy <- sweep(apply(x$F, c(2, 3), mean), 2, sys$r / 2, "/"); colnames(FFmsy) <- sys$species
  attr(out, "BBmsy") <- BBmsy; attr(out, "FFmsy") <- FFmsy; attr(out, "catch") <- Cbar
  out
}

# -----------------------------------------------------------------
# 5. Constrained optimisation: MMSY / MMEY with a Bpa safeguard
# -----------------------------------------------------------------

#' Maximise an equilibrium objective over segment effort caps
#'
#' Finds the effort caps that maximise total yield (or revenue, or profit)
#' at equilibrium subject to every constrained stock staying at or above
#' Bpa: B_i / Bmsy_i >= Bpa for all i in `species`.
#'
#' The search is a vectorised zooming grid in log1p(effort) space: a grid of
#' n^d candidates is run in one call, the best feasible one is kept, and the
#' box shrinks around it. This is robust to the kinks and small cycles of a
#' behavioural fleet, which defeat gradient-based optimisers.
#'
#' @param levers  segments whose caps are optimised, or "total" for one global
#'   cap split in proportion `ratio` across all segments
#' @param fixed   named caps for segments not in `levers` (e.g. the free
#'   agent: mixed = 40)
#' @param ratio   split of the global cap when levers = "total"
#' @param upper   upper bound of each cap
#' @param species names of the stocks under the safeguard (default all)
#' @param min_share viability floor: every segment's cap must be at least this
#'   share of the total cap (e.g. 0.2)
#' @param ...     passed to ms_simulate() (gamma, behaviour, phi, open, ...)
#' @return list: caps, summary row, B/Bmsy and F/Fmsy, effort by segment and
#'   metier, the binding stock(s), and the search history
ms_optimise <- function(sys, objective = c("yield", "revenue", "profit"),
                        levers = names(sys$segments), fixed = NULL,
                        ratio = rep(1 / sys$G, sys$G), upper = 3000,
                        Bpa = 0.5, species = sys$species, min_share = 0, n = NULL, zoom = 5,
                        years = 150, ...) {
  objective <- match.arg(objective)
  seg <- names(sys$segments)
  d <- if (identical(levers, "total")) 1 else length(levers)
  n <- n %||% c(41, 19, 11)[min(d, 3)]
  lo <- rep(0, d); hi <- rep(log1p(upper), d)
  to_caps <- function(X) {                      # X: k x d in log1p space
    E <- expm1(X)
    if (identical(levers, "total")) return(outer(E[, 1], ratio))
    caps <- matrix(0, nrow(E), sys$G, dimnames = list(NULL, seg))
    caps[, levers] <- E
    for (g in names(fixed)) caps[, g] <- fixed[[g]]
    caps
  }
  isp <- match(species, sys$species)
  best <- NULL; hist <- list()
  for (z in seq_len(zoom)) {
    axes <- lapply(seq_len(d), function(k) seq(lo[k], hi[k], length.out = n))
    X <- as.matrix(expand.grid(axes))
    caps <- to_caps(X)
    sm <- ms_summary(ms_simulate(sys, caps, years = years, ...), Bpa = Bpa)
    BB <- attr(sm, "BBmsy")
    feas <- apply(BB[, isp, drop = FALSE], 1, min) >= Bpa &
      apply(caps, 1, min) >= min_share * rowSums(caps)
    val <- sm[[objective]]; val[!feas] <- -Inf
    i <- which.max(val)
    hist[[z]] <- data.frame(zoom = z, caps, value = sm[[objective]], feasible = feas)
    if (!is.finite(val[i])) break             # nothing feasible (e.g. free agent too large)
    best <- list(X = X[i, ], caps = caps[i, ], summary = sm[i, ], BBmsy = BB[i, ],
                 FFmsy = attr(sm, "FFmsy")[i, ])
    step <- (hi - lo) / (n - 1)
    lo <- pmax(0, X[i, ] - 2 * step); hi <- pmin(log1p(upper), X[i, ] + 2 * step)
  }
  if (is.null(best)) return(list(feasible = FALSE, history = do.call(rbind, hist)))
  # rerun the optimum to keep effort by segment and metier
  x <- ms_simulate(sys, matrix(best$caps, 1), years = years, ...)
  best$E_seg <- setNames(colMeans(x$E_seg[, 1, , drop = FALSE][, 1, ]), seg)
  best$E_seg_metier <- matrix(apply(x$E_seg_metier[, 1, , , drop = FALSE], c(3, 4), mean),
                              sys$G, sys$M, dimnames = list(seg, sys$metiers$name))
  best$binding <- sys$species[isp][abs(best$BBmsy[isp] - Bpa) < 0.05]
  best$feasible <- TRUE; best$objective <- objective
  best$history <- do.call(rbind, hist)
  best
}

#' Open-access (status quo) equilibrium: caps non-binding, effort set by profit
ms_open_access <- function(sys, years = 500, cap = 1000, ...) {
  x <- ms_simulate(sys, matrix(rep_len(cap, sys$G), 1, sys$G), years = years, ...)
  s <- ms_summary(x)
  s$E_seg <- NULL
  list(E_seg = setNames(colMeans(x$E_seg[, 1, , drop = FALSE][, 1, ]), names(sys$segments)),
       B = colMeans(x$B[, 1, ]), summary = s, BBmsy = attr(s, "BBmsy")[1, ], FFmsy = attr(s, "FFmsy")[1, ])
}

#' Status quo: capacity-limited equilibrium. Each segment has `capacity` days
#' (vector or scalar); profit decides how many are used (economic exit).
ms_status_quo <- function(sys, capacity = 1000, ...) ms_open_access(sys, cap = capacity, ...)

# -----------------------------------------------------------------
# 6. Tuning: fleet composition as an open-access equilibrium
# -----------------------------------------------------------------

#' Tune daily costs so that a target fleet composition is an open-access equilibrium
#'
#' At open access each segment's active days are the share of its vessels
#' that break even (heterogeneous vessel efficiency), so the composition of
#' the active fleet is set by the daily costs. The G log cost multipliers are
#' solved so that G residuals vanish: the shares of days of the first G - 1
#' segments and one status target for the stocks,
#'   r = (log share_g / target_g, g < G;  log status / target).
#' Newton steps with a forward-difference Jacobian; the base point and the G
#' perturbed points run in one vectorised call (warm-started).
#'
#' @param shares  target share of active days by segment (sums to 1)
#' @param status  c(median = x): median B/Bmsy; c(min = x): the weakest stock;
#'   c(E_total = x): total active days
#' @param cap     fleet capacity (registered days) by segment
#' @param ...     passed to ms_simulate() (gamma, ...)
#' @return the tuned system (cost_seg replaced), with attributes "tuning"
#'   (history), "status_quo" (the open-access equilibrium) and "converged"
ms_tune <- function(sys, shares = c(demersal = 0.35, pelagic = 0.35, mixed = 0.30),
                    status = c(median = 0.9), cap = 1000, max_iter = 25, tol = 0.02,
                    h = 0.05, max_step = 0.5, years_first = 400, years = 250, avg_years = 30,
                    verbose = TRUE, ...) {
  G <- sys$G; seg <- names(sys$segments)
  shares <- shares[seg] / sum(shares)
  stat_fun <- switch(names(status),
    median  = function(BB, E) apply(BB, 1, stats::median),
    min     = function(BB, E) apply(BB, 1, min),
    E_total = function(BB, E) rowSums(E),
    stop("status must be named median, min or E_total"))
  resid <- function(x) {
    E <- apply(x$E_seg, c(2, 3), mean); BB <- attr(ms_summary(x), "BBmsy")
    sh <- E / rowSums(E)
    cbind(log(pmax(sh[, -G, drop = FALSE], 1e-4) / matrix(shares[-G], nrow(E), G - 1, byrow = TRUE)),
          log(pmax(stat_fun(BB, E), 1e-4) / status))
  }
  lc <- log(sys$cost_seg); B0 <- NULL; E0 <- NULL; hist <- list(); conv <- FALSE
  best <- list(err = Inf, lc = lc)
  for (it in seq_len(max_iter)) {
    LC <- rbind(lc, t(lc + diag(h, G)))                    # base + G perturbations
    nr <- nrow(LC)
    x <- ms_simulate(sys, matrix(cap, nr, G, byrow = TRUE), years = if (it == 1) years_first else years,
                     cost_seg = exp(LC), capacity = matrix(cap, nr, G, byrow = TRUE),
                     B_init = if (is.null(B0)) NULL else B0[rep(1, nr), , drop = FALSE],
                     E_init = if (is.null(E0)) NULL else E0[rep(1, nr), , drop = FALSE],
                     avg_years = avg_years, ...)
    R <- resid(x); r0 <- R[1, ]
    E <- apply(x$E_seg, c(2, 3), mean)
    hist[[it]] <- data.frame(iter = it, t(setNames(round(exp(lc), 4), paste0("cost_", seg))),
                             t(setNames(round(E[1, ] / sum(E[1, ]), 3), paste0("share_", seg))),
                             status = round(exp(r0[G]) * status, 3), check.names = FALSE)
    if (verbose) cat(sprintf("iter %2d | costs %s | shares %s | %s %.3f\n", it,
                             paste(sprintf("%.3f", exp(lc)), collapse = "/"),
                             paste(sprintf("%.3f", E[1, ] / sum(E[1, ])), collapse = "/"),
                             names(status), exp(r0[G]) * status))
    if (max(abs(r0)) < best$err) best <- list(err = max(abs(r0)), lc = lc)
    if (max(abs(r0)) < tol) { conv <- TRUE; break }
    J <- t((R[-1, , drop = FALSE] - matrix(r0, G, G, byrow = TRUE)) / h)   # d r / d log c
    step <- tryCatch(-solve(J, r0), error = function(e) -0.3 * sign(r0))
    if (max(abs(step)) > max_step) step <- step * max_step / max(abs(step))
    lc <- lc + step
    B0 <- matrix(apply(x$B[, 1, , drop = FALSE], 3, mean), 1); E0 <- matrix(pmax(E[1, ], 1), 1)
  }
  sys$cost_seg <- setNames(exp(best$lc), seg)                 # best iterate
  attr(sys, "tuning") <- do.call(rbind, hist)
  attr(sys, "max_error") <- best$err
  attr(sys, "status_quo") <- ms_open_access(sys, cap = cap, avg_years = avg_years, ...)
  attr(sys, "converged") <- conv
  sys
}

# -----------------------------------------------------------------
# 7. Feedback control rules
# -----------------------------------------------------------------

#' Segment effort rule with an optional choke-stock closure
#'
#' Each year the assessed B/Bmsy of every stock is the true value times a
#' lognormal assessment error (CV `assess_cv`, AR(1) `assess_rho`). Each
#' segment's days are a fraction of its reference days (the caps passed to
#' ms_simulate(), e.g. the status-quo days): the fraction ramps from 0 at
#' Blim to 1 at Btrig (both relative to Bmsy) of the most depleted stock the
#' segment can catch substantially (any usable metier has >= `own_q` of the
#' stock's highest catchability).
#'
#' Stocks in `close` are protected by a closure instead: the open share of
#' the metiers that drive them (catchability >= `close_q` of the stock's
#' highest) follows the same ramp, and they are left out of the effort rule.
#' `max_change` limits the annual change of each segment's fraction.
#'
#' `Btrig` (effort rule) and `Btrig_close` (closure) may be vectors with one
#' value per row of the simulation, so a grid of triggers runs in one call.
#'
#' `assessed`: names of the stocks with an assessment (default all). The rule
#' sees only these; a stock in `close` must be assessed.
#'
#' Create a new rule for every simulation run (it keeps the last fraction and
#' the assessment error between years).
hcr_effort <- function(sys, Btrig = 0.8, Blim = 0.25, assess_cv = 0.2, assess_rho = 0,
                       own_q = 0.2, close = NULL, close_q = 0.5, max_change = Inf,
                       Btrig_close = Btrig, assessed = NULL) {
  S <- sys$S; G <- sys$G; M <- sys$M
  qmax <- apply(sys$Q, 2, max)
  own <- sapply(seq_len(G), function(g) {
    acc <- if (sys$loyalty[g] > 0) sys$home[g, ] > 0 else rep(TRUE, M)
    apply(sys$Q[acc, , drop = FALSE], 2, max) >= own_q * qmax })     # S x G
  if (!is.null(assessed)) {
    if (!all(close %in% assessed)) stop("a stock in `close` must be assessed")
    own[!sys$species %in% assessed, ] <- FALSE
  }
  ic <- match(close, sys$species)
  if (length(ic)) own[ic, ] <- FALSE
  drive <- sys$Q >= close_q * matrix(qmax, M, S, byrow = TRUE)          # M x S
  sig <- if (assess_cv > 0) sqrt(log(1 + assess_cv^2)) else 0
  ramp <- function(d, bt) pmin(1, pmax(0, (d - Blim) / (bt - Blim)))
  prev <- NULL; eps <- NULL
  rule <- function(st) {
    n <- nrow(st$B)
    if (st$year == 1) { prev <<- NULL; eps <<- NULL }
    nd <- st$nd %||% n; idx <- st$idx %||% seq_len(n)
    z <- matrix(stats::rnorm(nd * S), nd, S)[idx, , drop = FALSE]
    bt <- rep_len(Btrig, n); btc <- rep_len(Btrig_close, n)
    eps <<- if (is.null(eps) || assess_rho == 0) z else assess_rho * eps + sqrt(1 - assess_rho^2) * z
    bb <- st$B / (st$K / 2) * exp(sig * eps - sig^2 / 2)              # assessed B/Bmsy
    mult <- matrix(sapply(seq_len(G), function(g)
      if (any(own[, g])) ramp(apply(bb[, own[, g], drop = FALSE], 1, min), bt) else rep(1, n)), n, G)
    if (is.finite(max_change) && !is.null(prev))
      mult <- pmin(pmax(mult, prev * (1 - max_change)), prev * (1 + max_change))
    prev <<- mult
    if (!length(ic)) return(mult)
    opm <- matrix(1, n, M)
    for (i in ic) for (m in which(drive[, i])) opm[, m] <- pmin(opm[, m], ramp(bb[, i], btc))
    list(mult = mult, open = opm)
  }
  attr(rule, "own") <- own; attr(rule, "drive") <- drive
  rule
}

#' Performance of projection runs (named list of ms_sim with keep = "all")
perf_rows <- function(x, rows = NULL, years = NULL, Bpa = 0.5, Blim = 0.25, assessed = NULL) {
  sys <- x$sys; ny <- dim(x$B)[1]; rows <- rows %||% seq_len(dim(x$B)[2])
  yr <- years %||% seq(max(1, ny - 19), ny)
  bb <- sweep(x$B[yr, rows, , drop = FALSE], 3, sys$K / 2, "/")
  Ey <- apply(x$E_seg[yr, rows, , drop = FALSE], 1:2, sum)              # years x sims
  aav <- colMeans(abs(diff(Ey)) / pmax(Ey[-nrow(Ey), , drop = FALSE], 1e-9))
  pst <- apply(bb < Bpa, 3, mean)
  data.frame(
    yield = round(mean(apply(x$C[yr, rows, , drop = FALSE], 1:2, sum))),
    profit = round(mean(x$profit[yr, rows])),
    days = round(mean(Ey)),
    P_any_below_Blim = round(mean(apply(bb < Blim, 1:2, any)), 3),
    P_any_below_Bpa = round(mean(apply(bb < Bpa, 1:2, any)), 3),
    worst_P_below_Bpa = round(max(pst), 3),
    worst_stock = sys$species[which.max(pst)],
    worst_P_assessed = round(max(pst[if (is.null(assessed)) seq_along(pst) else sys$species %in% assessed]), 3),
    median_min_BBmsy = round(stats::median(apply(bb, 1:2, min)), 2),
    AAV_effort = round(stats::median(aav), 3))
}

ms_performance <- function(runs, years = NULL, Bpa = 0.5, Blim = 0.25) {
  do.call(rbind, lapply(names(runs), function(nm)
    cbind(run = nm, perf_rows(runs[[nm]], years = years, Bpa = Bpa, Blim = Blim))))
}

#' Tune the triggers of hcr_effort() to a risk criterion
#'
#' Runs a grid of effort-rule triggers (and closure triggers, if `close` is
#' given) in one vectorised call with common random numbers: every grid point
#' sees the same process and assessment errors. The tuned rule maximises the
#' objective among the grid points that meet the risk limit.
#'
#' @param triggers       effort-rule triggers (B/Bmsy)
#' @param close_triggers closure triggers; NULL = same as the effort trigger
#' @param risk           named limit on a perf_rows() column, e.g.
#'   c(worst_P_below_Bpa = 0.05) (every stock below Bpa at most 5% of the time)
#'   or c(worst_P_assessed = 0.05) (only the assessed stocks)
#' @param assessed       stocks with an assessment (default all)
#' @param objective      perf_rows() column to maximise ("yield" or "profit")
#' @param E_ref          reference days by segment (the caps the rule scales)
#' @param ...            passed to ms_simulate() (gamma, capacity, proc_cv, ...)
ms_tune_hcr <- function(sys, triggers = seq(0.5, 1.2, by = 0.05), close = NULL,
                        close_triggers = NULL, risk = c(worst_P_below_Bpa = 0.05),
                        objective = "yield", E_ref, B_init, E_init = E_ref,
                        nsim = 200, years = 50, eval_years = 31:50, seed = 1,
                        Blim = 0.25, assess_cv = 0.2, max_change = Inf, assessed = NULL, ...) {
  grid <- if (is.null(close) || is.null(close_triggers)) data.frame(Btrig = triggers, Btrig_close = triggers)
          else expand.grid(Btrig = triggers, Btrig_close = close_triggers)
  k <- nrow(grid); n <- k * nsim
  rep_rows <- function(v, m) matrix(v, n, m, byrow = TRUE)
  rule <- hcr_effort(sys, Btrig = rep(grid$Btrig, each = nsim), Btrig_close = rep(grid$Btrig_close, each = nsim),
                     Blim = Blim, close = close, assess_cv = assess_cv, max_change = max_change,
                     assessed = assessed)
  x <- ms_simulate(sys, rep_rows(E_ref, sys$G), years, B_init = rep_rows(B_init, sys$S),
                   E_init = rep_rows(E_init, sys$G), hcr = rule, crn = nsim, seed = seed, keep = "all", ...)
  tab <- cbind(grid, do.call(rbind, lapply(seq_len(k), function(j)
    perf_rows(x, rows = (j - 1) * nsim + seq_len(nsim), years = eval_years, assessed = assessed))))
  if (is.null(close)) tab$Btrig_close <- NULL
  ok <- tab[[names(risk)]] <= risk
  best <- if (any(ok)) tab[ok, ][which.max(tab[ok, objective]), ] else NULL
  list(table = tab, best = best, risk = risk, objective = objective)
}
