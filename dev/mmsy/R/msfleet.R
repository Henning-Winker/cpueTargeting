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
                      price = NULL, cost_share = 0.3,
                      cost_ground = c(shallow = 1, deep = 1.5),
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

  # costs per unit effort: cost_share of mean unfished revenue rate, by ground
  rev0 <- as.vector(Q %*% (price * sp$B0))
  cost <- cost_share * mean(rev0) * cost_ground[met$ground]

  structure(list(
    species = sp$species, group = sp$group, r = sp$r, K = sp$B0, price = price,
    metiers = met, Q = Q, avail = avail, month_w = month_w,
    segments = segments, home = home_mat, loyalty = loyalty[names(segments)],
    cost = setNames(cost, met$name), rev0 = setNames(rev0, met$name),
    S = ns, M = M, G = G), class = "ms_system")
}

# -----------------------------------------------------------------
# 2. Effort allocation (preference function)
# -----------------------------------------------------------------

#' Allocate segment effort to metiers
#'
#' share_gm = loyalty_g * P(m | home_g) + (1 - loyalty_g) * P(m | all)
#' P(m | set) proportional to open_m * w0_m * pmax(profit_m, 0)^gamma
#' gamma = 0: profit-blind (equal shares, or w0); larger gamma: stronger
#' preference for the most profitable metier. Returns an n x M matrix of
#' effort by metier.
#' @param profit n x M expected profit per unit effort
#' @param E_seg  n x G effort by segment
#' @param open   n x M or length-M 0/1 (closures)
allocate_effort <- function(sys, profit, E_seg, gamma, open = 1, w0 = NULL) {
  n <- nrow(profit); M <- sys$M
  pr <- pmax(profit, 0)
  util <- if (gamma == 0) matrix(1, n, M) else pr^gamma
  if (!is.null(w0)) util <- util * (if (is.matrix(w0)) w0 else matrix(w0, n, M, byrow = TRUE))
  util <- util * (if (is.matrix(open)) open else matrix(open, n, M, byrow = TRUE))
  E_m <- matrix(0, n, M)
  all_set <- rep(1, M)
  for (g in seq_len(sys$G)) {
    for (part in c("home", "free")) {
      lam <- if (part == "home") sys$loyalty[g] else 1 - sys$loyalty[g]
      if (lam == 0) next
      acc <- if (part == "home") sys$home[g, ] else all_set
      u <- util * matrix(acc, n, M, byrow = TRUE)
      tot <- rowSums(u)
      p <- u / ifelse(tot > 0, tot, 1)          # rows with no profitable option: no effort
      E_m <- E_m + p * (lam * E_seg[, g])
    }
  }
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
#' @param phi     entry/exit rate: effort_t+1 = min(cap, E_t * exp(phi * profit margin));
#'   0 = effort always at the cap
#' @param open    metier access (length M, or function(year) -> length M)
#' @param hcr     optional function(state) -> n x G multiplier on caps
#' @param B_init  n x S initial biomass (default K)
#' @param proc_cv process error CV on growth; ncv year draws per row
#' @param keep    "final" (equilibrium summaries) or "all" (trajectories)
ms_simulate <- function(sys, caps, years = 100, gamma = 1.5,
                        behaviour = c("dynamic", "fixed"),
                        info = c("current", "lagged"),
                        phi = 0, open = 1, hcr = NULL, B_init = NULL,
                        proc_cv = 0, rho = 0, seed = NULL, keep = c("final", "all"),
                        avg_years = 10) {
  behaviour <- match.arg(behaviour); info <- match.arg(info); keep <- match.arg(keep)
  if (!is.null(seed)) set.seed(seed)
  caps <- as.matrix(caps); n <- nrow(caps); S <- sys$S; M <- sys$M; G <- sys$G
  K <- matrix(sys$K, n, S, byrow = TRUE); r <- matrix(sys$r, n, S, byrow = TRUE)
  price <- sys$price; Q <- sys$Q; tQ <- t(Q)
  cost_m <- matrix(sys$cost, n, M, byrow = TRUE)
  B <- B_init %||% K
  E_seg <- caps
  open_fun <- if (is.function(open)) open else function(y) open
  w_fixed <- NULL
  if (behaviour == "fixed") {
    # shares frozen at the unfished allocation (per unit segment effort)
    v0 <- (K * matrix(price, n, S, byrow = TRUE)) %*% tQ - cost_m
    w_fixed <- v0  # used through allocate_effort with gamma at unfished profit
  }
  cpue_last <- (B * matrix(price, n, S, byrow = TRUE)) %*% tQ   # value per effort, last year
  eps <- 0
  sig <- if (proc_cv > 0) sqrt(log(1 + proc_cv^2)) else 0

  nkeep <- if (keep == "all") years else min(avg_years, years)
  Bk <- array(NA_real_, c(nkeep, n, S)); Ck <- Bk; Fk <- Bk
  Ek <- array(NA_real_, c(nkeep, n, M)); Rk <- matrix(NA_real_, nkeep, n); Pk <- Rk
  Esk <- array(NA_real_, c(nkeep, n, G))

  for (y in seq_len(years)) {
    op <- open_fun(y)
    if (!is.null(hcr)) {
      E_cap_y <- caps * hcr(list(B = B, K = K, year = y, sys = sys))
    } else E_cap_y <- caps
    E_seg <- pmin(E_seg, E_cap_y); if (phi == 0) E_seg <- E_cap_y
    Bm <- B
    C_y <- matrix(0, n, S); F_y <- C_y; Em_y <- matrix(0, n, M); rev_m <- matrix(0, n, M)
    for (mo in 1:12) {
      av <- matrix(sys$avail[mo, ], n, S, byrow = TRUE)
      Bav <- Bm * av
      if (behaviour == "fixed") {
        E_m <- allocate_effort(sys, w_fixed, E_seg * sys$month_w[mo], gamma, op)
      } else {
        exp_val <- if (info == "current") (Bav * matrix(price, n, S, byrow = TRUE)) %*% tQ
                   else cpue_last
        E_m <- allocate_effort(sys, exp_val - cost_m, E_seg * sys$month_w[mo], gamma, op)
      }
      Fm <- (E_m %*% Q) * av                       # n x S fishing mortality this month
      Cm <- Bm * (1 - exp(-Fm))
      # catch by metier (proportional split of each species' catch by F share)
      share_m <- array(0, c(n, M))
      Fs <- pmax(Fm, 1e-12)
      for (m in seq_len(M)) {
        Fim <- outer(E_m[, m], Q[m, ]) * av
        rev_m[, m] <- rev_m[, m] + rowSums((Cm * Fim / Fs) %*% diag(price, S))
      }
      Bm <- Bm - Cm
      C_y <- C_y + Cm; F_y <- F_y + Fm; Em_y <- Em_y + E_m
    }
    cpue_last <- ifelse(Em_y > 0, rev_m / pmax(Em_y, 1e-12), cpue_last)
    revenue <- as.vector(C_y %*% price)
    costs <- rowSums(Em_y * cost_m)
    # growth on start-of-year biomass with process error
    if (sig > 0) {
      z <- matrix(rnorm(n * S), n, S)
      eps <- if (y == 1 || rho == 0) z else rho * eps + sqrt(1 - rho^2) * z
      pe <- exp(sig * eps - sig^2 / 2)
    } else pe <- 1
    keep_idx <- if (keep == "all") y else y - (years - nkeep)
    if (keep_idx >= 1) {
      Bk[keep_idx, , ] <- B; Ck[keep_idx, , ] <- C_y; Fk[keep_idx, , ] <- F_y; Ek[keep_idx, , ] <- Em_y
      Rk[keep_idx, ] <- revenue; Pk[keep_idx, ] <- revenue - costs; Esk[keep_idx, , ] <- E_seg
    }
    B <- pmax(B + r * B * (1 - B / K) * pe - C_y, 1e-6 * K)
    # entry / exit of effort within the cap
    if (phi > 0) {
      margin <- ifelse(revenue > 0, (revenue - costs) / revenue, -1)
      E_seg <- pmin(E_cap_y, E_seg * exp(phi * margin))
    }
  }
  out <- list(B = Bk, C = Ck, F = Fk, E_metier = Ek, E_seg = Esk, revenue = Rk, profit = Pk,
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
