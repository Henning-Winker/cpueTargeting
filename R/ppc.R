#' Predictive check of the error distribution
#'
#' Compares the observed CPUE with replicate data sets simulated from the
#' fitted models, to check whether the error distribution (Tweedie,
#' delta-lognormal or delta-gamma) reproduces the share of zeros, the upper
#' tail and the shape of the positive catches, overall and by year (or any
#' other grouping).
#'
#' **How it works.** For each of `n_sim` replicates, the expected value of
#' every operation is computed from the fitted model and a new observation is
#' drawn from the fitted error distribution (Tweedie with the estimated power
#' and dispersion; for delta models a Bernoulli draw for presence and a
#' lognormal or gamma draw for the positive catch). With
#' `param_uncertainty = TRUE` (default), the model coefficients of each
#' replicate are first drawn from their approximate posterior (multivariate
#' normal with the Bayesian covariance matrix `Vp`), which makes this an
#' approximate posterior predictive check; otherwise it is a plug-in
#' (parametric bootstrap) check at the estimated values. Random vessel
#' effects are kept at their estimated values.
#'
#' Each replicate is treated as a complete data set: the summaries (share of
#' zeros, mean, upper quantiles of positive CPUE, histogram of log CPUE) are
#' computed per replicate, and the spread across replicates gives the
#' envelope (2.5-97.5%) against which the observed value is judged. The
#' `tail_prob` column is the share of replicates with a value at least as
#' large as the observed one (a posterior predictive p-value; values near 0
#' or 1 indicate misfit).
#'
#' **Reading the plots.**
#' * `zeros`: observed share of zero catches above the envelope means the
#'   model predicts too few zeros (e.g. a lognormal or Poisson-like error, or
#'   a missing targeting term in years when the fleet fished elsewhere);
#'   below means too many.
#' * `tails`: observed upper quantiles below the envelope mean the
#'   distribution's tail is too heavy (e.g. a negative binomial or a
#'   lognormal with large variance); above means too light.
#' * `distribution`: the share of operations in each bin of log10 CPUE of
#'   positive catches; systematic departures show a wrong shape (e.g. the
#'   wrong Tweedie power).
#'
#' Pooled over all operations, errors in opposite directions can cancel, so
#' the by-year summaries are usually more informative. Models without
#' targeting covariates often get the yearly share of zeros wrong in exactly
#' the years when targeting changes.
#'
#' @param fit Output of [fit_cpue()].
#' @param models Models to check (default: all).
#' @param n_sim Number of replicate data sets.
#' @param by Grouping for the summaries: `"year"` (default), any column of
#'   `fit$data` (e.g. `"FT"`, `"month"`) or `NULL` for all operations
#'   together.
#' @param param_uncertainty Logical; draw the coefficients of each replicate
#'   from their approximate posterior.
#' @param probs Upper quantiles of positive CPUE to compare.
#' @param bins Number of bins of the histogram of log10 positive CPUE.
#' @param seed Random seed.
#' @return An object of class `cpue_ppc`: a list with `summary` (one row per
#'   model, group and statistic: observed value, median and 2.5-97.5% range
#'   of the replicates, `tail_prob`), `overall` (the same for all operations
#'   together), `histogram` (observed and simulated shares per bin) and
#'   `plots` (`zeros`, `tails`, `distribution`).
#' @references Gelman, A., Meng, X.-L. and Stern, H. (1996) Posterior
#'   predictive assessment of model fitness via realized discrepancies.
#'   Statistica Sinica 6: 733-807.
#'
#'   Kleiber, C. and Zeileis, A. (2016) Visualizing count data regressions
#'   using rootograms. The American Statistician 70: 296-303.
#' @seealso [plot_residuals()] for randomised quantile residuals, which check
#'   each observation against its own predicted distribution.
#' @export
#' @examples
#' \dontrun{
#' pp <- ppc_cpue(fit, models = c("Base", "DPC+S"))
#' pp$plots$zeros
#' pp$plots$tails
#' pp$plots$distribution
#' pp$overall
#' }
ppc_cpue <- function(fit, models = NULL, n_sim = 50, by = "year",
                     param_uncertainty = TRUE, probs = c(0.9, 0.99), bins = 40,
                     seed = 1) {
  if (!inherits(fit, "cpue_fit")) stop("`fit` must come from fit_cpue().", call. = FALSE)
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the `mgcv` package.", call. = FALSE)
  }
  if (is.null(models)) models <- names(fit$models)
  for (m in models) check_model(fit, m)
  if (n_sim < 2) stop("`n_sim` must be at least 2.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  d <- fit$data
  eff <- if (is.null(fit$roles$effort)) rep(1, nrow(d)) else d[[fit$roles$effort]]
  obs <- if (isTRUE(fit$settings$offset)) d$y / eff else d$y
  grp <- if (is.null(by)) {
    factor(rep("All", nrow(d)))
  } else if (identical(by, "year")) {
    d$fyear
  } else {
    if (!by %in% names(d)) stop("Column `", by, "` not found.", call. = FALSE)
    factor(d[[by]])
  }
  all_grp <- factor(rep("All", nrow(d)))
  stat_names <- c("zeros", "mean", paste0("q", round(100 * probs)))

  pos_obs <- log10(obs[obs > 0])
  sims <- lapply(models, function(m) ppc_simulate(fit, m, n_sim, param_uncertainty, eff))
  names(sims) <- models
  pos_all <- unlist(lapply(sims, function(s) log10(s[s > 0])))
  rng <- range(c(pos_obs, stats::quantile(pos_all, c(0.001, 0.999), na.rm = TRUE)))
  breaks <- seq(rng[1], rng[2], length.out = bins + 1)

  summary_rows <- list()
  overall_rows <- list()
  hist_rows <- list()
  for (m in models) {
    s <- sims[[m]]
    summary_rows[[m]] <- ppc_compare(obs, s, grp, probs, stat_names, m)
    overall_rows[[m]] <- ppc_compare(obs, s, all_grp, probs, stat_names, m)
    h_obs <- ppc_hist(obs, breaks)
    h_sim <- apply(s, 2, ppc_hist, breaks = breaks)
    hist_rows[[m]] <- data.frame(
      model = m, mid = (breaks[-1] + breaks[-length(breaks)]) / 2,
      observed = h_obs, median = apply(h_sim, 1, stats::median),
      lower = apply(h_sim, 1, stats::quantile, 0.025),
      upper = apply(h_sim, 1, stats::quantile, 0.975))
  }
  summary <- do.call(rbind, summary_rows)
  overall <- do.call(rbind, overall_rows)
  histogram <- do.call(rbind, hist_rows)
  rownames(summary) <- rownames(overall) <- rownames(histogram) <- NULL
  for (nm in c("summary", "overall", "histogram")) {
    z <- get(nm)
    z$model <- factor(z$model, levels = models)
    assign(nm, z)
  }

  by_lab <- if (is.null(by)) "" else if (identical(by, "year")) "Year" else by
  zero_lab <- overall[overall$stat == "zeros", ]
  facet_lab <- sprintf("%s\nzeros: observed %.0f%%, simulated %.0f%% (%.0f-%.0f%%)",
                       as.character(zero_lab$model), 100 * zero_lab$observed,
                       100 * zero_lab$median,
                       100 * zero_lab$lower, 100 * zero_lab$upper)
  names(facet_lab) <- as.character(zero_lab$model)
  rot <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))

  zs <- summary[summary$stat == "zeros", ]
  zeros <- ggplot(zs, aes(x = .data$group)) +
    geom_linerange(aes(ymin = .data$lower, ymax = .data$upper), colour = "grey55",
                   linewidth = 2.5) +
    geom_point(aes(y = .data$median), colour = "grey30", shape = 95, size = 4) +
    geom_point(aes(y = .data$observed), colour = "#D55E00", size = 1.8) +
    facet_wrap(~model) +
    labs(x = by_lab, y = "Share of zero catches",
         title = "Zeros: observed (orange) and simulated (grey, 95% range)") +
    theme_bw() + rot

  ts <- summary[summary$stat %in% paste0("q", round(100 * probs)), ]
  ts$stat <- factor(ts$stat, levels = paste0("q", round(100 * probs)),
                    labels = paste0(round(100 * probs), "% quantile of positive CPUE"))
  tails <- ggplot(ts, aes(x = .data$group)) +
    geom_linerange(aes(ymin = .data$lower, ymax = .data$upper), colour = "grey55",
                   linewidth = 2.5) +
    geom_point(aes(y = .data$median), colour = "grey30", shape = 95, size = 4) +
    geom_point(aes(y = .data$observed), colour = "#D55E00", size = 1.8) +
    facet_grid(stat ~ model, scales = "free_y") +
    scale_y_log10() +
    labs(x = by_lab, y = "CPUE (log scale)",
         title = "Upper tail: observed (orange) and simulated (grey, 95% range)") +
    theme_bw() + rot

  distribution <- ggplot(histogram, aes(x = .data$mid)) +
    geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper), fill = "grey75") +
    geom_line(aes(y = .data$median), colour = "grey30", linetype = 2) +
    geom_step(aes(y = .data$observed), colour = "#D55E00", direction = "mid") +
    facet_wrap(~model, labeller = labeller(model = facet_lab)) +
    labs(x = "log10 CPUE of positive catches", y = "Share of all operations",
         title = "Positive catches: observed (orange) and simulated (grey, 95% range)") +
    theme_bw()

  structure(list(summary = summary, overall = overall, histogram = histogram,
                 plots = list(zeros = zeros, tails = tails, distribution = distribution),
                 settings = list(n_sim = n_sim, by = by, param_uncertainty = param_uncertainty,
                                 probs = probs, family = fit$settings$family)),
            class = "cpue_ppc")
}

#' @export
print.cpue_ppc <- function(x, ...) {
  cat("Predictive check:", x$settings$n_sim, "replicates,",
      if (x$settings$param_uncertainty) "with" else "without",
      "parameter uncertainty,", x$settings$family, "\n\n")
  o <- x$overall
  o$group <- NULL
  num <- c("observed", "median", "lower", "upper", "tail_prob")
  o[num] <- lapply(o[num], signif, 3)
  print(o, row.names = FALSE)
  invisible(x)
}

#' Simulate replicate CPUE data sets from one fitted model
#' @return Matrix: operations x replicates.
#' @noRd
ppc_simulate <- function(fit, model, n_sim, param_uncertainty, eff) {
  fits <- fit$models[[model]]$fits
  d <- fit$data
  n <- nrow(d)
  offset <- isTRUE(fit$settings$offset)
  if (fit$settings$family == "tweedie") {
    g <- fits$main
    mu <- exp(ppc_eta(g, d, n_sim, param_uncertainty))
    p <- as.numeric(g$family$getTheta(TRUE))
    y <- matrix(mgcv::rTweedie(as.vector(mu), p = p, phi = g$sig2), n, n_sim)
  } else {
    pr <- stats::plogis(ppc_eta(fits$binomial, d, n_sim, param_uncertainty))
    present <- matrix(stats::rbinom(n * n_sim, 1, as.vector(pr)), n, n_sim)
    gp <- fits$positive
    eta <- ppc_eta(gp, d, n_sim, param_uncertainty)
    if (fit$settings$family == "delta_lognormal") {
      val <- exp(eta + stats::rnorm(n * n_sim, 0, sqrt(gp$sig2)))
    } else {
      mu <- exp(eta)
      val <- stats::rgamma(n * n_sim, shape = 1 / gp$sig2, scale = as.vector(mu) * gp$sig2)
    }
    y <- present * matrix(val, n, n_sim)
  }
  if (offset) y <- y / eff
  y
}

#' Linear predictor for all operations under n_sim coefficient draws
#'
#' Random vessel effects of vessels not in the model (possible for the
#' positive part of delta models) are set to zero.
#' @noRd
ppc_eta <- function(g, d, n_sim, param_uncertainty) {
  nd <- d
  unseen <- rep(FALSE, nrow(d))
  if (!is.null(nd$vessel_re) && !is.null(g$model$vessel_re)) {
    known <- levels(droplevels(g$model$vessel_re))
    unseen <- !as.character(nd$vessel_re) %in% known
    if (any(unseen)) {
      nd$vessel_re <- as.character(nd$vessel_re)
      nd$vessel_re[unseen] <- known[1]
      nd$vessel_re <- factor(nd$vessel_re, levels = levels(g$model$vessel_re))
    }
  }
  X <- stats::predict(g, newdata = nd, type = "lpmatrix")
  off <- attr(X, "model.offset")
  if (is.null(off)) off <- 0
  re <- re_columns(g)
  if (length(re) && any(unseen)) X[unseen, re] <- 0
  beta <- stats::coef(g)
  B <- if (param_uncertainty) {
    matrix(mgcv::rmvn(n_sim, beta, g$Vp), nrow = n_sim)
  } else {
    matrix(beta, n_sim, length(beta), byrow = TRUE)
  }
  eta <- X %*% t(B)
  eta + off
}

#' Summary statistics of CPUE by group
#' @noRd
ppc_stats <- function(v, grp, probs) {
  sp <- split(v, grp)
  t(vapply(sp, function(z) {
    pos <- z[z > 0]
    q <- if (length(pos)) stats::quantile(pos, probs, names = FALSE) else rep(NA_real_, length(probs))
    c(mean(z == 0), mean(z), q)
  }, numeric(2 + length(probs))))
}

#' Observed against simulated summaries
#' @noRd
ppc_compare <- function(obs, sims, grp, probs, stat_names, model) {
  so <- ppc_stats(obs, grp, probs)
  ss <- vapply(seq_len(ncol(sims)), function(j) ppc_stats(sims[, j], grp, probs),
               so)
  out <- list()
  for (k in seq_along(stat_names)) {
    sk <- matrix(ss[, k, ], nrow = nrow(so))
    out[[k]] <- data.frame(
      model = model, group = factor(rownames(so), levels = rownames(so)),
      stat = stat_names[k], observed = so[, k],
      median = apply(sk, 1, stats::median, na.rm = TRUE),
      lower = apply(sk, 1, stats::quantile, 0.025, na.rm = TRUE),
      upper = apply(sk, 1, stats::quantile, 0.975, na.rm = TRUE),
      tail_prob = rowMeans(sk >= so[, k], na.rm = TRUE))
  }
  do.call(rbind, out)
}

#' Share of all operations in each bin of log10 positive CPUE
#' @noRd
ppc_hist <- function(v, breaks) {
  lv <- log10(v[v > 0])
  lv <- lv[lv >= breaks[1] & lv <= breaks[length(breaks)]]
  as.numeric(table(cut(lv, breaks, include.lowest = TRUE))) / length(v)
}
