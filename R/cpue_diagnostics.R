#' Compare fitted CPUE models
#'
#' One row per model (and per component for delta models) with the number of
#' observations, AIC, BIC, deviance explained, effective degrees of freedom,
#' the estimated Tweedie power `p` and the SD of the random vessel effect.
#'
#' Information criteria are only comparable between models fitted to the
#' same data with the same family. They tend to favour many targeting terms,
#' because the catch composition is not independent of the response (Winker
#' et al. 2014); use them together with the influence plots and the change
#' in the index.
#'
#' @param fit Output of [fit_cpue()].
#' @return A data frame.
#' @export
cpue_diagnostics <- function(fit) {
  check_model(fit, NULL)
  rows <- list()
  for (m in names(fit$models)) {
    for (comp in names(fit$models[[m]]$fits)) {
      g <- fit$models[[m]]$fits[[comp]]
      s <- summary(g)
      p <- tryCatch(as.numeric(g$family$getTheta(TRUE)), error = function(e) NA_real_)
      if (!grepl("Tweedie", g$family$family)) p <- NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        model = m, component = comp, n = s$n,
        AIC = stats::AIC(g), BIC = stats::BIC(g),
        dev_explained = round(s$dev.expl, 4), edf = round(sum(g$edf), 1),
        tweedie_p = round(p, 3), vessel_sd = round(vessel_sd(g), 3))
    }
  }
  out <- do.call(rbind, rows)
  out$dAIC <- stats::ave(out$AIC, out$component, FUN = function(a) a - min(a))
  out$AIC <- round(out$AIC, 1)
  out$BIC <- round(out$BIC, 1)
  out$dAIC <- round(out$dAIC, 1)
  out
}

#' SD of the random vessel effect
#' @noRd
vessel_sd <- function(g) {
  if (!length(re_columns(g))) return(NA_real_)
  vc <- tryCatch({
    utils::capture.output(v <- mgcv::gam.vcomp(g))
    v
  }, error = function(e) NULL)
  if (is.null(vc)) return(NA_real_)
  if (is.list(vc) && !is.null(vc$vc)) vc <- vc$vc
  if (is.matrix(vc) && "s(vessel_re)" %in% rownames(vc)) return(vc["s(vessel_re)", 1])
  if (!is.matrix(vc) && "s(vessel_re)" %in% names(vc)) return(unname(vc["s(vessel_re)"]))
  NA_real_
}

#' Plot standardised indices
#'
#' @param fit Output of [fit_cpue()].
#' @param models Models to show (default: all).
#' @param ci Logical; show 95% confidence intervals.
#' @param nominal Logical; show the nominal CPUE (black dashed line).
#' @param true Optional data frame with columns `year` and `true` (e.g. the
#'   true biomass of a simulation), scaled to mean 1 and drawn in black.
#' @return A ggplot object (relative indices, each scaled to its mean).
#' @export
plot_index <- function(fit, models = NULL, ci = TRUE, nominal = TRUE, true = NULL) {
  check_model(fit, NULL)
  idx <- fit$indices
  if (!is.null(models)) {
    idx <- idx[idx$model %in% models, ]
    idx$model <- factor(as.character(idx$model), levels = models)
  }
  p <- ggplot(idx, aes(x = .data$year, y = .data$relative, colour = .data$model))
  if (ci) {
    p <- p + geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper, fill = .data$model),
                         alpha = 0.15, colour = NA, show.legend = FALSE)
  }
  p <- p + geom_line(linewidth = 0.8) + geom_point(size = 1.2)
  if (nominal) {
    p <- p + geom_line(data = fit$nominal, aes(x = .data$year, y = .data$relative),
                       inherit.aes = FALSE, linetype = 2, colour = "grey30") +
      geom_point(data = fit$nominal, aes(x = .data$year, y = .data$relative),
                 inherit.aes = FALSE, shape = 21, fill = "white", colour = "grey30")
  }
  if (!is.null(true)) {
    true$true <- true$true / mean(true$true)
    p <- p + geom_line(data = true, aes(x = .data$year, y = .data$true),
                       inherit.aes = FALSE, colour = "black", linewidth = 1.1)
  }
  p + geom_hline(yintercept = 1, linetype = 3) +
    labs(x = "Year", y = "Relative CPUE", colour = "Model",
         title = fit$label,
         subtitle = paste(c(if (nominal) "dashed: nominal", if (!is.null(true)) "black: true"),
                          collapse = "; ")) +
    theme_bw()
}

#' Step plot of the standardisation
#'
#' Bentley et al. (2012): each panel shows the index of a model (in the
#' order fitted) together with the index of the previous model (grey,
#' dashed), so that the effect of adding each term on the index is visible.
#' Fit the models in the order in which terms are added, e.g.
#' `models = c("Base", "Covariates", "DPC")`.
#'
#' @param fit Output of [fit_cpue()].
#' @param models Models in order (default: all, as fitted).
#' @return A ggplot object.
#' @export
plot_step <- function(fit, models = NULL) {
  check_model(fit, NULL)
  if (is.null(models)) models <- names(fit$models)
  idx <- fit$indices[fit$indices$model %in% models, ]
  first <- data.frame(year = fit$nominal$year, relative = fit$nominal$relative)
  prev <- do.call(rbind, lapply(seq_along(models), function(i) {
    p <- if (i == 1) first else idx[idx$model == models[i - 1], c("year", "relative")]
    data.frame(year = p$year, relative = p$relative, model = models[i])
  }))
  idx$model <- factor(as.character(idx$model), levels = models)
  prev$model <- factor(prev$model, levels = models)
  ggplot(idx, aes(x = .data$year, y = .data$relative)) +
    geom_line(data = prev, colour = "grey55", linetype = 2) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    geom_hline(yintercept = 1, linetype = 3) +
    facet_wrap(~model, ncol = 1) +
    labs(x = "Year", y = "Relative CPUE",
         subtitle = "grey dashed: previous step (nominal CPUE for the first)") +
    theme_bw()
}

#' Influence of each model term on the index
#'
#' Bentley et al. (2012): for term \eqn{t} with effect \eqn{f_t(x_i)} on
#' operation \eqn{i}, the influence in year \eqn{y} is
#' \eqn{\exp(\overline{f_t}_y - \overline{f_t})}, the ratio between the
#' expected CPUE given the covariate values of that year and given those of
#' all years. A term has influence only if its covariate changes between
#' years: values above 1 mean the term raised the nominal CPUE of that year
#' (and the standardisation lowers it).
#'
#' @param fit Output of [fit_cpue()].
#' @param model Model name (default: the last model).
#' @param component For delta models, `"positive"` (default) or `"binomial"`.
#' @return A data frame with `term`, `year` and `influence`.
#' @export
influence_cpue <- function(fit, model = NULL, component = NULL) {
  g <- get_gam(fit, model, component)
  tm <- stats::predict(g, type = "terms")
  yr <- as.numeric(as.character(g$model$fyear))
  keep <- colnames(tm) != "fyear"
  if (!any(keep)) stop("The model has no terms besides year.", call. = FALSE)
  out <- do.call(rbind, lapply(colnames(tm)[keep], function(t) {
    f <- tm[, t]
    by_year <- tapply(f, yr, mean)
    data.frame(term = t, year = as.numeric(names(by_year)),
               influence = exp(as.numeric(by_year) - mean(f)))
  }))
  out$term <- factor(out$term, levels = colnames(tm)[keep])
  out
}

#' Plot the influence of each model term
#'
#' @inheritParams influence_cpue
#' @return A ggplot object.
#' @export
plot_influence <- function(fit, model = NULL, component = NULL) {
  model <- check_model(fit, model)
  inf <- influence_cpue(fit, model, component)
  ggplot(inf, aes(x = .data$year, y = .data$influence, colour = .data$term)) +
    geom_hline(yintercept = 1) +
    geom_line(linewidth = 0.9) + geom_point(size = 1) +
    labs(x = "Year", y = "Influence", colour = "Term",
         title = paste("Influence:", model)) +
    theme_bw()
}

#' Residual diagnostics
#'
#' For the Tweedie family, randomised quantile residuals (Dunn and Smyth
#' 1996) are computed by simulation from the fitted distribution, so that they
#' are standard normal if the model is correct, also with many zeros. For
#' delta models, deviance residuals of the chosen component are shown.
#'
#' @inheritParams influence_cpue
#' @param by Optional column of the data to plot residuals against (e.g.
#'   `"FT"`, `"month"`, `"vessel"`).
#' @param nsim Number of simulations per observation for quantile residuals.
#' @param max_points Maximum number of points in scatter plots.
#' @return A list of ggplots: `qq`, `fitted` (residuals against log fitted
#'   values with a smoother), `year` (by year) and, if `by` is given, `by`.
#' @references Dunn, P.K. and Smyth, G.K. (1996) Randomized quantile
#'   residuals. Journal of Computational and Graphical Statistics 5: 236-244.
#' @export
plot_residuals <- function(fit, model = NULL, component = NULL, by = NULL, nsim = 100,
                           max_points = 5000) {
  model <- check_model(fit, model)
  g <- get_gam(fit, model, component)
  res <- cpue_residuals(g, nsim = nsim)
  type <- attr(res, "type")
  dd <- data.frame(res = res, fitted = log(stats::fitted(g)),
                   year = factor(g$model$fyear))
  if (!is.null(by)) {
    rows <- if (identical(component, "binomial") || fit$settings$family == "tweedie") {
      seq_len(nrow(fit$data))
    } else which(fit$data$y > 0)
    if (!by %in% names(fit$data)) stop("Column `", by, "` not found.", call. = FALSE)
    dd$by <- factor(fit$data[[by]][rows])
  }
  sc <- if (nrow(dd) > max_points) dd[sample.int(nrow(dd), max_points), ] else dd
  ylab <- paste(type, "residuals")
  out <- list(
    qq = ggplot(sc, aes(sample = .data$res)) + stat_qq(size = 0.6, alpha = 0.5) +
      geom_abline(intercept = 0, slope = if (type == "quantile") 1 else stats::sd(dd$res),
                  colour = "red") +
      labs(x = "Theoretical quantiles", y = ylab, title = paste("Residuals:", model)) +
      theme_bw(),
    fitted = ggplot(sc, aes(x = .data$fitted, y = .data$res)) +
      geom_point(size = 0.6, alpha = 0.4) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "red") +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(x = "log fitted value", y = ylab) + theme_bw(),
    year = ggplot(dd, aes(x = .data$year, y = .data$res)) +
      geom_boxplot(outlier.size = 0.3, fill = "grey85") + geom_hline(yintercept = 0, linetype = 2) +
      labs(x = "Year", y = ylab) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7)))
  if (!is.null(by)) {
    out$by <- ggplot(dd, aes(x = .data$by, y = .data$res)) +
      geom_boxplot(outlier.size = 0.3, fill = "grey85") + geom_hline(yintercept = 0, linetype = 2) +
      labs(x = by, y = ylab) + theme_bw()
  }
  out
}

#' Randomised quantile residuals (Tweedie) or deviance residuals
#' @noRd
cpue_residuals <- function(g, nsim = 100) {
  if (!grepl("Tweedie", g$family$family)) {
    r <- stats::residuals(g, type = "deviance")
    attr(r, "type") <- "deviance"
    return(r)
  }
  y <- g$y
  mu <- stats::fitted(g)
  p <- as.numeric(g$family$getTheta(TRUE))
  phi <- g$sig2
  n <- length(y)
  lo <- hi <- numeric(n)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / 10000))
  for (ch in chunks) {
    sims <- matrix(mgcv::rTweedie(rep(mu[ch], nsim), p = p, phi = phi), nrow = length(ch))
    lo[ch] <- rowMeans(sims < y[ch])
    hi[ch] <- rowMeans(sims <= y[ch])
  }
  u <- stats::runif(n, lo, hi)
  eps <- 1 / (2 * nsim)
  r <- stats::qnorm(pmin(pmax(u, eps), 1 - eps))
  attr(r, "type") <- "quantile"
  r
}

#' Partial effects of the covariates
#'
#' Effect of each covariate on the expected CPUE (relative to its mean
#' effect), with the other covariates at their reference values and the
#' random vessel effect at zero.
#'
#' @inheritParams influence_cpue
#' @return A list of ggplots: `numeric` (smooth effects, e.g. PCs and month,
#'   with 95% CIs), `factor` (factor effects, e.g. fishing tactic), and
#'   `space` (map of the spatial effect over the observed cells), each only
#'   if present.
#' @export
plot_effects <- function(fit, model = NULL, component = NULL) {
  model <- check_model(fit, model)
  g <- get_gam(fit, model, component)
  d <- fit$data
  vars <- setdiff(all.vars(stats::as.formula(paste("~", fit$models[[model]]$rhs))),
                  c("fyear", "vessel_re", "log_effort"))
  base <- d[1, , drop = FALSE]
  for (v in names(d)) {
    z <- d[[v]]
    base[[v]] <- if (is.numeric(z)) stats::median(z) else mode_value(z)
  }
  base$log_effort <- 0
  excl <- if (length(re_columns(g))) "s(vessel_re)" else NULL
  pred <- function(nd) {
    pr <- stats::predict(g, newdata = nd, se.fit = TRUE, exclude = excl)
    data.frame(fit = as.numeric(pr$fit), se = as.numeric(pr$se.fit))
  }
  num <- list()
  fac <- list()
  for (v in setdiff(vars, c("pos_x", "pos_y"))) {
    z <- d[[v]]
    if (is.numeric(z) && length(unique(z)) > 12) {
      grid <- seq(stats::quantile(z, 0.01), stats::quantile(z, 0.99), length.out = 50)
    } else {
      grid <- sort(unique(z))
    }
    nd <- base[rep(1, length(grid)), , drop = FALSE]
    nd[[v]] <- grid
    pr <- pred(nd)
    pr$fit <- pr$fit - mean(pr$fit)
    eff <- data.frame(variable = v, value = grid, effect = exp(pr$fit),
                      lower = exp(pr$fit - 1.96 * pr$se), upper = exp(pr$fit + 1.96 * pr$se))
    if (is.numeric(z) && length(unique(z)) > 12) num[[v]] <- eff else {
      eff$value <- as.character(eff$value)
      fac[[v]] <- eff
    }
  }
  out <- list()
  if (length(num)) {
    out$numeric <- ggplot(do.call(rbind, num), aes(x = .data$value, y = .data$effect)) +
      geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper), fill = "grey80") +
      geom_line() + geom_hline(yintercept = 1, linetype = 2) +
      geom_rug(data = data.frame(variable = rep(names(num), each = min(nrow(d), 2000)),
                                 value = unlist(lapply(names(num), function(v) {
                                   utils::head(d[[v]], 2000)
                                 }))),
               aes(x = .data$value), inherit.aes = FALSE, alpha = 0.1, sides = "b") +
      facet_wrap(~variable, scales = "free_x") +
      labs(x = NULL, y = "Relative effect on CPUE", title = paste("Effects:", model)) +
      theme_bw()
  }
  if (length(fac)) {
    out$factor <- ggplot(do.call(rbind, fac), aes(x = .data$value, y = .data$effect)) +
      geom_hline(yintercept = 1, linetype = 2) +
      geom_errorbar(aes(ymin = .data$lower, ymax = .data$upper), width = 0.2) +
      geom_point(size = 2) + facet_wrap(~variable, scales = "free_x") +
      labs(x = NULL, y = "Relative effect on CPUE", title = paste("Effects:", model)) +
      theme_bw()
  }
  if (all(c("pos_x", "pos_y") %in% vars)) {
    cells <- unique(d[, c("pos_x", "pos_y")])
    nd <- base[rep(1, nrow(cells)), , drop = FALSE]
    nd$pos_x <- cells$pos_x
    nd$pos_y <- cells$pos_y
    pr <- pred(nd)
    cells$effect <- exp(pr$fit - mean(pr$fit))
    out$space <- ggplot(cells, aes(x = .data$pos_x, y = .data$pos_y, fill = log(.data$effect))) +
      geom_tile() +
      scale_fill_gradient2(low = "#2c7fb8", mid = "white", high = "#d95f0e") +
      coord_quickmap() +
      labs(x = fit$roles$lon, y = fit$roles$lat, fill = "log relative\neffect",
           title = paste("Spatial effect:", model)) + theme_bw()
  }
  out
}

#' Plot the random vessel effects
#'
#' @inheritParams influence_cpue
#' @return A ggplot of the estimated vessel effects (log scale) ordered by
#'   size, with bands of one and 1.96 estimated SDs.
#' @export
plot_vessel_effects <- function(fit, model = NULL, component = NULL) {
  model <- check_model(fit, model)
  g <- get_gam(fit, model, component)
  cols <- re_columns(g)
  if (!length(cols)) stop("The model has no random vessel effect (use vessel = TRUE).",
                          call. = FALSE)
  re <- data.frame(vessel = levels(fit$data$vessel_re), effect = stats::coef(g)[cols])
  re <- re[order(re$effect), ]
  re$rank <- seq_len(nrow(re))
  s <- vessel_sd(g)
  if (is.na(s)) s <- stats::sd(re$effect)
  ggplot(re, aes(x = .data$rank, y = .data$effect)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.96 * s, ymax = 1.96 * s,
             fill = "grey90") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -s, ymax = s, fill = "grey75") +
    geom_hline(yintercept = 0) + geom_point() +
    scale_x_continuous(breaks = re$rank, labels = re$vessel) +
    labs(x = "Vessel", y = "Random effect (log scale)",
         title = paste("Vessel effects:", model),
         subtitle = sprintf("estimated SD = %.2f", s)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
}

#' The GAM of a model (and component)
#' @noRd
get_gam <- function(fit, model = NULL, component = NULL) {
  model <- check_model(fit, model)
  fits <- fit$models[[model]]$fits
  if (is.null(component)) component <- if ("main" %in% names(fits)) "main" else "positive"
  if (!component %in% names(fits)) {
    stop("Component `", component, "` not found; available: ",
         paste(names(fits), collapse = ", "), call. = FALSE)
  }
  fits[[component]]
}
