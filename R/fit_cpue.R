#' Standardise CPUE for one species with mgcv
#'
#' Fits a set of GAMs or GAMMs (`mgcv`) to the CPUE of one species and
#' derives a standardised annual index from each. This is the model step of
#' the real-data workflow, after [cpue_data()], [select_species()] and
#' [add_targeting()].
#'
#' **Models.** Every model contains the year effect (`fyear`, a factor) and
#' differs in the targeting term:
#' * `Base`: year only;
#' * `Covariates`: year + `terms` (e.g. month, port, gear);
#' * `DPC`: year + `terms` + `s(PC1) + ... + s(PCn)` (direct principal
#'   components, Winker et al. 2013, 2014);
#' * `Cluster`: year + `terms` + `FT` (fishing tactic as a factor).
#'
#' Append `+S` to a name to add the spatial smoother `s(pos_x, pos_y)` of the
#' positions given in [cpue_data()] (e.g. `"DPC+S"`). Further models can be
#' given as right-hand sides in `add_models`. With `vessel = TRUE`, every model
#' gets a random vessel intercept `s(vessel_re, bs = "re")`.
#'
#' **Families.**
#' * `"tweedie"` (default): `mgcv::tw()`, a compound Poisson-gamma
#'   distribution with a point mass at zero. The power parameter `p` is
#'   estimated jointly with the model, also with random vessel effects, so
#'   the two-step approach (estimate `p` in a GAM, then fix it in `gamm()`)
#'   is not needed.
#' * `"delta_lognormal"`, `"delta_gamma"`: a binomial GAM for the probability
#'   of a positive catch and a lognormal (Gaussian on log CPUE) or gamma GAM
#'   for positive CPUE, with the same right-hand side. The index is the
#'   product of both.
#'
#' **Index.** For each year, the expected CPUE is predicted for a reference
#' set of covariate values and averaged: positions are averaged over all
#' observed cells, the variables in `average_over` over their observed
#' values (by default fishing tactic and month), other numeric covariates
#' (e.g. PC scores) are set to their median and other factors to their most
#' common level. Random vessel effects are set to zero (the average vessel).
#' Standard errors of the log index come from the delta method using the
#' Bayesian covariance matrix of the coefficients (`Vp`).
#'
#' @param x A `cpue_data` object (after [add_targeting()] for the `DPC` and
#'   `Cluster` models).
#' @param species The species to standardise (column or original name).
#' @param models Model names (see Details). Default: `Base`, `Covariates`
#'   (if `terms` is given), and `DPC` and `Cluster` (if targeting covariates
#'   are available).
#' @param terms Character vector of additional model terms used in all
#'   models except `Base`, e.g. `c('s(month, bs = "cc", k = 8)', "port")`.
#'   They can use any column of `x$data`.
#' @param vessel Logical; add a random vessel effect to every model.
#' @param family `"tweedie"`, `"delta_lognormal"` or `"delta_gamma"`.
#' @param n_pc Number of PCs in the `DPC` models (default: the number used by
#'   [add_targeting()]).
#' @param pc_k Basis dimension of each PC smoother.
#' @param spatial_k Basis dimension of the spatial smoother.
#' @param offset Logical. `FALSE` (default): the response is CPUE (catch /
#'   effort). `TRUE`: the response is catch with `offset(log(effort))`
#'   (Tweedie and the positive part of delta models).
#' @param average_over Variables whose observed values are averaged over to
#'   compute the index (see Details).
#' @param max_ref Maximum number of reference combinations per year (a
#'   random sample is used if there are more).
#' @param add_models Optional named character vector of further right-hand
#'   sides; they must contain `fyear`, e.g.
#'   `c(Interaction = "fyear + s(PC1, by = FT) + FT")`.
#' @param engine `"gam"` (default) or `"bam"` (with `discrete = TRUE`, much
#'   faster for large data).
#' @param seed Seed for sampling reference combinations.
#' @param ... Further arguments for `mgcv::gam()` or `mgcv::bam()`, e.g.
#'   `knots = list(month = c(0.5, 12.5))` for a cyclic month smoother.
#'
#' @return An object of class `cpue_fit`: a list with `models` (for each
#'   model the right-hand side and the fitted GAM(s)), `indices` (annual
#'   index of every model: absolute and relative to its mean, log SE and 95%
#'   CI), `nominal` (nominal CPUE by year), `data` (the data used, with the
#'   model columns), `reference` (the reference set used for the index),
#'   `species`, `label` and `settings`. See [cpue_index()],
#'   [cpue_diagnostics()], [plot_index()], [plot_influence()].
#' @references Bentley, N., Kendrick, T.H., Starr, P.J. and Breen, P.A.
#'   (2012) Influence plots and metrics: tools for better understanding
#'   fisheries catch-per-unit-effort standardizations. ICES Journal of
#'   Marine Science 69: 84-88.
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_cpue(x, species = "Pagellus bogaraveo",
#'                 terms = 's(month, bs = "cc", k = 8)', vessel = TRUE,
#'                 knots = list(month = c(0.5, 12.5)))
#' cpue_diagnostics(fit)
#' plot_index(fit)
#' plot_influence(fit, "DPC")
#' cpue_index(fit, "DPC")
#' }
fit_cpue <- function(x, species, models = NULL, terms = NULL, vessel = FALSE,
                     family = c("tweedie", "delta_lognormal", "delta_gamma"),
                     n_pc = NULL, pc_k = 4, spatial_k = 30, offset = FALSE,
                     average_over = c("FT", "month"), max_ref = 5000,
                     add_models = NULL, engine = c("gam", "bam"), seed = 1, ...) {
  check_cpue_data(x)
  family <- match.arg(family)
  engine <- match.arg(engine)
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the `mgcv` package to fit CPUE models.", call. = FALSE)
  }
  if (length(species) != 1) stop("Give one species.", call. = FALSE)
  sp <- species_codes(x, species)
  r <- x$roles
  d <- x$data
  years <- sort(unique(d[[r$year]]))
  d$fyear <- factor(d[[r$year]], levels = years)
  if (offset && is.null(r$effort)) stop("`offset = TRUE` needs an effort column.", call. = FALSE)
  eff <- if (is.null(r$effort)) rep(1, nrow(d)) else d[[r$effort]]
  d$y <- if (offset) d[[sp]] else d[[sp]] / eff
  d$log_effort <- log(eff)

  has_tg <- !is.null(x$targeting) && !is.null(d$FT)
  if (has_tg && nrow(x$targeting$pca$scores) != nrow(d)) {
    stop("The data were subset after add_targeting(); run add_targeting() again.",
         call. = FALSE)
  }
  if (has_tg && is.null(n_pc)) n_pc <- x$targeting$n_pc
  has_pos <- !is.null(r$lon) && !is.null(r$lat)
  if (has_pos) {
    d$pos_x <- d[[r$lon]]
    d$pos_y <- d[[r$lat]]
  }
  if (vessel) {
    if (is.null(r$vessel)) stop("No vessel column was defined in cpue_data().", call. = FALSE)
    d$vessel_re <- factor(d[[r$vessel]])
  }

  # Model right-hand sides
  pieces <- list(Base = character(0))
  if (length(terms)) pieces$Covariates <- terms
  if (has_tg) {
    pieces$DPC <- c(terms, sprintf("s(PC%d, k = %d)", seq_len(n_pc), pc_k))
    pieces$Cluster <- c(terms, "FT")
  }
  if (is.null(models)) models <- names(pieces)
  models <- unique(models)
  with_space <- grepl("\\+S$", models)
  base_name <- sub("\\+S$", "", models)
  unknown <- models[!base_name %in% names(pieces)]
  if (length(unknown)) {
    hint <- if (any(sub("\\+S$", "", unknown) %in% c("DPC", "Cluster")) && !has_tg) {
      " (run add_targeting() first)"
    } else if (any(sub("\\+S$", "", unknown) == "Covariates")) " (give `terms`)" else ""
    stop("Unknown or unavailable models: ", paste(unknown, collapse = ", "), hint, call. = FALSE)
  }
  if (any(with_space) && !has_pos) {
    stop("`+S` models need positions (lon, lat) in cpue_data().", call. = FALSE)
  }
  spatial_term <- sprintf("s(pos_x, pos_y, k = %d)", spatial_k)
  rhs <- lapply(seq_along(models), function(i) {
    c("fyear", if (with_space[i]) spatial_term, pieces[[base_name[i]]])
  })
  names(rhs) <- models
  if (!is.null(add_models)) {
    if (is.null(names(add_models)) || any(names(add_models) == "")) {
      stop("`add_models` must be a named character vector.", call. = FALSE)
    }
    if (!all(grepl("\\bfyear\\b", add_models))) {
      stop("Every model in `add_models` must contain `fyear`.", call. = FALSE)
    }
    rhs <- c(rhs, as.list(add_models))
  }
  if (vessel) rhs <- lapply(rhs, c, 's(vessel_re, bs = "re")')
  rhs <- vapply(rhs, paste, character(1), collapse = " + ")
  model_vars <- setdiff(unique(unlist(lapply(rhs, function(f) {
    all.vars(stats::as.formula(paste("~", f)))
  }))), c("fyear", "vessel_re", "log_effort"))
  missing_vars <- setdiff(model_vars, names(d))
  if (length(missing_vars)) {
    stop("Model variables not found in the data: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }

  if (family != "tweedie") {
    d$positive <- as.numeric(d$y > 0)
    share <- tapply(d$positive, d$fyear, mean)
    if (any(share %in% c(0, 1))) {
      warning("Years with no or only positive catches: ",
              paste(names(share)[share %in% c(0, 1)], collapse = ", "),
              "; the binomial year effects are not estimable there.", call. = FALSE)
    }
  }

  fitter <- function(form, data, fam) {
    if (engine == "gam") {
      mgcv::gam(form, data = data, family = fam, method = "REML", ...)
    } else {
      mgcv::bam(form, data = data, family = fam, method = "fREML", discrete = TRUE, ...)
    }
  }
  off <- if (offset) " + offset(log_effort)" else ""
  reference <- reference_grid(d, model_vars, average_over, max_ref, seed)

  out_models <- list()
  indices <- list()
  for (m in names(rhs)) {
    fits <- tryCatch({
      if (family == "tweedie") {
        list(main = fitter(stats::as.formula(paste("y ~", rhs[[m]], off)), d,
                           tweedie_family()))
      } else {
        pos <- d[d$y > 0, , drop = FALSE]
        pos_resp <- if (family == "delta_lognormal") "log(y) ~" else "y ~"
        pos_fam <- if (family == "delta_lognormal") stats::gaussian() else
          stats::Gamma(link = "log")
        list(binomial = fitter(stats::as.formula(paste("positive ~", rhs[[m]])), d,
                               stats::binomial()),
             positive = fitter(stats::as.formula(paste(pos_resp, rhs[[m]], off)), pos, pos_fam))
      }
    }, error = function(e) {
      stop("Model `", m, "` failed: ", conditionMessage(e), call. = FALSE)
    })
    idx <- model_index(fits, family, reference, years)
    idx$model <- m
    indices[[m]] <- idx
    out_models[[m]] <- list(rhs = rhs[[m]], fits = fits)
  }
  indices <- do.call(rbind, indices)
  rownames(indices) <- NULL
  indices$model <- factor(indices$model, levels = names(rhs))

  nominal <- data.frame(year = years,
                        n = as.numeric(table(d$fyear)),
                        positive = as.numeric(tapply(d$y > 0, d$fyear, mean)),
                        nominal = as.numeric(tapply(d[[sp]] / eff, d$fyear, mean)))
  nominal$relative <- nominal$nominal / mean(nominal$nominal)

  structure(list(models = out_models, indices = indices, nominal = nominal, data = d,
                 reference = reference, species = sp, label = sp_label(x, sp),
                 roles = r,
                 settings = list(family = family, vessel = vessel, offset = offset,
                                 n_pc = n_pc, engine = engine, terms = terms,
                                 average_over = intersect(average_over, model_vars))),
            class = "cpue_fit")
}

#' Reference covariate values over which the annual index is averaged
#' @noRd
reference_grid <- function(d, vars, average_over, max_ref, seed) {
  parts <- list()
  if (all(c("pos_x", "pos_y") %in% vars)) {
    parts$space <- unique(d[, c("pos_x", "pos_y")])
    vars <- setdiff(vars, c("pos_x", "pos_y"))
  }
  for (v in vars) {
    z <- d[[v]]
    val <- if (v %in% average_over) {
      u <- unique(z)
      u[order(u)]
    } else if (is.numeric(z)) {
      stats::median(z)
    } else {
      mode_value(z)
    }
    df <- data.frame(val, stringsAsFactors = FALSE)
    names(df) <- v
    parts[[v]] <- df
  }
  ref <- data.frame(.one = 1)
  for (p in parts) ref <- merge(ref, p, by = NULL)
  ref$.one <- NULL
  if (nrow(ref) > max_ref) {
    if (!is.null(seed)) set.seed(seed)
    ref <- ref[sample.int(nrow(ref), max_ref), , drop = FALSE]
  }
  ref$log_effort <- 0
  if (!is.null(d$vessel_re)) ref$vessel_re <- factor(levels(d$vessel_re)[1],
                                                     levels = levels(d$vessel_re))
  rownames(ref) <- NULL
  ref
}

#' Most common value, keeping factor levels
#' @noRd
mode_value <- function(z) {
  tab <- table(z)
  top <- names(tab)[which.max(tab)]
  if (is.factor(z)) factor(top, levels = levels(z)) else {
    if (is.logical(z)) as.logical(top) else top
  }
}

#' Columns of the random vessel effect in a GAM
#' @noRd
re_columns <- function(fit) {
  labs <- vapply(fit$smooth, function(s) s$label, character(1))
  i <- which(labs == "s(vessel_re)")
  if (!length(i)) return(integer(0))
  fit$smooth[[i]]$first.para:fit$smooth[[i]]$last.para
}

#' Mean prediction over the reference set and its log SE, by year
#' @noRd
year_means <- function(fit, reference, years, type = c("log", "logit"), add = 0) {
  type <- match.arg(type)
  beta <- stats::coef(fit)
  vp <- fit$Vp
  re <- re_columns(fit)
  levels_y <- levels(fit$model$fyear)
  out <- lapply(years, function(yr) {
    nd <- reference
    nd$fyear <- factor(as.character(yr), levels = levels_y)
    X <- stats::predict(fit, newdata = nd, type = "lpmatrix")
    if (length(re)) X[, re] <- 0
    eta <- drop(X %*% beta) + add
    if (type == "log") {
      mu <- exp(eta)
      g <- colMeans(X * mu)
    } else {
      mu <- stats::plogis(eta)
      g <- colMeans(X * (mu * (1 - mu)))
    }
    m <- mean(mu)
    se <- sqrt(max(0, drop(t(g) %*% vp %*% g)))
    c(mean = m, se_log = se / m)
  })
  out <- do.call(rbind, out)
  data.frame(year = years, mean = out[, "mean"], se_log = out[, "se_log"])
}

#' Annual index of one model
#' @noRd
model_index <- function(fits, family, reference, years) {
  if (family == "tweedie") {
    ym <- year_means(fits$main, reference, years, "log")
    idx <- data.frame(year = years, index = ym$mean, se_log = ym$se_log, prob = NA_real_)
  } else {
    pb <- year_means(fits$binomial, reference, years, "logit")
    add <- if (family == "delta_lognormal") fits$positive$sig2 / 2 else 0
    pm <- year_means(fits$positive, reference, years, "log", add = add)
    idx <- data.frame(year = years, index = pb$mean * pm$mean,
                      se_log = sqrt(pb$se_log^2 + pm$se_log^2), prob = pb$mean)
  }
  idx$relative <- idx$index / mean(idx$index)
  idx$lower <- idx$relative * exp(-1.96 * idx$se_log)
  idx$upper <- idx$relative * exp(1.96 * idx$se_log)
  idx
}

#' @export
print.cpue_fit <- function(x, ...) {
  cat("<cpue_fit>", x$label, "-", x$settings$family,
      if (x$settings$vessel) "with random vessel effects" else "", "\n")
  cat("  operations:", nrow(x$data), "; years:", min(x$nominal$year), "-",
      max(x$nominal$year), "\n")
  if (length(x$settings$average_over)) {
    cat("  index averaged over:", paste(x$settings$average_over, collapse = ", "),
        if (!is.null(x$data$pos_x)) "(and positions in spatial models)" else "", "\n")
  }
  cat("Models:\n")
  for (m in names(x$models)) cat(sprintf("  %-12s y ~ %s\n", m, x$models[[m]]$rhs))
  invisible(x)
}

#' Standardised index table
#'
#' @param fit Output of [fit_cpue()].
#' @param model Model name (default: the last model).
#' @return A data frame with `year`, `n` (operations), `nominal` (mean
#'   CPUE), `index` (standardised CPUE in the units of the data, at the
#'   reference covariate values), `relative` (index / its mean), `se_log`
#'   (standard error of the log index, the usual input to stock
#'   assessments), `cv`, and the 95% CI of the relative index.
#' @export
cpue_index <- function(fit, model = NULL) {
  model <- check_model(fit, model)
  idx <- fit$indices[fit$indices$model == model, ]
  data.frame(year = idx$year, n = fit$nominal$n, nominal = fit$nominal$nominal,
             index = idx$index, relative = idx$relative, se_log = idx$se_log,
             cv = sqrt(exp(idx$se_log^2) - 1), lower = idx$lower, upper = idx$upper)
}

#' Write the standardised index to a CSV file
#'
#' @inheritParams cpue_index
#' @param file Path of the CSV file.
#' @param digits Number of significant digits.
#' @return The index table, invisibly.
#' @export
export_index <- function(fit, file, model = NULL, digits = 4) {
  tab <- cpue_index(fit, model)
  num <- vapply(tab, is.numeric, logical(1)) & names(tab) != "year" & names(tab) != "n"
  tab[num] <- lapply(tab[num], signif, digits = digits)
  utils::write.csv(tab, file, row.names = FALSE)
  invisible(tab)
}

#' Check a model name, defaulting to the last model
#' @noRd
check_model <- function(fit, model) {
  if (!inherits(fit, "cpue_fit")) stop("`fit` must come from fit_cpue().", call. = FALSE)
  if (is.null(model)) model <- utils::tail(names(fit$models), 1)
  if (!model %in% names(fit$models)) {
    stop("Model `", model, "` not found; available: ",
         paste(names(fit$models), collapse = ", "), call. = FALSE)
  }
  model
}
