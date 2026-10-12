#' Standardise simulated CPUE with alternative targeting covariates
#'
#' Fits Tweedie GAMs (`mgcv::gam` with `tw()`, as in the CPUE guidelines) to
#' each species and compares the annual indices with the true biomass trend.
#'
#' Each model is `year` plus a targeting term:
#' * `Base`: no targeting term.
#' * `DPC`: `s(PC1) + ... + s(PCn)` from the catch composition (direct
#'   principal components, Winker et al. 2013, 2014), with `n = n_pc`.
#' * `Cluster`: fishing tactic from k-means clustering of the catch
#'   composition (as a factor).
#' * `HBF`: HBF as a factor (tuna model only).
#' * `Cluster_HBF`: Cluster + HBF (tuna model only).
#' * `DPC_HBF`: DPC + HBF (tuna model only).
#' * `HBS`: `offset(log(hbs_avail))`, a specialised habitat-based
#'   standardisation (Bigelow et al. 2002): effort is weighted by the overlap
#'   between the expected hook profile of the HBF and the species' depth
#'   profile under the average thermocline (tuna model only; not fitted by
#'   default).
#'
#' Append `+S` to any model name to add a spatial smoother `s(Xs, Ys)`, e.g.
#' `"DPC+S"`; `"Base+S"` is the spatial-only model. `Xs`, `Ys` are the
#' centres of `spatial_res` x `spatial_res` blocks of cells, i.e. the
#' resolution at which positions are reported.
#'
#' The latent `hook_depth` and realised `vertical_avail` are never used as
#' predictors.
#'
#' @param simulation Output of [simulate_tuna()] or [simulate_bathymetric()].
#' @param models Models to fit, e.g. `c("Base", "Base+S", "DPC+S", "HBF+S")`.
#'   Default: `Base`, plus `Base`, `Cluster`, `DPC`, `HBF`, `Cluster_HBF` and
#'   `DPC_HBF` with `+S` (the HBF models only for the tuna model).
#' @param spatial_res Size of the reporting blocks in cells (5 = 5 x 5 blocks
#'   as in the CPUE guidelines; 1 = full resolution).
#' @param spatial_k Basis dimension of the spatial smoother. By default it
#'   grows with the number of blocks (8 for 100 blocks, at most 60).
#' @param add_models Optional named character vector of extra right-hand
#'   sides, e.g. `c(HBF_period = "year + s(Xs, Ys, k = 8) + HBF_period")`.
#'   They can use any column of the returned `data`, including `Xs`, `Ys`,
#'   `period` (10-year block as a factor) and `HBF_period` (HBF x period).
#' @param vessel_re Logical; add a random vessel effect `s(vessel, bs = "re")`.
#' @param n_pc Number of principal components used by the DPC models
#'   (default 3; the 2023 CPUE guidelines used 2; at most the number of
#'   species). More PCs capture more of the targeting signal, but also more
#'   risk of absorbing real changes in relative abundance.
#' @param n_clusters Number of fishing tactics for the `Cluster` models. By
#'   default `n_pc + 1`, the rule of Winker et al. (2014).
#' @param seed Seed for the k-means clustering.
#' @param engine `"gam"` (default) or `"bam"` with `discrete = TRUE`, which is
#'   much faster for large data sets.
#' @param label Optional text appended to every model name, e.g. `" (1x1)"`,
#'   to tell fits at different resolutions apart after [bind_comparisons()].
#'
#' @return A list with `indices` (species, year, model, index with 95% CI and
#'   the true index), `diagnostics` (root mean squared error of the log index
#'   by species and model), `formulas` (the right-hand side of each model),
#'   `data` (the set-level data used, with PCs, tactics and block centres) and
#'   `pca` (loadings and variances of the catch-composition PCA; see
#'   [plot_pcs()]).
#' @export
#'
#' @examples
#' \dontrun{
#' om <- simulate_tuna(seed = 42)
#' comparison <- standardize_tuna_indices(om, models = c("Base", "Base+S",
#'                                         "DPC+S", "Cluster+S", "HBF+S"))
#' comparison$diagnostics
#' plot_indices(comparison)
#' }
standardize_tuna_indices <- function(simulation, models = NULL,
                                     spatial_res = 5, spatial_k = NULL,
                                     add_models = NULL, vessel_re = FALSE,
                                     n_pc = 3, n_clusters = n_pc + 1, seed = 1,
                                     engine = c("gam", "bam"), label = NULL) {
  engine <- match.arg(engine)
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the `mgcv` package to fit CPUE models.", call. = FALSE)
  }
  if (is.null(simulation$Data) || is.null(simulation$Biomass)) {
    stop("`simulation` must be the output of simulate_tuna() or simulate_bathymetric().",
         call. = FALSE)
  }
  species <- levels(simulation$Biomass$Sp)
  d <- as_set_data(simulation$Data)
  has_hbf <- "HBF" %in% names(d)
  if (has_hbf) {
    hbs <- as_set_data(simulation$Data, value = "hbs_avail")
    for (sp in species) d[[paste0("hbs_", sp)]] <- hbs[[sp]]
  }

  # Catch composition: proportions, fourth-root, PCA (as in the guidelines)
  catch <- as.matrix(d[, species])
  keep <- rowSums(catch) > 0
  d <- d[keep, , drop = FALSE]
  composition <- (catch[keep, , drop = FALSE] / rowSums(catch[keep, , drop = FALSE]))^0.25
  pca <- stats::prcomp(composition, scale. = TRUE)
  if (n_pc < 1 || n_pc > ncol(pca$x)) {
    stop("`n_pc` must be between 1 and ", ncol(pca$x), call. = FALSE)
  }
  pc_names <- colnames(pca$x)          # PC1, PC2, ... (all kept in `data`)
  for (j in seq_along(pc_names)) d[[pc_names[j]]] <- pca$x[, j]
  if (!is.null(seed)) set.seed(seed)
  d$FT <- factor(stats::kmeans(composition, centers = n_clusters, nstart = 20)$cluster)

  # Reporting resolution of positions
  d$Xs <- coarse_coord(d$X, spatial_res)
  d$Ys <- coarse_coord(d$Y, spatial_res)
  if (is.null(spatial_k)) {
    n_blocks <- nrow(unique(d[, c("Xs", "Ys")]))
    spatial_k <- min(60, max(8, round(n_blocks / 12)), n_blocks - 1)
  }

  years <- sort(unique(d$year))
  d$period <- factor((d$year - 1) %/% 10 + 1)
  d$year <- factor(d$year, levels = years)
  d$vessel <- factor(d$vessel)
  if (has_hbf) {
    d$HBF <- factor(d$HBF)
    d$HBF_period <- interaction(d$HBF, d$period, drop = TRUE)
  }

  # Targeting terms
  dpc <- sprintf("s(PC%d, k = 3)", seq_len(n_pc))
  targeting <- list(Base = character(0), DPC = dpc, Cluster = "FT")
  if (has_hbf) {
    targeting$HBF <- "HBF"
    targeting$Cluster_HBF <- c("FT", "HBF")
    targeting$DPC_HBF <- c(dpc, "HBF")
    targeting$HBS <- "offset(log(hbs))"
  }
  if (is.null(models)) {
    models <- c("Base", paste0(setdiff(names(targeting), "HBS"), "+S"))
  }
  models <- unique(models)
  with_space <- grepl("\\+S$", models)
  base_name <- sub("\\+S$", "", models)
  unknown <- models[!base_name %in% names(targeting)]
  if (length(unknown)) stop("Unknown or inapplicable models: ",
                            paste(unknown, collapse = ", "), call. = FALSE)

  # Model terms: year + (space) + targeting term
  spatial_term <- sprintf("s(Xs, Ys, k = %d)", spatial_k)
  terms <- lapply(seq_along(models), function(i) {
    c("year", if (with_space[i]) spatial_term, targeting[[base_name[i]]])
  })
  names(terms) <- models
  if (!is.null(add_models)) {
    if (is.null(names(add_models)) || any(names(add_models) == "")) {
      stop("`add_models` must be a named character vector.", call. = FALSE)
    }
    terms <- c(terms, as.list(add_models))
  }
  if (vessel_re) terms <- lapply(terms, c, 's(vessel, bs = "re")')
  formulas <- vapply(terms, paste, character(1), collapse = " + ")

  # Reference levels for prediction: only `year` changes
  mode_level <- function(f) names(which.max(table(f)))
  newdata <- data.frame(year = factor(years, levels = years),
                        Xs = stats::median(d$Xs), Ys = stats::median(d$Ys),
                        Xobs = stats::median(d$Xobs), Yobs = stats::median(d$Yobs),
                        period = factor(1, levels = levels(d$period)),
                        FT = factor(mode_level(d$FT), levels = levels(d$FT)),
                        vessel = factor(levels(d$vessel)[1], levels = levels(d$vessel)),
                        hbs = 1)
  for (pc in pc_names) newdata[[pc]] <- 0
  if (has_hbf) {
    newdata$HBF <- factor(mode_level(d$HBF), levels = levels(d$HBF))
    newdata$HBF_period <- factor(mode_level(d$HBF_period[d$period == "1"]),
                                 levels = levels(d$HBF_period))
  }

  results <- list()
  for (sp in species) {
    d$y <- d[[sp]]
    if (has_hbf) d$hbs <- pmax(d[[paste0("hbs_", sp)]], 1e-6)
    truth <- simulation$Biomass$True[simulation$Biomass$Sp == sp]
    for (m in names(formulas)) {
      form <- stats::as.formula(paste("y ~", formulas[[m]]))
      fit <- if (engine == "gam") {
        mgcv::gam(form, data = d, family = tweedie_family(), method = "REML")
      } else {
        mgcv::bam(form, data = d, family = tweedie_family(), method = "fREML",
                  discrete = TRUE)
      }
      pred <- stats::predict(fit, newdata = newdata, type = "link", se.fit = TRUE,
                             exclude = if (vessel_re) "s(vessel)" else NULL)
      fit_link <- as.numeric(pred$fit)
      se <- as.numeric(pred$se.fit)
      scale <- mean(exp(fit_link))
      results[[length(results) + 1]] <- data.frame(
        species = sp, year = years, model = m,
        index = exp(fit_link) / scale,
        lower = exp(fit_link - 1.96 * se) / scale,
        upper = exp(fit_link + 1.96 * se) / scale,
        true = truth[years])
    }
  }
  indices <- do.call(rbind, results)
  if (!is.null(label)) {
    indices$model <- paste0(indices$model, label)
    names(formulas) <- paste0(names(formulas), label)
  }
  indices$species <- factor(indices$species, levels = species)
  indices$model <- factor(indices$model, levels = names(formulas))
  err <- (log(indices$index) - log(indices$true))^2
  diagnostics <- stats::aggregate(list(rmse_log = err),
                                  by = list(species = indices$species, model = indices$model),
                                  FUN = function(e) sqrt(mean(e)))
  list(indices = indices, diagnostics = diagnostics, formulas = formulas, data = d,
       pca = list(loadings = pca$rotation, variance = pca$sdev^2))
}

#' mgcv's tw() family that also works when mgcv is not attached
#'
#' The functions inside `mgcv::tw()` are defined in an environment whose
#' parent is the global environment, so they cannot find mgcv's internal
#' functions (e.g. `ldTweedie()`) unless mgcv is attached with `library()`.
#' Re-parenting that environment to the mgcv namespace fixes this.
#' @noRd
tweedie_family <- function() {
  fam <- mgcv::tw()
  envs <- unique(lapply(Filter(is.function, fam), environment))
  for (env in envs) {
    if (is.environment(env) && identical(parent.env(env), globalenv())) {
      parent.env(env) <- asNamespace("mgcv")
    }
  }
  fam
}

#' Combine several standardisation results into one
#'
#' Useful to compare fits made with different settings (e.g. spatial
#' resolutions) in one plot. Give each call a different `label` so that the
#' model names are unique.
#'
#' @param ... Outputs of [standardize_tuna_indices()].
#' @return An object like the output of [standardize_tuna_indices()] with the
#'   `indices`, `diagnostics` and `formulas` of all inputs (`data` and `pca`
#'   are taken from the first).
#' @export
#' @examples
#' \dontrun{
#' a <- standardize_tuna_indices(om, models = c("Base", "HBF+S"))
#' b <- standardize_tuna_indices(om, models = "Base+S", spatial_res = 1,
#'                               label = " (1x1)")
#' plot_diagnostics(bind_comparisons(a, b))
#' }
bind_comparisons <- function(...) {
  fits <- list(...)
  formulas <- unlist(lapply(fits, `[[`, "formulas"))
  if (anyDuplicated(names(formulas))) {
    stop("Model names are not unique; use `label` in standardize_tuna_indices().",
         call. = FALSE)
  }
  stack <- function(what) {
    out <- do.call(rbind, lapply(fits, function(f) {
      x <- f[[what]]
      x$model <- as.character(x$model)
      x
    }))
    out$model <- factor(out$model, levels = names(formulas))
    out
  }
  list(indices = stack("indices"), diagnostics = stack("diagnostics"),
       formulas = formulas, data = fits[[1]]$data, pca = fits[[1]]$pca)
}

#' Keep only some models of a standardisation result
#'
#' Useful to reuse an expensive fit (e.g. a 1 x 1 spatial control) and attach
#' it to new fits with [bind_comparisons()] instead of refitting it.
#'
#' @param comparison Output of [standardize_tuna_indices()] or
#'   [bind_comparisons()].
#' @param models Model names to keep.
#' @return The same structure, restricted to `models`.
#' @export
#' @examples
#' \dontrun{
#' control <- select_models(old_fit, "Base+S (1x1)")
#' plot_indices(bind_comparisons(new_fit, control))
#' }
select_models <- function(comparison, models) {
  missing <- setdiff(models, names(comparison$formulas))
  if (length(missing)) stop("Models not found: ", paste(missing, collapse = ", "),
                            call. = FALSE)
  keep <- function(x) {
    x <- x[x$model %in% models, , drop = FALSE]
    x$model <- factor(as.character(x$model), levels = models)
    x
  }
  comparison$indices <- keep(comparison$indices)
  comparison$diagnostics <- keep(comparison$diagnostics)
  comparison$formulas <- comparison$formulas[models]
  comparison
}
