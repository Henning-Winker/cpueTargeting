#' Standardise simulated CPUE with alternative targeting covariates
#'
#' Fits Tweedie GAMs (`mgcv::gam` with `tw()`, as in the CPUE guidelines) to
#' each species and compares the annual indices with the true biomass trend.
#'
#' Candidate models (all include `year`):
#' * `Base`: year only.
#' * `Spatial`: + `s(Xobs, Yobs)` on the coarse reporting grid.
#' * `DPC`: Spatial + `s(PC1)` + `s(PC2)` from the catch composition
#'   (direct principal components, Winker et al. 2013, 2014).
#' * `Cluster`: Spatial + fishing tactic from k-means clustering of the
#'   catch composition.
#' * `HBF`: Spatial + HBF as a factor (tuna model only).
#' * `DPC_HBF`: DPC + HBF (tuna model only).
#' * `HBS`: Spatial + `offset(log(hbs_avail))`, i.e. habitat-based
#'   standardisation with the HBF-predicted vertical overlap as effective
#'   effort (Bigelow et al. 2002; tuna model only, assumes the habitat model
#'   is known).
#'
#' The latent `hook_depth` and realised `vertical_avail` are never used as
#' predictors.
#'
#' @param simulation Output of [simulate_tuna()] or [simulate_bathymetric()].
#' @param models Names of the models to fit (default: all that apply).
#' @param vessel_re Logical; add a random vessel effect `s(vessel, bs = "re")`.
#' @param n_clusters Number of fishing tactics for the `Cluster` model.
#' @param seed Seed for the k-means clustering.
#'
#' @return A list with `indices` (species, year, model, index with 95% CI and
#'   the true index), `diagnostics` (root mean squared error of the log index
#'   by species and model), and `data` (the set-level data used, with PCs and
#'   tactics).
#' @export
#'
#' @examples
#' \dontrun{
#' om <- simulate_tuna(seed = 42)
#' comparison <- standardize_tuna_indices(om)
#' comparison$diagnostics
#' plot_indices(comparison)
#' }
standardize_tuna_indices <- function(simulation, models = NULL, vessel_re = FALSE,
                                     n_clusters = 3, seed = 1) {
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
  d$PC1 <- pca$x[, 1]
  d$PC2 <- pca$x[, 2]
  if (!is.null(seed)) set.seed(seed)
  d$FT <- factor(stats::kmeans(composition, centers = n_clusters, nstart = 20)$cluster)
  years <- sort(unique(d$year))
  d$year <- factor(d$year, levels = years)
  d$vessel <- factor(d$vessel)
  if (has_hbf) d$HBF <- factor(d$HBF)

  spatial <- "s(Xobs, Yobs, k = 8)"
  dpc <- "s(PC1, k = 3) + s(PC2, k = 3)"
  terms <- list(Base = "year",
                Spatial = c("year", spatial),
                DPC = c("year", spatial, dpc),
                Cluster = c("year", spatial, "FT"))
  if (has_hbf) {
    terms$HBF <- c("year", spatial, "HBF")
    terms$DPC_HBF <- c("year", spatial, dpc, "HBF")
    terms$HBS <- c("year", spatial, "offset(log(hbs))")
  }
  if (vessel_re) terms <- lapply(terms, c, 's(vessel, bs = "re")')
  if (!is.null(models)) {
    unknown <- setdiff(models, names(terms))
    if (length(unknown)) stop("Unknown or inapplicable models: ",
                              paste(unknown, collapse = ", "), call. = FALSE)
    terms <- terms[models]
  }

  # Reference levels for prediction: only `year` changes
  mode_level <- function(f) names(which.max(table(f)))
  newdata <- data.frame(year = factor(years, levels = years),
                        Xobs = stats::median(d$Xobs), Yobs = stats::median(d$Yobs),
                        PC1 = 0, PC2 = 0,
                        FT = factor(mode_level(d$FT), levels = levels(d$FT)),
                        vessel = factor(levels(d$vessel)[1], levels = levels(d$vessel)),
                        hbs = 1)
  if (has_hbf) newdata$HBF <- factor(mode_level(d$HBF), levels = levels(d$HBF))

  results <- list()
  for (sp in species) {
    d$y <- d[[sp]]
    if (has_hbf) d$hbs <- pmax(d[[paste0("hbs_", sp)]], 1e-6)
    truth <- simulation$Biomass$True[simulation$Biomass$Sp == sp]
    for (m in names(terms)) {
      form <- stats::as.formula(paste("y ~", paste(terms[[m]], collapse = " + ")))
      fit <- mgcv::gam(form, data = d, family = mgcv::tw(), method = "REML")
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
  indices$species <- factor(indices$species, levels = species)
  indices$model <- factor(indices$model, levels = names(terms))
  err <- (log(indices$index) - log(indices$true))^2
  diagnostics <- stats::aggregate(list(rmse_log = err),
                                  by = list(species = indices$species, model = indices$model),
                                  FUN = function(e) sqrt(mean(e)))
  list(indices = indices, diagnostics = diagnostics, data = d)
}
