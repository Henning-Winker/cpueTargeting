#' Compare tuna CPUE indices with and without HBF adjustment
#'
#' Fits base, spatial, direct-principal-component, cluster, spatial-plus-HBF,
#' and spatial-plus-DPC-plus-HBF models to a tuna simulation. The latent actual
#' hook depth and vertical overlap are not used as predictors. Model indices
#' are standardized to the same reference settings and returned with the
#' operating-model biomass index for comparison.
#'
#' @param simulation Output from [simulate_tuna()].
#' @param species_names Labels for the four species, in simulator order.
#' @param seed Optional seed for the fishing-tactic k-means classification.
#'
#' @return A list with `indices` (annual relative indices and the true biomass
#'   index) and `diagnostics` (mean squared log error by species and model).
#' @export
#'
#' @examples
#' \dontrun{
#' om <- simulate_tuna(seed = 42)
#' comparison <- standardize_tuna_indices(om)
#' }
standardize_tuna_indices <- function(
    simulation, species_names = c("BET", "YFT", "ALB", "SWO"), seed = 22062017) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Install the suggested `mgcv` package to fit CPUE models.",
         call. = FALSE)
  }
  if (is.null(simulation$Data) || is.null(simulation$Biomass)) {
    stop("`simulation` must be the output from `simulate_tuna()`.",
         call. = FALSE)
  }
  n_species <- length(species_names)
  set_data <- simulation$Data
  d <- reshape2::dcast(
    set_data,
    SetID + year + vessel + Xobs + Yobs + HBF ~ Sp,
    value.var = "CPUE", fun.aggregate = mean
  )
  sp_cols <- paste0("Sp", seq_len(n_species))
  if (!all(sp_cols %in% names(d))) {
    stop("Could not find the expected species catch columns in simulation data.",
         call. = FALSE)
  }
  catch <- as.matrix(d[, sp_cols, drop = FALSE])
  catch[!is.finite(catch)] <- 0
  composition <- catch + 1e-8
  composition <- composition / rowSums(composition)
  pca_data <- as.data.frame(composition^(1 / 4))
  pca <- stats::prcomp(pca_data, center = TRUE, scale. = TRUE)
  pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  names(pcs) <- c("PC1", "PC2")
  d <- cbind(d, pcs)
  if (!is.null(seed)) set.seed(seed)
  d$FT <- factor(stats::kmeans(pca_data, centers = 3, nstart = 20)$cluster)
  years <- sort(unique(d$year))
  d$year <- factor(d$year, levels = years)
  d$vessel <- factor(d$vessel)

  forms <- list(
    Base = y ~ year,
    Spatial = y ~ year + mgcv::s(Xobs, Yobs, k = 8),
    DPC = y ~ year + mgcv::s(Xobs, Yobs, k = 8) +
      mgcv::s(PC1, k = 3) + mgcv::s(PC2, k = 3),
    Cluster = y ~ year + mgcv::s(Xobs, Yobs, k = 8) + FT,
    Spatial_HBF = y ~ year + mgcv::s(Xobs, Yobs, k = 8) + HBF,
    Spatial_DPC_HBF = y ~ year + mgcv::s(Xobs, Yobs, k = 8) + HBF +
      mgcv::s(PC1, k = 3) + mgcv::s(PC2, k = 3)
  )
  prediction <- data.frame(
    year = factor(years, levels = levels(d$year)),
    Xobs = mean(d$Xobs), Yobs = mean(d$Yobs), HBF = 8,
    PC1 = 0, PC2 = 0, FT = levels(d$FT)[1]
  )

  results <- list()
  for (s in seq_len(n_species)) {
    d$y <- d[[sp_cols[s]]]
    for (model_name in names(forms)) {
      fit <- mgcv::gam(forms[[model_name]], data = d,
        family = mgcv::Tweedie(p = 1.2), method = "REML")
      pred <- stats::predict(fit, newdata = prediction,
                             type = "link", se.fit = TRUE)
      raw <- exp(as.numeric(pred$fit))
      truth_rows <- simulation$Biomass$Sp == paste0("Sp", s)
      truth <- simulation$Biomass$True[truth_rows]
      if (length(truth) != length(years)) {
        stop("True biomass index does not align with simulated years.",
             call. = FALSE)
      }
      results[[length(results) + 1L]] <- data.frame(
        species = species_names[s], year = as.integer(as.character(years)),
        model = model_name, index = raw / mean(raw),
        lower = exp(as.numeric(pred$fit) - 1.96 * as.numeric(pred$se.fit)) /
          mean(raw),
        upper = exp(as.numeric(pred$fit) + 1.96 * as.numeric(pred$se.fit)) /
          mean(raw),
        true = truth
      )
    }
  }
  indices <- do.call(rbind, results)
  diagnostics <- stats::aggregate(
    (log(index + 1e-8) - log(true + 1e-8))^2 ~ species + model,
    data = indices, FUN = mean
  )
  names(diagnostics)[3] <- "mean_squared_log_error"
  list(indices = indices, diagnostics = diagnostics)
}
