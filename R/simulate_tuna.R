#' Default settings for the illustrative tuna longline scenario
#'
#' Returns the settings used by the current tuna example. The vertical habitat
#' and HBF-to-depth parameters are illustrative placeholders and need
#' fishery-specific calibration before scientific interpretation.
#'
#' @return A named list accepted by [simulate_tuna()].
#' @export
tuna_settings <- function() {
  list(
    n_years = 30, n_species = 4,
    price_fish = matrix(c(
      rep(c(1.5, 1.0, 1.0, 0.5), 10),
      rep(c(0.5, 1.5, 1.0, 1.0), 10),
      rep(c(0.5, 1.0, 0.5, 2.0), 10)
    ), ncol = 4, byrow = TRUE),
    func_mvt_dist = rep("Exp", 4), func_mvt_depth = rep("Lognorm", 4),
    func_mvt_lat = rep("Unif", 4),
    Fish_dist_par1 = c(0.4, 0.4, 0.5, 0.6), Fish_dist_par2 = rep(0, 4),
    Fish_depth_par1 = c(550, 250, 400, 450),
    Fish_depth_par2 = rep(0.5, 4),
    Fish_range_par1 = rep(0, 4), Fish_range_par2 = rep(26, 4),
    Fish_vertical_mu = c(260, 140, 200, 210),
    Fish_vertical_sd = c(90, 75, 90, 105),
    HBF_levels = c(5, 8, 11),
    HBF_depth_mean = c("5" = 125, "8" = 185, "11" = 245),
    HBF_depth_sd = c("5" = 30, "8" = 42, "11" = 58),
    HBF_hook_sd = c("5" = 45, "8" = 60, "11" = 78),
    HBF_choice_temp = 0.20,
    B0 = c(320000, 420000, 260000, 350000),
    r = c(0.16, 0.18, 0.12, 0.10), no.mvt = FALSE, Nvessels = 30,
    Tot_effort = 500, CV_effort = 0.2, qq_original = 0.1,
    CV_vessel = 0.2, do.tweedie = TRUE, xi = 1.2, phi = 0.1,
    Preference = 1.25, Changing_preference = FALSE,
    Range_X = 1:25, Range_Y = 1:25, SD_O = 150, SpatialScale = 4,
    year.depth.restriction = 1000, depth.restriction = c(0, 1),
    plotting = FALSE, Interactive = FALSE
  )
}

#' Simulate spatially distributed tuna with longline hook-depth overlap
#'
#' Generates the spatial operating model using the tuna settings, species
#' vertical profiles and noisy HBF-to-hook-depth relationship. Bathymetry
#' governs horizontal habitat association; the separate vertical profiles
#' govern overlap between fish and longline hooks.
#'
#' @param settings A named settings list, usually from [tuna_settings()].
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with `Data` (set-by-species observations, including HBF and
#'   latent hook-depth/overlap truth), `Biomass`, spatial biomass information,
#'   and diagnostic `Plots` created by the operating model.
#' @import ggplot2
#' @import grDevices
#' @import stats
#' @export
#'
#' @examples
#' \dontrun{
#' om <- simulate_tuna(seed = 42)
#' head(om$Data)
#' }
simulate_tuna <- function(settings = tuna_settings(), seed = NULL) {
  if (!is.list(settings)) stop("`settings` must be a named list.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  generate_tuna_scenario_data(settings)
}
