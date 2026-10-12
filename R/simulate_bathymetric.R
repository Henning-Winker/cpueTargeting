#' Default settings for the four-species bathymetric scenario
#'
#' Returns the Run 22 settings from the original `OMcpueSIM22.R` example.
#' Species distributions and movement are associated with bathymetry and
#' horizontal distance; the model has no vertical availability component.
#'
#' @return A named list accepted by [simulate_bathymetric()].
#' @export
bathymetric_settings <- function() {
  list(
    n_years = 30, n_species = 4,
    price_fish = matrix(c(
      rep(c(1.5, 1, 1, 0.5), 10),
      rep(c(0.5, 1.5, 1, 1), 10),
      rep(c(0.5, 1, 0.5, 2), 10)
    ), ncol = 4, byrow = TRUE),
    func_mvt_dist = rep("Exp", 4),
    func_mvt_depth = rep("Lognorm", 4),
    func_mvt_lat = rep("Unif", 4),
    Fish_dist_par1 = c(0.4, 0.4, 0.4, 1),
    Fish_dist_par2 = rep(0, 4),
    Fish_depth_par1 = c(220, 550, 250, 1000),
    Fish_depth_par2 = c(0.5, 0.5, 0.3, 0.5),
    Fish_range_par1 = rep(0, 4), Fish_range_par2 = rep(50, 4),
    B0 = c(320000, 420000, 56000, 540000),
    r = c(0.3, 0.17, 0.12, 0.08), no.mvt = FALSE,
    Nvessels = 30, Tot_effort = 500, CV_effort = 0.2,
    qq_original = 0.1, CV_vessel = 0.2, do.tweedie = TRUE,
    xi = 1.2, phi = 0.1, Preference = 1,
    Changing_preference = FALSE,
    Range_X = 1:50, Range_Y = 1:50, SD_O = 150, SpatialScale = 4,
    year.depth.restriction = 1000, depth.restriction = c(0, 1),
    plotting = FALSE, Interactive = FALSE
  )
}

#' Simulate the four-species bathymetric Run 22 scenario
#'
#' Reproduces the spatial structure of the original Run 22 operating model:
#' four species are distributed across a bathymetric landscape, and fishing
#' effort responds to expected multispecies catch. The original files are
#' preserved in `inst/extdata`. This model represents horizontal bathymetric
#' habitat only and does not include the vertical tuna availability process.
#'
#' @param settings A named settings list, usually from
#'   [bathymetric_settings()].
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with `Data` (set-by-species observations), `Biomass` (annual
#'   relative nominal and true biomass indices), `Biomass.XY` (spatial density
#'   data), and `Plots` (ggplot diagnostics).
#' @import ggplot2
#' @import grDevices
#' @import stats
#' @export
#'
#' @examples
#' \dontrun{
#' om <- simulate_bathymetric(seed = 42)
#' head(om$Data)
#' om$Plots$Depth
#' }
simulate_bathymetric <- function(settings = bathymetric_settings(), seed = NULL) {
  if (!is.list(settings) || is.null(names(settings))) {
    stop("`settings` must be a named list.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  generate_bathymetric_scenario_data(settings)
}
