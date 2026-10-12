#' Plot the 2014 cluster simulation's biomass, effort, and catches
#'
#' Creates the main diagnostic plots for a [simulate_cluster()] result. All
#' returned plots are ggplot objects: annual biomass trajectories, annual
#' habitat effort allocation, and trip-level catch distributions by species.
#'
#' @param simulation Output from [simulate_cluster()].
#'
#' @return A named list of three ggplot objects: `biomass`, `effort`, and
#'   `catches`.
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- simulate_cluster(seed = 1)
#' plots <- plot_cluster_2014(sim)
#' plots$biomass
#' }
plot_cluster_2014 <- function(simulation) {
  if (is.null(simulation$data) || is.null(simulation$biomass) ||
      is.null(simulation$effort)) {
    stop("`simulation` must be the output from `simulate_cluster()`.",
         call. = FALSE)
  }

  to_long <- function(df, id, names_to, values_to) {
    cols <- setdiff(names(df), id)
    out <- data.frame(rep(df[[id]], length(cols)),
                      factor(rep(cols, each = nrow(df)), levels = cols),
                      unlist(df[cols], use.names = FALSE))
    names(out) <- c(id, names_to, values_to)
    out
  }
  biomass <- to_long(simulation$biomass, "year", "species", "biomass")
  p_biomass <- ggplot2::ggplot(biomass,
      ggplot2::aes(x = .data$year, y = .data$biomass, colour = .data$species,
                   group = .data$species)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(title = "True biomass trajectories", x = "Year",
                  y = "Biomass", colour = "Species") +
    ggplot2::theme_bw()

  effort <- to_long(simulation$effort, "year", "habitat", "share")
  p_effort <- ggplot2::ggplot(effort,
      ggplot2::aes(x = .data$year, y = .data$share, colour = .data$habitat,
                   group = .data$habitat)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Annual effort allocation", x = "Year",
                  y = "Proportion of effort", colour = "Habitat") +
    ggplot2::theme_bw()

  species <- setdiff(names(simulation$biomass), "year")
  catches <- to_long(simulation$data[, c("trip", species)], "trip", "species", "cpue")
  p_catches <- ggplot2::ggplot(catches, ggplot2::aes(x = log10(.data$cpue + 1))) +
    ggplot2::geom_histogram(bins = 40, fill = "grey65", colour = "white") +
    ggplot2::facet_wrap(~species, scales = "free_y") +
    ggplot2::labs(title = "Trip-level catch distributions",
                  x = expression(log[10](CPUE + 1)), y = "Trips") +
    ggplot2::theme_bw()

  list(biomass = p_biomass, effort = p_effort, catches = p_catches)
}
