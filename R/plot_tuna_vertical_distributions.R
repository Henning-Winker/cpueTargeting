#' Plot fish vertical profiles and HBF hook-depth profiles
#'
#' The fish curves show the scenario's species-specific vertical availability.
#' The hook curves show marginal distributions of hook depth for each HBF
#' level, combining between-set depth uncertainty and within-set hook spread.
#'
#' @param settings Settings list, usually returned by [tuna_settings()].
#' @param max_depth Maximum depth, in metres, shown in the plot.
#' @param file Optional file path. If supplied, saves a 300-dpi PNG and still
#'   returns the plot object.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' plot_tuna_vertical_distributions()
plot_tuna_vertical_distributions <- function(settings = tuna_settings(),
                                             max_depth = 500, file = NULL) {
  required <- c("Fish_vertical_mu", "Fish_vertical_sd", "HBF_levels",
                "HBF_depth_mean", "HBF_depth_sd", "HBF_hook_sd")
  if (!all(required %in% names(settings))) {
    stop("Settings are missing vertical-profile parameters: ",
         paste(setdiff(required, names(settings)), collapse = ", "),
         call. = FALSE)
  }
  if (max_depth <= 0) stop("`max_depth` must be positive.", call. = FALSE)
  depth <- seq(0, max_depth, length.out = 501)
  species <- c("BET", "YFT", "ALB", "SWO")
  fish <- do.call(rbind, lapply(seq_along(settings$Fish_vertical_mu), function(i) {
    data.frame(depth_m = depth,
      density = stats::dnorm(depth, settings$Fish_vertical_mu[i],
                             settings$Fish_vertical_sd[i]),
      profile = species[i], panel = "Species vertical availability")
  }))
  gear <- do.call(rbind, lapply(settings$HBF_levels, function(h) {
    key <- as.character(h)
    sd <- sqrt(settings$HBF_depth_sd[key]^2 + settings$HBF_hook_sd[key]^2)
    data.frame(depth_m = depth,
      density = stats::dnorm(depth, settings$HBF_depth_mean[key], sd),
      profile = paste("HBF", h), panel = "Marginal hook-depth distribution")
  }))
  profiles <- rbind(fish, gear)
  profiles$panel <- factor(profiles$panel,
    levels = c("Species vertical availability", "Marginal hook-depth distribution"))
  p <- ggplot2::ggplot(profiles,
      ggplot2::aes(x = density, y = depth_m, colour = profile)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::facet_grid(panel ~ ., scales = "free_x") +
    ggplot2::scale_y_reverse(limits = c(max_depth, 0)) +
    ggplot2::labs(title = "Illustrative tuna and longline vertical distributions",
      x = "Probability density per metre", y = "Depth (m)", colour = NULL,
      caption = "Hook profiles combine between-set HBF depth uncertainty and within-set hook spread.") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  if (!is.null(file)) ggplot2::ggsave(file, p, width = 9, height = 8, dpi = 300)
  p
}
