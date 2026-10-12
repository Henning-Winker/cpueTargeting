jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#' Diagnostic plots for an operating-model run
#'
#' @param sim Output of [simulate_bathymetric()] or [simulate_tuna()].
#' @return A named list of ggplots: `Env` (bathymetry, or SST and thermocline
#'   anomaly for the tuna model), `SpDens` (initial densities), `Hist` (CPUE
#'   frequencies), `CPUE` (nominal versus true index) and, for the tuna model,
#'   `HBF` (annual use of each HBF level). `Depth` is kept as an alias of
#'   `Env` for the bathymetric model.
#' @export
plot_om <- function(sim) {
  env <- sim$Env
  raster <- function(data, fill, title) {
    ggplot(data, aes(x = .data$x, y = .data$y, fill = .data[[fill]])) +
      geom_raster() +
      scale_fill_gradientn(colours = jet_colors(10)) +
      coord_quickmap(expand = FALSE) +
      labs(title = title, x = "x", y = "y") +
      theme_bw()
  }
  plots <- list()
  if ("Depth" %in% names(env)) {
    plots$Env <- raster(env, "Depth", "Bathymetry")
    plots$Depth <- plots$Env
  } else {
    long <- rbind(data.frame(x = env$x, y = env$y, value = env$SST, field = "SST (deg C)"),
                  data.frame(x = env$x, y = env$y, value = env$thermocline,
                             field = "Thermocline anomaly (m)"))
    plots$Env <- ggplot(long, aes(x = .data$x, y = .data$y, fill = .data$value)) +
      geom_raster() +
      scale_fill_gradientn(colours = jet_colors(10), name = NULL) +
      facet_wrap(~field) +
      coord_quickmap(expand = FALSE) +
      labs(title = "Environment (north = top)") +
      theme_bw()
  }
  dens <- sim$Biomass.XY
  dens$log_density <- log(dens$Density + 1)
  plots$SpDens <- raster(dens, "log_density", "True spatial densities") +
    labs(fill = "log(Density + 1)") +
    facet_wrap(~variable)
  plots$Hist <- ggplot(sim$Data, aes(x = .data$CPUE)) +
    geom_histogram(bins = 80, fill = "grey", colour = 1) +
    facet_wrap(~Sp, scales = "free") +
    labs(x = "CPUE", y = "Frequency") +
    theme_bw()
  plots$CPUE <- ggplot(sim$Biomass, aes(x = .data$year)) +
    geom_line(aes(y = .data$Nominal)) +
    geom_point(aes(y = .data$Nominal), shape = 21, fill = "grey") +
    geom_line(aes(y = .data$True), linewidth = 1, colour = "darkblue") +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~Sp, scales = "free") +
    labs(x = "Year", y = "Relative index (nominal: points, true: blue)") +
    theme_bw()
  if ("HBF" %in% names(sim$Data)) {
    sets <- sim$Data[!duplicated(sim$Data$SetID), c("year", "HBF")]
    plots$HBF <- ggplot(sets, aes(x = .data$year, fill = factor(.data$HBF))) +
      geom_bar(position = "fill") +
      labs(x = "Year", y = "Proportion of sets", fill = "HBF") +
      theme_bw()
  }
  plots
}

#' Plot species vertical profiles and hook-depth profiles by HBF
#'
#' Fish profiles are shown for a cell with zero thermocline anomaly; hook
#' profiles are the expected distributions for each HBF level at the mean
#' mainline angle without current shoaling.
#'
#' @param settings Settings list, usually from [tuna_settings()].
#' @param file Optional file path; if supplied, a 300-dpi PNG is saved.
#' @return A ggplot object.
#' @export
#' @examples
#' plot_tuna_vertical_distributions()
plot_tuna_vertical_distributions <- function(settings = tuna_settings(), file = NULL) {
  s <- utils::modifyList(tuna_settings(), settings)
  edges <- seq(0, s$max_depth, by = s$depth_bin)
  mids <- edges[-1] - s$depth_bin / 2
  fish <- do.call(rbind, lapply(seq_along(s$species), function(i) {
    p <- stats::dnorm(mids, s$depth_mu[i], s$depth_sd[i])
    data.frame(depth = mids, proportion = p / sum(p), profile = s$species[i],
               panel = "Species")
  }))
  hooks <- nominal_hook_profiles(s, edges)
  gear <- do.call(rbind, lapply(seq_along(s$HBF_levels), function(h) {
    data.frame(depth = mids, proportion = hooks[h, ],
               profile = paste("HBF", s$HBF_levels[h]), panel = "Hooks")
  }))
  profiles <- rbind(fish, gear)
  profiles$panel <- factor(profiles$panel, levels = c("Species", "Hooks"))
  profiles$profile <- factor(profiles$profile, levels = unique(profiles$profile))
  p <- ggplot(profiles, aes(x = .data$proportion, y = .data$depth,
                            colour = .data$profile)) +
    geom_path(linewidth = 1) +
    facet_wrap(~panel, scales = "free_x") +
    scale_y_reverse() +
    labs(x = paste0("Proportion per ", s$depth_bin, "-m bin"), y = "Depth (m)",
         colour = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  if (!is.null(file)) ggsave(file, p, width = 9, height = 6, dpi = 300)
  p
}

#' Plot standardised indices against the true biomass trend
#'
#' @param comparison Output of [standardize_tuna_indices()].
#' @param models Optional subset of model names to show.
#' @return A ggplot object.
#' @export
plot_indices <- function(comparison, models = NULL) {
  idx <- comparison$indices
  if (!is.null(models)) idx <- idx[idx$model %in% models, ]
  truth <- idx[!duplicated(idx[, c("species", "year")]), ]
  ggplot(idx, aes(x = .data$year, y = .data$index, colour = .data$model)) +
    geom_line() +
    geom_line(data = truth, mapping = aes(x = .data$year, y = .data$true),
              colour = "black", linewidth = 1.1, inherit.aes = FALSE) +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~species, scales = "free_y") +
    labs(x = "Year", y = "Relative index (true: black)", colour = "Model") +
    theme_bw()
}
