jet_colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#' Diagnostic plots for an operating-model run
#'
#' @param sim Output of [simulate_bathymetric()] or [simulate_tuna()].
#' @return A named list of ggplots: `Env` (bathymetry, or SST for the tuna
#'   model), `Thermocline` (fixed thermocline anomaly; tuna model only),
#'   `SpDens` (initial densities), `Hist` (CPUE
#'   frequencies), `CPUE` (nominal versus true index) and, for the tuna model,
#'   `HBF` (annual use of each HBF level) and `Avail` (see
#'   [plot_availability()]). `Depth` is kept as an alias of
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
    # separate maps so that each field has its own colour scale
    plots$Env <- raster(env, "SST", "SST (north = top): horizontal habitat") +
      labs(fill = "SST (deg C)")
    plots$Thermocline <- ggplot(env, aes(x = .data$x, y = .data$y,
                                         fill = .data$thermocline)) +
      geom_raster() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                           name = "Anomaly (m)\n+ = deeper") +
      coord_quickmap(expand = FALSE) +
      labs(title = "Fixed thermocline anomaly: shifts species depth", x = "x", y = "y") +
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
    plots$Avail <- plot_availability(sim)
  }
  plots
}

#' Plot species vertical profiles and hook-depth profiles by HBF
#'
#' One panel per HBF level. Grey bars show where the hooks of a set with that
#' HBF fish (the expected hook-depth distribution, drawn like Fig. 2 of
#' Bigelow et al. 2002); the coloured translucent areas show where each
#' species lives (thermocline anomaly of zero). Where a species' area overlaps
#' the grey bars, that set can catch it. The number next to each species is
#' its vertical availability to that HBF, i.e. the overlap that multiplies
#' catchability in the operating model (see [plot_availability()]).
#'
#' @param settings Settings list, usually from [tuna_settings()].
#' @param overlay Logical; `TRUE` (default) overlays the species on each HBF
#'   panel; `FALSE` shows the species in a separate first panel.
#' @param file Optional file path; if supplied, a 300-dpi PNG is saved.
#' @return A ggplot object.
#' @export
#' @examples
#' plot_tuna_vertical_distributions()
plot_tuna_vertical_distributions <- function(settings = tuna_settings(), overlay = TRUE,
                                             file = NULL) {
  s <- utils::modifyList(tuna_settings(), settings)
  edges <- seq(0, s$max_depth, by = s$depth_bin)
  mids <- edges[-1] - s$depth_bin / 2
  hbf_panels <- paste("HBF", s$HBF_levels)
  panels <- if (overlay) hbf_panels else c("Species", hbf_panels)
  fish_prof <- fish_profiles(0, s, edges)
  hooks <- nominal_hook_profiles(s, edges)

  fish_one <- function(panel) {
    out <- do.call(rbind, lapply(seq_along(s$species), function(i) {
      data.frame(depth = mids, percent = 100 * fish_prof[[i]][1, ],
                 species = s$species[i], panel = panel)
    }))
    out
  }
  fish <- if (overlay) do.call(rbind, lapply(hbf_panels, fish_one)) else fish_one("Species")
  fish$species <- factor(fish$species, levels = s$species)
  fish$panel <- factor(fish$panel, levels = panels)
  gear <- do.call(rbind, lapply(seq_along(s$HBF_levels), function(h) {
    data.frame(depth = mids, percent = 100 * hooks[h, ], panel = hbf_panels[h])
  }))
  gear$panel <- factor(gear$panel, levels = panels)

  # vertical availability of each species to each HBF (as in the OM)
  avail <- do.call(rbind, lapply(seq_along(s$HBF_levels), function(h) {
    a <- vapply(fish_prof, function(f) sum(hooks[h, ] * f[1, ]), numeric(1)) *
      s$avail_ref / s$depth_bin
    data.frame(panel = hbf_panels[h], species = s$species,
               label = sprintf("%s %.2f", s$species, a))
  }))
  avail$panel <- factor(avail$panel, levels = panels)
  avail$species <- factor(avail$species, levels = s$species)
  x_max <- max(c(fish$percent, gear$percent))
  avail$x <- x_max
  avail$y <- s$max_depth * (0.62 + 0.08 * (as.integer(avail$species) - 1))

  p <- ggplot() +
    geom_col(data = gear, aes(x = .data$percent, y = .data$depth),
             orientation = "y", width = s$depth_bin, fill = "grey55", colour = NA)
  if (overlay) {
    p <- p +
      geom_ribbon(data = fish, aes(y = .data$depth, xmin = 0, xmax = .data$percent,
                                   fill = .data$species),
                  orientation = "y", alpha = 0.35) +
      geom_path(data = fish, aes(x = .data$percent, y = .data$depth,
                                 colour = .data$species), linewidth = 0.7) +
      geom_text(data = avail, aes(x = .data$x, y = .data$y, label = .data$label,
                                  colour = .data$species),
                hjust = 1, size = 3, show.legend = FALSE) +
      labs(title = "Where the hooks fish (grey) and where the fish are (colour)",
           subtitle = "Numbers: vertical availability of each species to that HBF")
  } else {
    p <- p + geom_path(data = fish, aes(x = .data$percent, y = .data$depth,
                                        colour = .data$species), linewidth = 1)
  }
  p <- p +
    facet_wrap(~panel, nrow = 1) +
    scale_y_reverse() +
    labs(x = paste0("Frequency (% per ", s$depth_bin, " m)"), y = "Depth (m)",
         colour = NULL, fill = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  if (!is.null(file)) ggsave(file, p, width = 11, height = 5, dpi = 300)
  p
}

#' Plot standardised indices against the true biomass trend
#'
#' @param comparison Output of [standardize_tuna_indices()], or a named list of
#'   such outputs (one per scenario), which are shown as rows.
#' @param models Optional subset of model names to show.
#' @param ci Logical; show 95% confidence intervals.
#' @param type How to draw the indices: `"lines"` (default) joins the annual
#'   values; `"smooth"` shows only a loess smoother of each index, which makes
#'   trends easier to compare; `"points+smooth"` shows the annual values as
#'   points with the smoother on top.
#' @param span Smoothing span of the loess smoother (larger = smoother).
#' @return A ggplot object.
#' @export
plot_indices <- function(comparison, models = NULL, ci = FALSE,
                         type = c("lines", "smooth", "points+smooth"), span = 0.5) {
  type <- match.arg(type)
  idx <- bind_scenarios(comparison, "indices")
  if (!is.null(models)) idx <- idx[idx$model %in% models, ]
  truth <- idx[!duplicated(idx[, c("scenario", "species", "year")]), ]
  p <- ggplot(idx, aes(x = .data$year, y = .data$index, colour = .data$model))
  if (ci) {
    p <- p + geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper, fill = .data$model),
                         alpha = 0.15, colour = NA, show.legend = FALSE)
  }
  smoother <- geom_smooth(method = "loess", formula = y ~ x, span = span, se = FALSE,
                          linewidth = 0.9)
  p <- p + switch(type,
    lines = list(geom_line()),
    smooth = list(smoother),
    `points+smooth` = list(geom_point(size = 1, alpha = 0.6), smoother))
  p <- p +
    geom_line(data = truth, mapping = aes(x = .data$year, y = .data$true),
              colour = "black", linewidth = 1.1, inherit.aes = FALSE) +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(x = "Year", y = "Relative index (true: black)", colour = "Model") +
    theme_bw()
  if (length(unique(idx$scenario)) > 1) {
    p + facet_grid(scenario ~ species, scales = "free_y")
  } else {
    p + facet_wrap(~species, scales = "free_y")
  }
}

#' Compare the accuracy of standardised indices
#'
#' Bar chart of the root mean squared error of the log index (relative to the
#' true biomass index) by model and species.
#'
#' @param comparison Output of [standardize_tuna_indices()], or a named list of
#'   such outputs (one per scenario), shown as differently coloured bars.
#' @param models Optional subset of model names to show.
#' @return A ggplot object.
#' @export
plot_diagnostics <- function(comparison, models = NULL) {
  diag <- bind_scenarios(comparison, "diagnostics")
  if (!is.null(models)) diag <- diag[diag$model %in% models, ]
  ggplot(diag, aes(x = .data$model, y = .data$rmse_log, fill = .data$scenario)) +
    geom_col(position = position_dodge(preserve = "single"), width = 0.75) +
    facet_wrap(~species) +
    labs(x = NULL, y = "RMSE of log index", fill = "Scenario") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Stack a component of one or several comparison results
#' @noRd
bind_scenarios <- function(comparison, what) {
  if (!is.null(comparison[[what]])) comparison <- list(" " = comparison)
  if (is.null(names(comparison))) names(comparison) <- paste("Scenario", seq_along(comparison))
  out <- do.call(rbind, lapply(names(comparison), function(n) {
    data.frame(scenario = n, comparison[[n]][[what]])
  }))
  out$scenario <- factor(out$scenario, levels = names(comparison))
  out
}

#' Plot the thermocline anomaly
#'
#' @param sim Output of [simulate_tuna()].
#' @param years Years to map (default: six years spread over the series).
#' @return A list of two ggplots: `series` (mean anomaly by year with the
#'   10-90% range across cells) and `maps` (anomaly fields in selected years).
#' @export
plot_thermocline <- function(sim, years = NULL) {
  th <- sim$Thermocline
  if (is.null(th)) stop("`sim` must be the output of simulate_tuna().", call. = FALSE)
  stats_by_year <- function(f) tapply(th$anomaly, th$year, f)
  series <- data.frame(year = sort(unique(th$year)),
                       mean = stats_by_year(mean),
                       lower = stats_by_year(function(a) stats::quantile(a, 0.1)),
                       upper = stats_by_year(function(a) stats::quantile(a, 0.9)))
  p_series <- ggplot(series, aes(x = .data$year)) +
    geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper), fill = "grey80") +
    geom_line(aes(y = .data$mean), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_reverse() +
    labs(x = "Year", y = "Thermocline anomaly (m, deeper = down)",
         title = "Thermocline: basin mean and 10-90% range across cells") +
    theme_bw()
  if (is.null(years)) years <- unique(round(seq(1, max(th$year), length.out = 6)))
  p_maps <- ggplot(th[th$year %in% years, ],
                   aes(x = .data$x, y = .data$y, fill = .data$anomaly)) +
    geom_raster() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                         name = "Anomaly (m)\n+ = deeper") +
    facet_wrap(~year, labeller = label_both) +
    coord_quickmap(expand = FALSE) +
    labs(title = "Thermocline anomaly by year") +
    theme_bw()
  list(series = p_series, maps = p_maps)
}

#' Show how a thermocline shift moves species relative to the hooks
#'
#' Species depth profiles for a shallow, average and deep thermocline, with
#' the expected hook-depth profiles of each HBF level in grey. When the
#' thermocline is shallow, deep-set hooks fish below the deep species; when
#' it is deep, shallow sets miss them.
#'
#' @param settings Settings list, usually from [tuna_settings()].
#' @param shifts Thermocline anomalies to show (m, positive = deeper).
#' @return A ggplot object.
#' @export
#' @examples
#' plot_thermocline_effect()
plot_thermocline_effect <- function(settings = tuna_settings(), shifts = c(-40, 0, 40)) {
  s <- utils::modifyList(tuna_settings(), settings)
  edges <- seq(0, s$max_depth, by = s$depth_bin)
  mids <- edges[-1] - s$depth_bin / 2
  fish <- do.call(rbind, lapply(shifts, function(sh) {
    prof <- fish_profiles(sh, s, edges)
    do.call(rbind, lapply(seq_along(s$species), function(i) {
      data.frame(depth = mids, proportion = prof[[i]][1, ], species = s$species[i],
                 shift = sprintf("%+d m", as.integer(sh)))
    }))
  }))
  fish$species <- factor(fish$species, levels = s$species)
  fish$shift <- factor(fish$shift, levels = unique(fish$shift))
  hooks <- nominal_hook_profiles(s, edges)
  gear <- do.call(rbind, lapply(seq_along(s$HBF_levels), function(h) {
    data.frame(depth = mids, proportion = hooks[h, ], HBF = paste("HBF", s$HBF_levels[h]))
  }))
  label_pos <- do.call(rbind, lapply(split(gear, gear$HBF), function(g) g[which.max(g$proportion), ]))
  ggplot(fish, aes(x = .data$proportion, y = .data$depth)) +
    geom_path(data = gear, aes(group = .data$HBF), colour = "grey60", linetype = 2) +
    geom_text(data = label_pos, aes(label = .data$HBF), colour = "grey40", size = 3,
              hjust = -0.1) +
    geom_path(aes(colour = .data$shift), linewidth = 1) +
    facet_wrap(~species, nrow = 1) +
    scale_y_reverse() +
    scale_colour_manual(values = shift_colours(length(shifts)),
                        name = "Thermocline\nanomaly") +
    labs(x = paste0("Proportion per ", s$depth_bin, "-m bin"), y = "Depth (m)") +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Colours from shallow (orange) through black to deep (blue)
#' @noRd
shift_colours <- function(n) {
  if (n == 1) return("black")
  grDevices::colorRampPalette(c("#D55E00", "black", "#0072B2"))(n)
}

#' Vertical availability by HBF over time
#'
#' Mean realised vertical availability (overlap between hooks and fish) of the
#' sets fishing with each HBF level, by year and species. With a fixed
#' thermocline the lines are flat apart from changes in where the fleet
#' fishes; with a moving thermocline the catchability of a given HBF changes
#' from year to year, which an HBF factor cannot absorb.
#'
#' @param sim Output of [simulate_tuna()].
#' @return A ggplot object.
#' @export
plot_availability <- function(sim) {
  d <- sim$Data
  if (!"vertical_avail" %in% names(d)) stop("`sim` must be the output of simulate_tuna().",
                                            call. = FALSE)
  avail <- stats::aggregate(list(avail = d$vertical_avail),
                            by = list(year = d$year, HBF = d$HBF, Sp = d$Sp), FUN = mean)
  avail$Sp <- factor(avail$Sp, levels = unique(d$Sp))
  avail$HBF <- factor(avail$HBF)
  ggplot(avail, aes(x = .data$year, y = .data$avail, colour = .data$HBF)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~Sp, nrow = 1) +
    labs(x = "Year", y = "Mean vertical availability", colour = "HBF") +
    theme_bw()
}

#' What do the principal components of the catch composition capture?
#'
#' Plots the PC scores used by the DPC models (Winker et al. 2013, 2014) in
#' four ways: the species loadings (what each PC means), mean scores by cell
#' (where each tactic is used), mean scores against the realised hook depth of
#' the set for each HBF (how the PCs track vertical targeting; tuna model
#' only), and mean scores by year (how targeting changes over time).
#'
#' `hook_depth` is latent truth from the operating model; it is used here only
#' to show what the PCs pick up, not in any standardisation model.
#'
#' @param comparison Output of [standardize_tuna_indices()] (or of
#'   [bind_comparisons()]).
#' @param pcs Principal components to show.
#' @param depth_bin Width of the hook-depth bins (m).
#' @param min_sets Minimum number of sets for a depth bin to be shown.
#' @return A list of ggplots: `loadings`, `space`, `depth` (tuna model only)
#'   and `time`.
#' @export
plot_pcs <- function(comparison, pcs = c("PC1", "PC2"), depth_bin = 20, min_sets = 30) {
  d <- comparison$data
  missing <- setdiff(pcs, names(d))
  if (length(missing)) stop("Not found in the data: ", paste(missing, collapse = ", "),
                            call. = FALSE)
  long <- do.call(rbind, lapply(pcs, function(pc) {
    data.frame(d[, intersect(c("X", "Y", "year", "HBF", "hook_depth"), names(d))],
               pc = pc, score = d[[pc]])
  }))
  long$pc <- factor(long$pc, levels = pcs)
  plots <- list()

  # Loadings: what each PC means in terms of species
  if (!is.null(comparison$pca)) {
    load <- comparison$pca$loadings[, pcs, drop = FALSE]
    var_share <- comparison$pca$variance / sum(comparison$pca$variance)
    names(var_share) <- colnames(comparison$pca$loadings)
    labels <- sprintf("%s (%.0f%%)", pcs, 100 * var_share[pcs])
    loadings <- data.frame(species = rep(rownames(load), length(pcs)),
                           pc = factor(rep(labels, each = nrow(load)), levels = labels),
                           loading = as.vector(load))
    loadings$species <- factor(loadings$species, levels = rownames(load))
    plots$loadings <- ggplot(loadings, aes(x = .data$species, y = .data$loading,
                                           fill = .data$loading > 0)) +
      geom_col(show.legend = FALSE) +
      geom_hline(yintercept = 0) +
      facet_wrap(~pc) +
      scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00")) +
      labs(x = NULL, y = "Loading", title = "What the PCs mean") +
      theme_bw()
  }

  # Space: mean score by cell
  space <- stats::aggregate(list(score = long$score),
                            by = list(x = long$X, y = long$Y, pc = long$pc), FUN = mean)
  plots$space <- ggplot(space, aes(x = .data$x, y = .data$y, fill = .data$score)) +
    geom_raster() +
    scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2",
                         name = "Mean score") +
    facet_wrap(~pc) +
    coord_quickmap(expand = FALSE) +
    labs(x = "x", y = "y (north = top)", title = "PC scores in space") +
    theme_bw()

  # Depth: mean score against realised hook depth, by HBF
  if (all(c("hook_depth", "HBF") %in% names(long))) {
    long$depth <- (floor(long$hook_depth / depth_bin) + 0.5) * depth_bin
    by_depth <- stats::aggregate(list(score = long$score),
                                 by = list(depth = long$depth, HBF = long$HBF, pc = long$pc),
                                 FUN = mean)
    n_sets <- stats::aggregate(list(n = long$score),
                               by = list(depth = long$depth, HBF = long$HBF, pc = long$pc),
                               FUN = length)
    by_depth <- by_depth[n_sets$n >= min_sets, , drop = FALSE]
    by_depth <- by_depth[order(by_depth$pc, by_depth$HBF, by_depth$depth), ]
    plots$depth <- ggplot(by_depth, aes(x = .data$score, y = .data$depth,
                                        colour = factor(.data$HBF))) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_path(linewidth = 0.8) +
      geom_point(size = 1.2) +
      facet_wrap(~pc) +
      scale_y_reverse() +
      labs(x = "Mean score", y = "Realised median hook depth of the set (m)",
           colour = "HBF", title = "PC scores against hook depth") +
      theme_bw()
  }

  # Time: mean score by year
  by_year <- stats::aggregate(list(score = long$score),
                              by = list(year = as.integer(as.character(long$year)),
                                        pc = long$pc),
                              FUN = mean)
  plots$time <- ggplot(by_year, aes(x = .data$year, y = .data$score)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    facet_wrap(~pc) +
    labs(x = "Year", y = "Mean score", title = "PC scores over time") +
    theme_bw()
  plots
}

#' Species distributions in three dimensions: sections and depth layers
#'
#' Combines the horizontal density of each species with its vertical profile
#' in each cell (including the thermocline anomaly of the chosen year) and
#' shows the result in two ways:
#' * `section`: a latitude x depth section (density averaged over x), with
#'   depth downwards and north to the right;
#' * `layers`: maps of the density in depth layers, species by layer.
#'
#' Horizontal densities are the initial (unfished) distribution, so the plots
#' show where each species lives rather than the effect of fishing.
#'
#' @param sim Output of [simulate_tuna()].
#' @param year Year whose thermocline anomaly is used (matters only when the
#'   thermocline moves between years).
#' @param layers Depth-layer boundaries (m) for the layer maps.
#' @return A list of two ggplots, `section` and `layers`.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_tuna(seed = 1)
#' p <- plot_sections(sim)
#' p$section
#' p$layers
#' }
plot_sections <- function(sim, year = 1, layers = c(0, 100, 200, 300, 500)) {
  if (is.null(sim$Thermocline)) stop("`sim` must be the output of simulate_tuna().",
                                     call. = FALSE)
  s <- sim$settings
  env <- sim$Env
  edges <- seq(0, s$max_depth, by = s$depth_bin)
  mids <- edges[-1] - s$depth_bin / 2
  th <- sim$Thermocline
  anomaly <- th$anomaly[th$year == year][match(paste(env$x, env$y),
                                               paste(th$x[th$year == year], th$y[th$year == year]))]
  profiles <- fish_profiles(anomaly, s, edges)
  dens <- sim$Biomass.XY

  section <- do.call(rbind, lapply(seq_along(s$species), function(i) {
    b <- dens$Density[dens$variable == s$species[i]]
    by_lat <- rowsum(b * profiles[[i]], env$y) / as.vector(table(env$y))
    data.frame(y = rep(as.numeric(rownames(by_lat)), ncol(by_lat)),
               depth = rep(mids, each = nrow(by_lat)),
               density = as.vector(by_lat) / max(by_lat),
               species = s$species[i])
  }))
  section$species <- factor(section$species, levels = s$species)
  p_section <- ggplot(section, aes(x = .data$y, y = .data$depth, fill = .data$density)) +
    geom_raster() +
    scale_fill_gradientn(colours = jet_colors(10), name = "Relative\ndensity") +
    scale_y_reverse(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    facet_wrap(~species, nrow = 1) +
    labs(x = "y (south to north)", y = "Depth (m)",
         title = paste("Latitude x depth sections, year", year)) +
    theme_bw()

  layer_bins <- cut(mids, breaks = layers, include.lowest = TRUE,
                    labels = paste0(utils::head(layers, -1), "-", layers[-1], " m"))
  maps <- do.call(rbind, lapply(seq_along(s$species), function(i) {
    b <- dens$Density[dens$variable == s$species[i]]
    share <- vapply(levels(layer_bins), function(l) {
      rowSums(profiles[[i]][, layer_bins == l, drop = FALSE])
    }, numeric(nrow(env)))
    d <- b * share
    data.frame(x = env$x, y = env$y,
               layer = factor(rep(levels(layer_bins), each = nrow(env)),
                              levels = levels(layer_bins)),
               density = as.vector(d) / max(d), species = s$species[i])
  }))
  maps$species <- factor(maps$species, levels = s$species)
  p_layers <- ggplot(maps, aes(x = .data$x, y = .data$y, fill = .data$density)) +
    geom_raster() +
    scale_fill_gradientn(colours = jet_colors(10), name = "Relative\ndensity") +
    facet_grid(layer ~ species) +
    coord_quickmap(expand = FALSE) +
    labs(x = "x", y = "y (north = top)",
         title = paste("Density by depth layer, year", year)) +
    theme_bw()

  list(section = p_section, layers = p_layers)
}
