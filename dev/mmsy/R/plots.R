# =================================================================
# plots.R: stock status, F, effort and yield (dev-mmsy prototype)
# =================================================================

library(ggplot2)

# categorical slots (validated: slots 1-3 all-pairs, 1-6 adjacent)
ms_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
ms_ink  <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#8a8984",
             grid = "#e8e7e3", ref = "#8a8984")

# FAO 3-alpha codes for the logbook species
ms_codes <- c("Dentex dentex" = "DEC", "Pagellus erythrinus" = "PAC",
              "Pagellus bogaraveo" = "SBR", "Merluccius merluccius" = "HKE",
              "Conger conger" = "COE", "Phycis blennoides" = "GFB",
              "Thunnus alalunga" = "ALB", "Coryphaena hippurus" = "DOL",
              "Seriola dumerili" = "AMB", "Scomber colias" = "VMA")

theme_ms <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace% theme(
    text = element_text(colour = ms_ink[["primary"]]),
    axis.text = element_text(colour = ms_ink[["secondary"]], size = rel(0.85)),
    axis.title = element_text(colour = ms_ink[["secondary"]], size = rel(0.9)),
    panel.grid.major = element_line(colour = ms_ink[["grid"]], linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = rel(0.9),
                              margin = margin(2, 0, 4, 0)),
    legend.position = "top", legend.justification = "left",
    legend.title = element_text(colour = ms_ink[["secondary"]], size = rel(0.9)),
    plot.title = element_text(face = "bold", hjust = 0, size = rel(1.15),
                              margin = margin(0, 0, 4, 0)),
    plot.subtitle = element_text(colour = ms_ink[["secondary"]], hjust = 0,
                                 margin = margin(0, 0, 8, 0)),
    plot.title.position = "plot", panel.spacing = unit(10, "pt"))
}

sp_label <- function(sys) {
  setNames(sprintf("%s (%s)  r = %.2f, p = %g", sys$species, ms_codes[sys$species],
                   sys$r, sys$price), sys$species)
}

# reference lines: Bmsy/Fmsy = 1 (solid), Bpa and Blim (dashed, dotted)
ref_B <- function(Bpa = 0.5, Blim = 0.25) list(
  geom_hline(yintercept = 1, colour = ms_ink[["ref"]], linewidth = 0.3),
  geom_hline(yintercept = Bpa, colour = ms_ink[["ref"]], linewidth = 0.3, linetype = "dashed"),
  geom_hline(yintercept = Blim, colour = ms_ink[["ref"]], linewidth = 0.3, linetype = "dotted"))
ref_F <- function() geom_hline(yintercept = 1, colour = ms_ink[["ref"]], linewidth = 0.3)

# -----------------------------------------------------------------
# 1. Equilibrium against total effort
# -----------------------------------------------------------------

#' Tidy equilibrium scan: list of ms_summary() outputs (named by scenario)
#' each with column E (total effort)
eq_long <- function(scans, sys) {
  do.call(rbind, lapply(names(scans), function(nm) {
    s <- scans[[nm]]; BB <- attr(s, "BBmsy"); FF <- attr(s, "FFmsy"); CC <- attr(s, "catch")
    data.frame(scenario = nm, E = rep(s$E, sys$S),
               species = rep(sys$species, each = nrow(s)),
               BBmsy = as.vector(BB), FFmsy = as.vector(FF), yield = as.vector(CC))
  }))
}

# effort at which the first species falls below Bpa (the constrained limit)
eq_limits <- function(scans, Bpa = 0.5) {
  do.call(rbind, lapply(names(scans), function(nm) {
    s <- scans[[nm]]; ok <- which(s$n_below_Bpa == 0 & s$E > 0)
    i <- ok[which.max(s$yield[ok])]
    data.frame(scenario = nm, E = s$E[i], yield = s$yield[i], profit = s$profit[i])
  }))
}

#' Small multiples by species of B/Bmsy, F/Fmsy or yield against total effort
plot_eq_species <- function(scans, sys, what = c("BBmsy", "FFmsy", "yield"),
                            Bpa = 0.5, Blim = 0.25, ncol = 5) {
  what <- match.arg(what)
  d <- eq_long(scans, sys); d <- d[d$E > 0, ]
  d$species <- factor(d$species, sys$species, labels = ms_codes[sys$species])
  lim <- eq_limits(scans, Bpa)
  lim$scenario <- factor(lim$scenario, names(scans)); d$scenario <- factor(d$scenario, names(scans))
  ylab <- c(BBmsy = "B / Bmsy", FFmsy = "F / Fmsy", yield = "Equilibrium yield (t)")[[what]]
  p <- ggplot(d, aes(E, .data[[what]], colour = scenario)) +
    geom_vline(data = lim, aes(xintercept = E, colour = scenario), linewidth = 0.3,
               linetype = "longdash", show.legend = FALSE) +
    (if (what == "BBmsy") ref_B(Bpa, Blim) else if (what == "FFmsy") ref_F() else NULL) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ species, ncol = ncol, scales = if (what == "yield") "free_y" else "fixed") +
    scale_x_log10(labels = scales::label_number(big.mark = ",")) +
    scale_colour_manual(values = ms_cols, name = NULL) +
    labs(x = "Total effort (equal split across segments, log scale)", y = ylab) + theme_ms()
  if (what == "BBmsy") p <- p + coord_cartesian(ylim = c(0, 2))
  if (what == "FFmsy") p <- p + coord_cartesian(ylim = c(0, 4))
  p
}

#' Aggregate equilibrium: total yield, profit, and the worst-off stock
plot_eq_total <- function(scans, Bpa = 0.5, Blim = 0.25, E_max = 600) {
  d <- do.call(rbind, lapply(names(scans), function(nm) {
    s <- scans[[nm]]; s <- s[s$E > 0 & s$E <= E_max, ]
    rbind(data.frame(scenario = nm, E = s$E, panel = "Total yield (t)", value = s$yield),
          data.frame(scenario = nm, E = s$E, panel = "Profit", value = s$profit),
          data.frame(scenario = nm, E = s$E, panel = "Lowest B / Bmsy", value = s$min_BBmsy))
  }))
  d$panel <- factor(d$panel, c("Total yield (t)", "Profit", "Lowest B / Bmsy"))
  d$scenario <- factor(d$scenario, names(scans))
  lim <- eq_limits(scans, Bpa); lim$scenario <- factor(lim$scenario, names(scans))
  refs <- data.frame(panel = factor("Lowest B / Bmsy", levels(d$panel)),
                     y = c(1, Bpa, Blim), lt = c("solid", "dashed", "dotted"))
  ggplot(d, aes(E, value, colour = scenario)) +
    geom_hline(data = refs, aes(yintercept = y, linetype = lt), colour = ms_ink[["ref"]],
               linewidth = 0.3, show.legend = FALSE) +
    geom_hline(data = data.frame(panel = factor("Profit", levels(d$panel))), aes(yintercept = 0),
               colour = ms_ink[["ref"]], linewidth = 0.3) +
    geom_vline(data = lim, aes(xintercept = E, colour = scenario), linewidth = 0.3,
               linetype = "longdash", show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_linetype_identity() +
    scale_x_log10(labels = scales::label_number(big.mark = ",")) +
    scale_colour_manual(values = ms_cols, name = NULL) +
    labs(x = "Total effort (equal split across segments, log scale)", y = NULL) + theme_ms()
}

#' Effort by metier against total effort (stacked shares), one panel per scenario
plot_eq_effort <- function(sims, sys, E) {
  d <- do.call(rbind, lapply(names(sims), function(nm) {
    Em <- apply(sims[[nm]]$E_metier, c(2, 3), mean)
    sh <- Em / pmax(rowSums(Em), 1e-12)
    data.frame(scenario = nm, E = rep(E, sys$M),
               metier = rep(gsub("_", " ", sys$metiers$name), each = length(E)),
               share = as.vector(sh))
  }))
  d <- d[d$E > 0, ]
  d$metier <- factor(d$metier, gsub("_", " ", sys$metiers$name[c(1, 4, 2, 5, 3, 6)]))
  d$scenario <- factor(d$scenario, names(sims))
  ggplot(d, aes(E, share, fill = metier)) +
    geom_area(colour = "#fcfcfb", linewidth = 0.25, position = "stack") +
    facet_wrap(~ scenario) +
    scale_x_log10(labels = scales::label_number(big.mark = ","), expand = c(0, 0)) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    scale_fill_manual(values = ms_cols, name = "Métier") +
    labs(x = "Total effort (equal split across segments, log scale)",
         y = "Share of effort") + theme_ms()
}

#' Kobe plot at equilibrium: one point per species, panels = scenarios
plot_kobe_eq <- function(points, sys, Bpa = 0.5, Blim = 0.25) {
  points$code <- ms_codes[points$species]
  points$group <- sys$group[match(points$species, sys$species)]
  ggplot(points, aes(BBmsy, FFmsy)) +
    annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = 1, fill = "#0ca30c", alpha = 0.07) +
    annotate("rect", xmin = -Inf, xmax = 1, ymin = 1, ymax = Inf, fill = "#d03b3b", alpha = 0.07) +
    geom_vline(xintercept = 1, colour = ms_ink[["ref"]], linewidth = 0.3) +
    geom_vline(xintercept = Bpa, colour = ms_ink[["ref"]], linewidth = 0.3, linetype = "dashed") +
    geom_vline(xintercept = Blim, colour = ms_ink[["ref"]], linewidth = 0.3, linetype = "dotted") +
    geom_hline(yintercept = 1, colour = ms_ink[["ref"]], linewidth = 0.3) +
    geom_point(aes(colour = group), size = 2.6, stroke = 0) +
    ggrepel::geom_text_repel(aes(label = code), size = 2.7, colour = ms_ink[["primary"]],
                             min.segment.length = 0.2, segment.colour = ms_ink[["muted"]],
                             segment.size = 0.25, box.padding = 0.25, max.overlaps = Inf, seed = 1) +
    facet_wrap(~ scenario) +
    scale_colour_manual(values = ms_cols[1:2], name = NULL) +
    scale_x_continuous(limits = c(0, 2)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = "B / Bmsy", y = "F / Fmsy") + theme_ms()
}

# -----------------------------------------------------------------
# 2. Projections over time
# -----------------------------------------------------------------

q90 <- function(a, f) apply(a, f, stats::quantile, c(0.05, 0.5, 0.95), na.rm = TRUE)

#' Tidy projection quantiles by species for several runs (named list of ms_sim)
traj_long <- function(runs, sys) {
  do.call(rbind, lapply(names(runs), function(nm) {
    x <- runs[[nm]]; ny <- dim(x$B)[1]
    bb <- sweep(x$B, 3, sys$K / 2, "/"); ff <- sweep(x$F, 3, sys$r / 2, "/")
    one <- function(a, var) {
      q <- q90(a, c(1, 3))            # 3 x years x species
      data.frame(run = nm, var = var, year = rep(seq_len(ny), sys$S),
                 species = rep(sys$species, each = ny),
                 lo = as.vector(q[1, , ]), med = as.vector(q[2, , ]), hi = as.vector(q[3, , ]))
    }
    rbind(one(bb, "BBmsy"), one(ff, "FFmsy"), one(x$C, "yield"))
  }))
}

#' Small multiples by species: B/Bmsy, F/Fmsy or catch over time (median, 90% band)
plot_traj_species <- function(runs, sys, what = c("BBmsy", "FFmsy", "yield"),
                              Bpa = 0.5, Blim = 0.25, ncol = 5) {
  what <- match.arg(what)
  d <- traj_long(runs, sys); d <- d[d$var == what, ]
  d$species <- factor(d$species, sys$species, labels = ms_codes[sys$species])
  d$run <- factor(d$run, names(runs))
  ylab <- c(BBmsy = "B / Bmsy", FFmsy = "F / Fmsy", yield = "Catch (t)")[[what]]
  p <- ggplot(d, aes(year, med, colour = run, fill = run)) +
    (if (what == "BBmsy") ref_B(Bpa, Blim) else if (what == "FFmsy") ref_F() else NULL) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ species, ncol = ncol, scales = if (what == "yield") "free_y" else "fixed") +
    scale_colour_manual(values = ms_cols, name = NULL) +
    scale_fill_manual(values = ms_cols, name = NULL) +
    labs(x = "Year", y = ylab) + theme_ms()
  if (what == "BBmsy") p <- p + coord_cartesian(ylim = c(0, 2.2))
  if (what == "FFmsy") p <- p + coord_cartesian(ylim = c(0, 5))
  p
}

#' Effort by segment and total yield over time (median, 90% band)
plot_traj_fleet <- function(runs, sys) {
  d <- do.call(rbind, lapply(names(runs), function(nm) {
    x <- runs[[nm]]; ny <- dim(x$B)[1]
    qe <- q90(x$E_seg, c(1, 3)); qy <- q90(apply(x$C, 1:2, sum), 1)
    rbind(data.frame(run = nm, panel = paste("Effort:", rep(names(sys$segments), each = ny)),
                     year = rep(seq_len(ny), sys$G), lo = as.vector(qe[1, , ]),
                     med = as.vector(qe[2, , ]), hi = as.vector(qe[3, , ])),
          data.frame(run = nm, panel = "Total yield (t)", year = seq_len(ny),
                     lo = qy[1, ], med = qy[2, ], hi = qy[3, ]))
  }))
  d$panel <- factor(d$panel, c(paste("Effort:", names(sys$segments)), "Total yield (t)"))
  d$run <- factor(d$run, names(runs))
  ggplot(d, aes(year, med, colour = run, fill = run)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    scale_y_continuous(limits = c(0, NA)) +
    scale_colour_manual(values = ms_cols, name = NULL) +
    scale_fill_manual(values = ms_cols, name = NULL) +
    labs(x = "Year", y = NULL) + theme_ms()
}
