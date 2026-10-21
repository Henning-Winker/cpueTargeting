# run from the repository root after 06_optimise.R and 07_viability_closure.R
source("dev/mmsy/R/plots.R"); library(patchwork)
a <- readRDS("dev/mmsy/output/optimise.rds"); b <- readRDS("dev/mmsy/output/viability_closure.rds")
opts <- data.frame(
  option = c("One global cap (equal split)", "Segment caps, each ≥ 20%",
             "Demersal + pelagic caps, no mixed", "Three caps, unrestricted",
             "Equal split + shallow reef closed", "Segment caps ≥ 20% + reef closed"),
  yield = c(b[[1]]$summary$yield, b[[2]]$summary$yield,
            a$res1[["Demersal + pelagic caps, mixed = 0 | yield"]]$summary$yield,
            a$res1[["Three segment caps | yield"]]$summary$yield, b[[4]]$summary$yield, b[[3]]$summary$yield),
  binding = c("DEC", "DEC", "DEC", "DEC", "AMB", "AMB"),
  type = c("Effort caps", "Effort caps", "Effort caps (demersal fleet closed)", "Effort caps (demersal fleet closed)",
           "Caps + métier closure", "Caps + métier closure"))
opts$option <- factor(opts$option, rev(opts$option))
opts$type <- factor(opts$type, unique(opts$type))
p1 <- ggplot(opts, aes(yield, option, colour = type)) +
  geom_segment(aes(x = 0, xend = yield, yend = option), linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.0f  (binding: %s)", yield, binding)), hjust = -0.15, size = 3,
            colour = ms_ink[["primary"]]) +
  scale_x_continuous(limits = c(0, 1900), expand = c(0, 0)) +
  scale_colour_manual(values = ms_cols[c(1, 2, 3)], name = NULL) +
  labs(x = "Maximum equilibrium yield with every stock ≥ Bpa (t)", y = NULL,
       title = "A. Safeguarded MMSY by management option") +
  theme_ms() + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
r2 <- a$res2; lim <- min(r2$mixed_cap[!r2$feasible])
p2 <- ggplot(r2[r2$feasible, ], aes(mixed_cap, yield)) +
  annotate("rect", xmin = lim, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#d03b3b", alpha = 0.07) +
  annotate("text", x = lim + 1, y = 820, label = "Infeasible:\nno caps keep\nall stocks ≥ Bpa", hjust = 0,
           vjust = 1, size = 3, colour = ms_ink[["secondary"]]) +
  geom_line(colour = ms_cols[1], linewidth = 0.7) + geom_point(colour = ms_cols[1], size = 2) +
  scale_y_continuous(limits = c(0, NA)) + scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
  labs(x = "Uncontrolled mixed-segment effort (free agent)", y = "Maximum safeguarded yield (t)",
       title = "B. Mixed segment as a free agent",
       subtitle = "Demersal and pelagic caps optimised for each level.\nAbout 87% of mixed effort goes to demersal métiers.") +
  theme_ms()
ggsave("dev/mmsy/figs/opt_options.png", p1 + p2 + plot_layout(widths = c(1.25, 1)),
       width = 12, height = 4.6, dpi = 150, bg = "#fcfcfb")
