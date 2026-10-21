# run from the repository root: equilibrium yield over effort under different conceptual treatments
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
library(patchwork)
sys <- readRDS("dev/mmsy/output/tuned_sys.rds"); sq <- attr(sys, "status_quo")
Esq <- sq$E_seg; Etot_sq <- sum(Esq)
mult <- c(0, 10^seq(-1.5, log10(8), length.out = 90))
caps <- outer(mult, Esq)                                     # status-quo composition of days
no_reef <- replace(rep(1, sys$M), sys$metiers$name == "Reef_shallow", 0)
eq <- function(nm, ...) {
  # caps = registered days; the days used are those of the vessels that break even
  s <- ms_summary(ms_simulate(sys, caps, 400, capacity = caps, ...)); BB <- attr(s, "BBmsy")
  data.frame(treatment = nm, E = s$E_total, yield = s$yield, min_BBmsy = apply(BB, 1, min),
             n_below_Bpa = rowSums(BB < 0.5))
}
# 1. shared q (mse.SEA): q_i = F_i / E at the status quo, F_i = (F/Fmsy) * r / 2
q <- sq$FFmsy * sys$r / 2 / Etot_sq
E <- mult[mult <= 4] * Etot_sq
Bk <- pmax(outer(E, q / sys$r, function(e, a) 1 - e * a), 0)  # B/K, n x S
Yq <- rowSums(outer(E, q * sys$K) * Bk)
d <- rbind(
  data.frame(treatment = "Shared q (mse.SEA, analytic)", E = E, yield = Yq, min_BBmsy = apply(2 * Bk, 1, min),
             n_below_Bpa = rowSums(2 * Bk < 0.5)),
  eq("Métiers, profit-blind (γ = 0)", gamma = 0),
  eq("Métiers, profit-seeking (γ = 1.5)", gamma = 1.5),
  eq("Profit-seeking + shallow reef closed", gamma = 1.5, open = no_reef))
d$treatment <- factor(d$treatment, unique(d$treatment))
msy_sum <- sum(sys$r * sys$K / 4)
pts <- do.call(rbind, lapply(split(d, d$treatment), function(x) {
  i <- which.max(x$yield); ok <- which(x$min_BBmsy >= 0.5 & x$E > 0); j <- ok[which.max(x$yield[ok])]
  rbind(data.frame(x[i, ], point = "MMSY"), data.frame(x[j, ], point = "Safeguarded (all ≥ Bpa)")) }))
print(pts[, c("treatment", "point", "E", "yield", "min_BBmsy", "n_below_Bpa")], row.names = FALSE, digits = 3)
cat("sum of single-species MSY:", round(msy_sum), "| status quo: days", round(Etot_sq), "yield", round(sq$summary$yield), "\n")
saveRDS(list(curves = d, points = pts), "dev/mmsy/output/yield_curves.rds")

sqp <- data.frame(E = Etot_sq, yield = sq$summary$yield)
p1 <- ggplot(d, aes(E, yield, colour = treatment)) +
  geom_hline(yintercept = msy_sum, colour = ms_ink[["ref"]], linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = max(d$E), y = msy_sum, label = "Sum of single-species MSY", hjust = 1, vjust = -0.5,
           size = 3, colour = ms_ink[["secondary"]]) +
  geom_line(linewidth = 0.8) +
  geom_point(data = pts, aes(shape = point), size = 2.8, stroke = 1) +
  geom_point(data = sqp, aes(E, yield), inherit.aes = FALSE, shape = 23, size = 3.2, fill = ms_ink[["primary"]]) +
  annotate("text", x = Etot_sq, y = sq$summary$yield, label = "  status quo", hjust = 0, vjust = 1.6, size = 3,
           colour = ms_ink[["primary"]]) +
  scale_shape_manual(values = c(MMSY = 16, `Safeguarded (all ≥ Bpa)` = 1), name = NULL) +
  scale_colour_manual(values = ms_cols, name = NULL) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Total days fished (status-quo segment mix of caps)", y = "Equilibrium yield (t)", title = "A. Multispecies yield curves") +
  theme_ms() + theme(legend.position = "right")
p2 <- ggplot(d[d$E > 0, ], aes(E, min_BBmsy, colour = treatment)) +
  geom_hline(yintercept = c(1, 0.5, 0.25), linetype = c("solid", "dashed", "dotted"), colour = ms_ink[["ref"]], linewidth = 0.3) +
  geom_line(linewidth = 0.8) + scale_colour_manual(values = ms_cols, guide = "none") +
  labs(x = "Total days fished", y = "Lowest B/Bmsy", title = "B. The weakest stock") + theme_ms()
p3 <- ggplot(d[d$E > 0, ], aes(E, n_below_Bpa, colour = treatment)) +
  geom_step(linewidth = 0.8) + scale_colour_manual(values = ms_cols, guide = "none") +
  scale_y_continuous(breaks = 0:10) +
  labs(x = "Total days fished", y = "Stocks below Bpa", title = "C. Stocks below Bpa") + theme_ms()
fig <- p1 / (p2 | p3) + plot_layout(heights = c(1.3, 1)) +
  plot_annotation(title = "Equilibrium yield over effort: conceptual treatments compared",
    subtitle = paste0("Tuned system, deterministic equilibrium. Caps scaled from the status-quo days; the days fished are those of the vessels\n",
                      "that break even, so the simulated curves end where effort is limited by profit. Shared q is calibrated to the status quo."),
    theme = theme_ms())
ggsave("dev/mmsy/figs/yield_curves.png", fig, width = 11, height = 8, dpi = 150, bg = "#fcfcfb")
