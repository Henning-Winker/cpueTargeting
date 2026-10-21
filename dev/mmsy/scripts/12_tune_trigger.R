# run from the repository root: tune the trigger of the effort rule (and of the
# Dentex reef closure) so that every stock is below Bpa at most 5% of the time
source("dev/mmsy/R/derive_landscape.R"); source("dev/mmsy/R/msfleet.R"); source("dev/mmsy/R/plots.R")
library(patchwork)
sys <- readRDS("dev/mmsy/output/tuned_sys.rds"); sq <- attr(sys, "status_quo")
args <- list(sys = sys, E_ref = sq$E_seg, B_init = sq$B, gamma = 1.5, proc_cv = 0.1, rho = 0.3, impl_cv = 0.1,
             capacity = 1000, risk = c(worst_P_below_Bpa = 0.05), objective = "yield")
tr <- seq(0.5, 1.5, by = 0.05)
t0 <- Sys.time()
A <- do.call(ms_tune_hcr, c(args, list(triggers = tr)))
B <- do.call(ms_tune_hcr, c(args, list(triggers = tr, close = "Dentex dentex")))
C <- do.call(ms_tune_hcr, c(args, list(triggers = seq(0.4, 1.0, by = 0.1), close = "Dentex dentex",
                                      close_triggers = c(0.75, 1, 1.5, 2, 3, 5, 1e6), nsim = 100)))
cat("runs:", format(round(Sys.time() - t0, 1)), "\n")
show <- function(o, nm) { cat("\n==", nm, "\n"); print(o$table[seq(1, nrow(o$table), by = 2), ], row.names = FALSE)
  cat("tuned:\n"); print(o$best, row.names = FALSE) }
show(A, "Effort rule only"); show(B, "Effort rule + reef closure (same trigger)")
cat("\n== Separate triggers, best feasible:\n"); print(C$best, row.names = FALSE)
cat("feasible share of grid:", mean(C$table$worst_P_below_Bpa <= 0.05), "\n")
saveRDS(list(A = A, B = B, C = C), "dev/mmsy/output/tune_trigger.rds")

# figure: risk, yield and profit against the trigger
d <- rbind(cbind(rule = "Effort rule only", A$table[, c("Btrig", "worst_P_below_Bpa", "yield", "profit")]),
           cbind(rule = "Effort rule + Dentex reef closure", B$table[, c("Btrig", "worst_P_below_Bpa", "yield", "profit")]))
dl <- rbind(data.frame(d[, c("rule", "Btrig")], panel = "Risk: worst stock P(B < Bpa)", value = d$worst_P_below_Bpa),
            data.frame(d[, c("rule", "Btrig")], panel = "Yield (t)", value = d$yield),
            data.frame(d[, c("rule", "Btrig")], panel = "Profit", value = d$profit))
dl$panel <- factor(dl$panel, unique(dl$panel)); dl$rule <- factor(dl$rule, unique(dl$rule))
best <- rbind(if (!is.null(A$best)) data.frame(rule = "Effort rule only", Btrig = A$best$Btrig),
              if (!is.null(B$best)) data.frame(rule = "Effort rule + Dentex reef closure", Btrig = B$best$Btrig))
best$rule <- factor(best$rule, levels(dl$rule))
p <- ggplot(dl, aes(Btrig, value, colour = rule)) +
  geom_hline(data = data.frame(panel = factor("Risk: worst stock P(B < Bpa)", levels(dl$panel)), y = 0.05),
             aes(yintercept = y), linetype = "dashed", colour = ms_ink[["ref"]], linewidth = 0.3) +
  geom_vline(data = best, aes(xintercept = Btrig, colour = rule), linetype = "longdash", linewidth = 0.3, show.legend = FALSE) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  facet_wrap(~ panel, scales = "free_y") + scale_colour_manual(values = ms_cols, name = NULL) +
  labs(x = "Trigger (B/Bmsy)", y = NULL, title = "Tuning the trigger",
       subtitle = paste0("Years 31-50, 200 simulations with common random numbers. Dashed: risk limit 5%; ",
                         "long dash: tuned trigger (lowest risk-compliant, highest yield)")) + theme_ms()
ggsave("dev/mmsy/figs/tune_trigger.png", p, width = 11, height = 3.9, dpi = 150, bg = "#fcfcfb")
# heatmap: separate triggers (closure trigger 1e6 = permanent full closure)
tc <- C$table; tc$ok <- tc$worst_P_below_Bpa <= 0.05
lev <- sort(unique(tc$Btrig_close)); lab <- ifelse(lev > 100, "Full closure", as.character(lev))
tc$close_f <- factor(tc$Btrig_close, lev, lab)
bc <- C$best; if (!is.null(bc)) bc$close_f <- factor(bc$Btrig_close, lev, lab)
h <- ggplot(tc, aes(factor(Btrig), close_f, fill = yield)) +
  geom_tile(aes(alpha = ok), colour = "#fcfcfb", linewidth = 0.6) +
  scale_alpha_manual(values = c(`FALSE` = 0.25, `TRUE` = 1), guide = "none") +
  geom_text(aes(label = round(yield)), size = 2.6,
            colour = ifelse(tc$ok & tc$yield > median(tc$yield), "#fcfcfb", ms_ink[["primary"]])) +
  (if (!is.null(C$best)) geom_tile(data = bc, fill = NA, colour = ms_ink[["primary"]], linewidth = 1.2)) +
  scale_fill_gradient(low = "#cde2fb", high = "#104281", name = "Yield (t)") +
  labs(x = "Effort-rule trigger (B/Bmsy)", y = "Reef-closure trigger (Dentex B/Bmsy)",
       title = "Separate triggers for effort and closure",
       subtitle = "Faded: misses the risk limit (worst stock P(B < Bpa) > 5%); outlined: tuned (highest yield)") + theme_ms() +
  theme(legend.position = "right", panel.grid.major = element_blank())
ggsave("dev/mmsy/figs/tune_trigger_2d.png", h, width = 7.5, height = 5, dpi = 150, bg = "#fcfcfb")
