#' Run scenarios over several seeds (replicates)
#'
#' Simulates each scenario for each seed, standardises the data and keeps the
#' indices and accuracy diagnostics. Results are saved to `file` after every
#' fit, and seeds/scenarios already in `file` are skipped, so a long (e.g.
#' overnight) run can be interrupted and resumed.
#'
#' @param scenarios Named list of settings lists (see [tuna_settings()]).
#' @param seeds Integer vector of seeds.
#' @param models Models passed to [standardize_tuna_indices()].
#' @param extra Optional list of extra standardisation calls, each a list of
#'   arguments for [standardize_tuna_indices()] (e.g. a 1 x 1 control with a
#'   `label`); their results are added with [bind_comparisons()].
#' @param file Optional `.rds` file to save and resume from.
#' @param ... Further arguments for [standardize_tuna_indices()] (e.g.
#'   `spatial_res`, `n_pc`, `engine`).
#' @return A list with `diagnostics` and `indices`, each with columns
#'   `scenario` and `seed` added.
#' @export
#' @examples
#' \dontrun{
#' base <- tuna_settings()
#' scenarios <- list(Fixed = base,
#'                   Moving = modifyList(base, list(thermocline_annual_sd = 15,
#'                                                  thermocline_annual_field_sd = 10)))
#' reps <- run_replicates(scenarios, seeds = 1:10,
#'                        models = c("Base+S", "Cluster+S", "DPC+S", "HBF+S"),
#'                        engine = "bam", file = "replicates.rds")
#' plot_replicates(reps)
#' }
run_replicates <- function(scenarios, seeds, models, extra = NULL, file = NULL, ...) {
  if (is.null(names(scenarios)) || any(names(scenarios) == "")) {
    stop("`scenarios` must be a named list of settings.", call. = FALSE)
  }
  out <- list(diagnostics = NULL, indices = NULL)
  if (!is.null(file) && file.exists(file)) out <- readRDS(file)
  done <- if (is.null(out$diagnostics)) character(0) else
    unique(paste(out$diagnostics$scenario, out$diagnostics$seed))

  for (seed in seeds) {
    for (sc in names(scenarios)) {
      if (paste(sc, seed) %in% done) next
      message("Scenario ", sc, ", seed ", seed)
      sim <- simulate_tuna(scenarios[[sc]], seed = seed, plots = FALSE)
      fit <- standardize_tuna_indices(sim, models = models, ...)
      for (args in extra) {
        fit <- bind_comparisons(fit, do.call(standardize_tuna_indices,
                                             c(list(simulation = sim), args)))
      }
      tag <- function(x) data.frame(scenario = sc, seed = seed, x)
      out$diagnostics <- rbind(out$diagnostics, tag(fit$diagnostics))
      out$indices <- rbind(out$indices, tag(fit$indices))
      if (!is.null(file)) saveRDS(out, file)
    }
  }
  model_levels <- unique(as.character(out$diagnostics$model))
  out$diagnostics$model <- factor(out$diagnostics$model, levels = model_levels)
  out$indices$model <- factor(out$indices$model, levels = model_levels)
  out$diagnostics$scenario <- factor(out$diagnostics$scenario, levels = names(scenarios))
  out$indices$scenario <- factor(out$indices$scenario, levels = names(scenarios))
  out
}

#' Accuracy of standardised indices across replicates
#'
#' Boxplots of the RMSE of the log index over seeds, by model, scenario and
#' species.
#'
#' @param replicates Output of [run_replicates()].
#' @param models Optional subset of model names to show.
#' @return A ggplot object.
#' @export
plot_replicates <- function(replicates, models = NULL) {
  diag <- replicates$diagnostics
  if (!is.null(models)) diag <- diag[diag$model %in% models, ]
  n_seeds <- length(unique(diag$seed))
  ggplot(diag, aes(x = .data$model, y = .data$rmse_log, fill = .data$scenario)) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(~species) +
    labs(x = NULL, y = "RMSE of log index", fill = "Scenario",
         title = paste0("Accuracy across ", n_seeds, " replicates")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
