#' Select operations likely to have fished the habitat of a target species
#'
#' Stephens and MacCall (2004): a binomial GLM predicts the presence of the
#' target species from the presence or absence of the other species in each
#' operation. Operations with a predicted probability above a critical value
#' are kept. The critical value is the one at which the number of predicted
#' presences equals the number of observed presences, i.e. it minimises
#' \eqn{|\sum obs - \sum (\hat p > p_{crit})|}, which balances false
#' positives and false negatives. This is the implementation used in the
#' simulation study of Winker et al. (2014).
#'
#' The species set matters: species that co-occur with the target (positive
#' coefficients) and species of other habitats (negative coefficients)
#' both carry information. Use [select_species()] first, or give `species`.
#'
#' @param x A `cpue_data` object.
#' @param target Target species (column or original name).
#' @param species Predictor species (default: all others in `x`).
#' @param thresholds Candidate critical probabilities.
#' @param min_keep Minimum share of operations to keep; if the critical value
#'   would keep fewer, the largest threshold that keeps at least this share
#'   is used and a message is given.
#' @return An object of class `stephens_maccall`: a list with `fit` (the
#'   GLM), `coefficients` (species, estimate, SE), `prob` (predicted
#'   probability for each operation), `curve` (threshold, difference between
#'   observed and predicted presences, and number of operations kept),
#'   `critical` (the critical probability), `selected` (logical, one per
#'   operation), `target`, and `year` (for the plots).
#' @references Stephens, A. and MacCall, A. (2004) A multispecies approach to
#'   subsetting logbook data for purposes of estimating CPUE. Fisheries
#'   Research 70: 299-310.
#' @seealso [plot_stephens_maccall()], [subset.cpue_data()]
#' @export
#' @examples
#' \dontrun{
#' sm <- stephens_maccall(x, target = "Pagellus bogaraveo")
#' plot_stephens_maccall(sm)
#' x_sm <- subset(x, sm$selected, note = "Stephens & MacCall")
#' }
stephens_maccall <- function(x, target, species = NULL,
                             thresholds = seq(0.01, 0.99, by = 0.01), min_keep = 0) {
  check_cpue_data(x)
  target <- species_codes(x, target)
  species <- if (is.null(species)) setdiff(x$species, target) else species_codes(x, species)
  species <- setdiff(species, target)
  if (!length(species)) stop("Need at least one predictor species.", call. = FALSE)
  pa <- as.data.frame((as.matrix(x$data[, species, drop = FALSE]) > 0) * 1)
  # species present everywhere or nowhere carry no information
  informative <- vapply(pa, function(z) length(unique(z)) > 1, logical(1))
  if (!all(informative)) {
    message("Dropped species without contrast in presence: ",
            paste(sp_label(x, species[!informative]), collapse = ", "))
    pa <- pa[, informative, drop = FALSE]
    species <- species[informative]
  }
  obs <- as.numeric(x$data[[target]] > 0)
  if (length(unique(obs)) < 2) {
    stop("The target species is present in all (or no) operations; ",
         "subsetting by habitat is not needed.", call. = FALSE)
  }
  form <- stats::as.formula(paste("obs ~", paste(species, collapse = " + ")))
  fit <- stats::glm(form, data = data.frame(obs = obs, pa), family = stats::binomial())
  prob <- as.numeric(stats::fitted(fit))
  kept <- vapply(thresholds, function(t) sum(prob > t), numeric(1))
  diff <- abs(sum(obs) - kept)
  crit <- thresholds[which.min(diff)]
  if (mean(prob > crit) < min_keep) {
    ok <- thresholds[kept / length(prob) >= min_keep]
    if (!length(ok)) stop("No threshold keeps `min_keep` of the operations.", call. = FALSE)
    message(sprintf("Critical value %.2f keeps only %.0f%% of operations; using %.2f.",
                    crit, 100 * mean(prob > crit), max(ok)))
    crit <- max(ok)
  }
  cf <- summary(fit)$coefficients
  cf <- cf[rownames(cf) != "(Intercept)", , drop = FALSE]
  coefs <- data.frame(species = sp_label(x, rownames(cf)), estimate = cf[, 1],
                      se = cf[, 2], row.names = NULL)
  structure(list(fit = fit, coefficients = coefs, prob = prob,
                 curve = data.frame(threshold = thresholds, difference = diff, kept = kept),
                 critical = crit, selected = prob > crit, observed = obs,
                 target = sp_label(x, target), year = x$data[[x$roles$year]]),
            class = "stephens_maccall")
}

#' @export
print.stephens_maccall <- function(x, ...) {
  cat("<stephens_maccall> target:", x$target, "\n")
  cat(sprintf("  critical probability: %.2f\n", x$critical))
  cat(sprintf("  operations kept: %d of %d (%.0f%%)\n", sum(x$selected), length(x$selected),
              100 * mean(x$selected)))
  cat(sprintf("  target present in %.0f%% of kept and %.0f%% of dropped operations\n",
              100 * mean(x$observed[x$selected]),
              100 * mean(x$observed[!x$selected])))
  invisible(x)
}

#' Diagnostic plots for Stephens and MacCall subsetting
#'
#' @param sm Output of [stephens_maccall()].
#' @return A list of ggplots: `threshold` (difference between observed and
#'   predicted presences against the threshold, with the critical value),
#'   `coefficients` (species effects on the presence of the target with
#'   95% CIs), `probability` (distribution of predicted probabilities for
#'   operations with and without the target) and `kept` (share of operations
#'   kept by year; a trend means the subset changes over time).
#' @export
plot_stephens_maccall <- function(sm) {
  if (!inherits(sm, "stephens_maccall")) stop("`sm` must come from stephens_maccall().",
                                              call. = FALSE)
  cf <- sm$coefficients
  cf$species <- factor(cf$species, levels = cf$species[order(cf$estimate)])
  pr <- data.frame(prob = sm$prob,
                   target = factor(ifelse(sm$observed == 1, "present", "absent"),
                                   levels = c("absent", "present")))
  kept <- stats::aggregate(list(kept = sm$selected), by = list(year = sm$year), FUN = mean)
  list(
    threshold = ggplot(sm$curve, aes(x = .data$threshold, y = .data$difference)) +
      geom_line() +
      geom_vline(xintercept = sm$critical, linetype = 2, colour = "red") +
      labs(x = "Probability threshold",
           y = "|observed - predicted presences|",
           title = paste("Stephens & MacCall:", sm$target),
           subtitle = sprintf("critical value = %.2f", sm$critical)) + theme_bw(),
    coefficients = ggplot(cf, aes(x = .data$estimate, y = .data$species)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_errorbar(aes(xmin = .data$estimate - 1.96 * .data$se,
                        xmax = .data$estimate + 1.96 * .data$se),
                    width = 0.2, orientation = "y") +
      geom_point(size = 2) +
      labs(x = "Effect on presence of target (logit)", y = NULL) + theme_bw(),
    probability = ggplot(pr, aes(x = .data$prob, fill = .data$target)) +
      geom_histogram(bins = 40, position = "identity", alpha = 0.6) +
      geom_vline(xintercept = sm$critical, linetype = 2, colour = "red") +
      labs(x = "Predicted probability of target", y = "Operations", fill = "Target") +
      theme_bw(),
    kept = ggplot(kept, aes(x = .data$year, y = .data$kept)) +
      geom_line() + geom_point() + coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Year", y = "Share of operations kept") + theme_bw()
  )
}
