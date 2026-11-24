# Consistency checks of species recording, vessel set and catch composition
# over the time series.

#' Check the consistency of species recording over time
#'
#' Long logbook and sales-note series often change their species coding:
#' codes are introduced, split, merged or retired when databases or
#' identification practices change. Because targeting covariates
#' ([add_targeting()]) are computed from the catch composition, such coding
#' breaks look like changes in targeting and can bias the standardised index
#' (Winker et al. 2013 discuss the need for consistent species groups).
#'
#' `species_consistency()` computes the share of operations recording each
#' species in each year and looks for abrupt changes:
#' * for every species and every pair of consecutive years it compares the
#'   mean occurrence in the `window` years before and after; a **step** is a
#'   ratio of at least `ratio` (with the larger side recorded in at least
#'   `min_occurrence` of operations). A species can have several steps;
#' * steps are classified as `appears`, `disappears` (mean occurrence below
#'   `absent` on one side), `increase` or `decrease`;
#' * years with steps in at least `min_species` species are candidate
#'   **coding breaks**; the years between breaks form **blocks** within which
#'   the recording is consistent;
#' * at each break, species that disappear are paired with species that
#'   appear or increase at the same time (negative correlation of their
#'   occurrence series around the break) as **candidate replacements**,
#'   ranked by the similarity of their occurrence level and of the vessels
#'   that record them before and after the break (a replacement code is
#'   usually recorded by the same vessels at a similar rate).
#'   The best candidate for each disappearing species is returned in
#'   `suggested`, ready for [harmonise_species()] or the `synonyms` argument
#'   of [cpue_data()].
#'
#' The suggestions are inferred only from the timing and must be confirmed
#' with the species code lists of the data provider. Run the check on the
#' full data, before [select_species()] removes species.
#'
#' @param x A `cpue_data` object.
#' @param min_occurrence Minimum yearly share of operations for a species to
#'   be checked, and for the larger side of a step.
#' @param ratio Minimum ratio of the mean occurrence after/before (or
#'   before/after) a split for a step.
#' @param window Number of years on each side of a split.
#' @param min_species Minimum number of species with a step in the same year
#'   (or in adjacent years) for a coding break.
#' @param absent Mean occurrence below which a species is treated as absent.
#' @param min_score Minimum `score` of a candidate pair for a suggestion.
#' @return An object of class `species_consistency`, a list with
#'   * `occurrence`: matrix (species x year) of the share of operations
#'     recording each species;
#'   * `steps`: one row per step (`species`, `label`, `year` = first year
#'     after the step, `type`, `before`, `after`, `log_ratio`);
#'   * `years`: per year, the number of species recorded, steps and the
#'     Bray-Curtis dissimilarity of the occurrence profile to the previous year;
#'   * `breaks`: candidate breaks with their `type` (species replaced, i.e. a
#'     coding change; species lost, e.g. an extraction for some species only;
#'     species added) and the species involved;
#'   * `blocks`: consistent blocks of years;
#'   * `candidates`: candidate replacement pairs (`from`, `to`, `break_year`,
#'     `correlation` of the occurrence series around the break,
#'     `level_similarity` of the loss of `from` and the gain of `to`, 0-1,
#'     `vessel_similarity` of the vessels recording them, cosine 0-1, and
#'     their mean `score`);
#'   * `suggested`: named character vector `c(from = "to")` of the best
#'     replacement for each disappearing species (original names).
#' @references Winker, H., Kerwath, S.E. and Attwood, C.G. (2013) Comparison
#'   of two approaches to standardize catch-per-unit-effort for targeting
#'   behaviour in a multispecies hand-line fishery. Fisheries Research 139,
#'   118-131.
#' @seealso [plot_species_consistency()], [harmonise_species()],
#'   [data_consistency()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' sc <- species_consistency(x)
#' sc
species_consistency <- function(x, min_occurrence = 0.02, ratio = 4, window = 3,
                                min_species = 2, absent = 0.002, min_score = 0.5) {
  check_cpue_data(x)
  yr <- x$data[[x$roles$year]]
  years <- sort(unique(yr))
  ny <- length(years)
  if (ny < 2 * window) stop("Need at least ", 2 * window, " years.", call. = FALSE)
  pos <- as.matrix(x$data[, x$species, drop = FALSE]) > 0
  occ <- t(rowsum(pos * 1, factor(yr, levels = years))) /
    rep(as.numeric(table(factor(yr, levels = years))), each = ncol(pos))
  dimnames(occ) <- list(x$species, years)
  checked <- rownames(occ)[apply(occ, 1, max) >= min_occurrence]
  eps <- min_occurrence / 2
  thr <- log(ratio)

  steps <- list()
  for (s in checked) {
    o <- occ[s, ]
    lr <- before <- after <- rep(NA_real_, ny - 1)
    for (k in seq_len(ny - 1)) {
      before[k] <- mean(o[max(1, k - window + 1):k])
      after[k] <- mean(o[(k + 1):min(ny, k + window)])
      lr[k] <- log((after[k] + eps) / (before[k] + eps))
    }
    cand <- which(abs(lr) >= thr & pmax(before, after) >= min_occurrence)
    # keep local maxima of |log ratio| (one step per break)
    peak <- vapply(cand, function(k) {
      nb <- max(1, k - window):min(ny - 1, k + window)
      abs(lr[k]) >= max(abs(lr[nb]))
    }, logical(1))
    for (k in cand[peak]) {
      type <- if (before[k] < absent) "appears" else if (after[k] < absent) "disappears" else
        if (lr[k] > 0) "increase" else "decrease"
      steps[[length(steps) + 1]] <- data.frame(
        species = s, label = sp_label(x, s), year = years[k + 1], type = type,
        before = before[k], after = after[k], log_ratio = lr[k], stringsAsFactors = FALSE)
    }
  }
  steps <- if (length(steps)) do.call(rbind, steps) else
    data.frame(species = character(0), label = character(0), year = numeric(0),
               type = character(0), before = numeric(0), after = numeric(0),
               log_ratio = numeric(0))
  steps <- steps[order(steps$year, steps$type, -abs(steps$log_ratio)), ]
  rownames(steps) <- NULL

  # Year table
  oc <- occ[checked, , drop = FALSE]
  bc <- c(NA, vapply(2:ny, function(j) {
    a <- oc[, j - 1]; b <- oc[, j]
    if (sum(a + b) == 0) return(0)
    sum(abs(a - b)) / sum(a + b)
  }, numeric(1)))
  n_steps <- as.numeric(table(factor(steps$year, levels = years)))
  years_tab <- data.frame(year = years, operations = as.numeric(table(factor(yr, levels = years))),
                          species_recorded = colSums(occ > 0), steps = n_steps,
                          dissimilarity = bc, row.names = NULL)

  # Coding breaks: years (merging adjacent years) with steps in >= min_species species
  breaks <- data.frame(year = numeric(0), n_species = integer(0), species = character(0),
                       stringsAsFactors = FALSE)
  if (nrow(steps)) {
    by_year <- tapply(steps$species, steps$year, function(z) length(unique(z)))
    cand_years <- as.numeric(names(by_year))
    groups <- cumsum(c(1, diff(cand_years) > 1))
    for (g in unique(groups)) {
      yy <- cand_years[groups == g]
      st <- steps[steps$year %in% yy, ]
      if (length(unique(st$species)) >= min_species) {
        main <- yy[which.max(by_year[as.character(yy)])]
        breaks <- rbind(breaks, data.frame(
          year = main, n_species = length(unique(st$species)),
          species = paste(unique(st$label), collapse = ", "), stringsAsFactors = FALSE))
      }
    }
  }
  # Blocks between breaks
  cuts <- sort(unique(breaks$year))
  starts <- c(min(years), cuts)
  ends <- c(cuts - 1, max(years))
  blocks <- data.frame(block = seq_along(starts), from = starts, to = ends)
  blocks <- blocks[blocks$from <= blocks$to, ]

  # Candidate replacements at each break: losers paired with gainers whose
  # occurrence series is negatively correlated around the break; ranked by
  # the similarity of the occurrence level (before vs after) and, with a
  # vessel column, of the vessel profile (which vessels record them).
  vcol <- x$roles$vessel
  vprof <- function(s, yy) {
    if (is.null(vcol)) return(NULL)
    k <- pos[, s] & yr %in% yy
    tab <- table(factor(x$data[[vcol]][k], levels = unique(x$data[[vcol]])))
    if (sum(tab) == 0) return(NULL)
    as.numeric(tab) / sum(tab)
  }
  candidates <- data.frame(from = character(0), to = character(0), break_year = numeric(0),
                           correlation = numeric(0), level_similarity = numeric(0),
                           vessel_similarity = numeric(0), score = numeric(0),
                           from_type = character(0), to_type = character(0),
                           stringsAsFactors = FALSE)
  for (b in breaks$year) {
    st <- steps[abs(steps$year - b) <= 1, ]
    losers <- unique(st$species[st$type %in% c("disappears", "decrease")])
    gainers <- unique(st$species[st$type %in% c("appears", "increase")])
    yy <- which(years >= b - 2 * window & years <= b + 2 * window - 1)
    pre <- years[years < b & years >= b - window]
    post <- years[years >= b & years < b + window]
    for (l in losers) for (g in gainers) {
      r <- suppressWarnings(stats::cor(occ[l, yy], occ[g, yy]))
      if (is.na(r) || r >= -0.5) next
      loss <- mean(occ[l, as.character(pre)]) - mean(occ[l, as.character(post)])
      gain <- mean(occ[g, as.character(post)]) - mean(occ[g, as.character(pre)])
      lev <- exp(-abs(log((max(loss, 0) + eps) / (max(gain, 0) + eps))))
      pl <- vprof(l, pre); pg <- vprof(g, post)
      vs <- if (is.null(pl) || is.null(pg)) NA_real_ else
        sum(pl * pg) / sqrt(sum(pl^2) * sum(pg^2))
      candidates <- rbind(candidates, data.frame(
        from = sp_label(x, l), to = sp_label(x, g), break_year = b, correlation = r,
        level_similarity = lev, vessel_similarity = vs,
        score = mean(c(lev, vs), na.rm = TRUE),
        from_type = st$type[st$species == l][1], to_type = st$type[st$species == g][1],
        stringsAsFactors = FALSE))
    }
  }
  candidates <- candidates[order(candidates$break_year, candidates$from, -candidates$score), ]
  rownames(candidates) <- NULL
  sug <- candidates[candidates$from_type == "disappears" & candidates$score >= min_score, ]
  sug <- sug[!duplicated(sug$from), ]
  suggested <- stats::setNames(sug$to, sug$from)
  # Type of each break
  if (nrow(breaks)) {
    breaks$type <- vapply(breaks$year, function(b) {
      st <- steps[abs(steps$year - b) <= 1, ]
      n_loss <- sum(st$type %in% c("disappears", "decrease"))
      n_gain <- sum(st$type %in% c("appears", "increase"))
      if (any(sug$break_year == b)) "species replaced (coding change)" else
        if (n_loss > n_gain) "species lost (extraction or recording change)" else
          if (n_gain > n_loss) "species added" else "mixed"
    }, character(1))
    breaks <- breaks[, c("year", "type", "n_species", "species")]
  }

  structure(list(occurrence = occ, labels = stats::setNames(sp_label(x, rownames(occ)), rownames(occ)),
                 steps = steps, years = years_tab, breaks = breaks,
                 blocks = blocks, candidates = candidates, suggested = suggested,
                 settings = list(min_occurrence = min_occurrence, ratio = ratio,
                                 window = window, min_species = min_species,
                                 absent = absent, checked = checked)),
            class = "species_consistency")
}

#' @export
print.species_consistency <- function(x, ...) {
  cat("<species_consistency>", length(x$settings$checked), "species checked over",
      ncol(x$occurrence), "years;", nrow(x$steps), "steps\n")
  if (nrow(x$breaks)) {
    cat("Candidate coding breaks:\n")
    for (i in seq_len(nrow(x$breaks))) {
      cat(sprintf("  %g: %s, %d species (%s)\n", x$breaks$year[i], x$breaks$type[i],
                  x$breaks$n_species[i], x$breaks$species[i]))
    }
    cat("Consistent blocks:", paste0(x$blocks$from, "-", x$blocks$to, collapse = ", "), "\n")
  } else {
    cat("No coding breaks found.\n")
  }
  if (length(x$suggested)) {
    cat("Suggested harmonisation (confirm with the code lists):\n")
    cat(paste0("  ", names(x$suggested), " -> ", x$suggested, collapse = "\n"), "\n")
  }
  invisible(x)
}

#' Plot the species-recording consistency check
#'
#' @param x A [species_consistency()] result.
#' @param type `"heatmap"`: occurrence by year of the species with steps
#'   (dashed lines: coding breaks); `"years"`: dissimilarity of the
#'   occurrence profile to the previous year and number of steps per year;
#'   `"pairs"`: occurrence series of the candidate replacement pairs.
#' @param species Optional species (column or original names) to add to the
#'   heatmap, e.g. the target.
#' @param max_species Show at most this many species in the heatmap.
#' @return A ggplot.
#' @export
plot_species_consistency <- function(x, type = c("heatmap", "years", "pairs"),
                                     species = NULL, max_species = 40) {
  if (!inherits(x, "species_consistency")) stop("`x` must come from species_consistency().", call. = FALSE)
  type <- match.arg(type)
  occ <- x$occurrence
  years <- as.numeric(colnames(occ))
  vlines <- x$breaks$year - 0.5
  if (type == "heatmap") {
    sp <- unique(x$steps$species)
    if (!is.null(species)) {
      extra <- ifelse(species %in% rownames(occ), species,
                      names(x$labels)[match(species, x$labels)])
      sp <- unique(c(sp, stats::na.omit(extra)))
    }
    if (length(sp) == 0) sp <- rownames(occ)[seq_len(min(max_species, nrow(occ)))]
    if (length(sp) > max_species) sp <- sp[seq_len(max_species)]
    first_step <- tapply(x$steps$year, x$steps$species, min)[sp]
    first_step[is.na(first_step)] <- Inf
    sp <- sp[order(first_step, -rowMeans(occ[sp, , drop = FALSE]))]
    lab <- unname(x$labels[sp])
    d <- data.frame(species = factor(rep(lab, times = length(years)), levels = rev(unique(lab))),
                    year = rep(years, each = length(sp)),
                    occurrence = as.vector(occ[sp, , drop = FALSE]))
    return(ggplot(d, aes(x = .data$year, y = .data$species, fill = .data$occurrence)) +
             geom_tile() + scale_fill_viridis_c(name = "Share of\noperations") +
             geom_vline(xintercept = vlines, colour = "red", linetype = 2) +
             labs(x = "Year", y = "Species") + theme_bw())
  }
  if (type == "years") {
    y <- x$years
    d <- rbind(data.frame(year = y$year, value = y$dissimilarity,
                          what = "Dissimilarity to previous year"),
               data.frame(year = y$year, value = y$steps, what = "Species with a step"))
    return(ggplot(d[!is.na(d$value), ], aes(x = .data$year, y = .data$value)) +
             geom_col(data = function(z) z[z$what == "Species with a step", ], fill = "grey60") +
             geom_line(data = function(z) z[z$what != "Species with a step", ]) +
             geom_point(data = function(z) z[z$what != "Species with a step", ]) +
             geom_vline(xintercept = vlines, colour = "red", linetype = 2) +
             facet_wrap(~what, ncol = 1, scales = "free_y") +
             labs(x = "Year", y = NULL) + theme_bw())
  }
  cand <- x$candidates
  if (!nrow(cand)) stop("No candidate replacement pairs.", call. = FALSE)
  lab2code <- stats::setNames(names(x$labels), x$labels)
  rows <- lapply(seq_len(nrow(cand)), function(i) {
    f <- lab2code[[cand$from[i]]]; g <- lab2code[[cand$to[i]]]
    rbind(data.frame(pair = sprintf("%s -> %s (%g)", cand$from[i], cand$to[i], cand$break_year[i]),
                     year = years, occurrence = occ[f, ], role = "from"),
          data.frame(pair = sprintf("%s -> %s (%g)", cand$from[i], cand$to[i], cand$break_year[i]),
                     year = years, occurrence = occ[g, ], role = "to"))
  })
  d <- do.call(rbind, rows)
  ggplot(d, aes(x = .data$year, y = .data$occurrence, colour = .data$role)) +
    geom_line() + geom_point(size = 0.8) +
    geom_vline(xintercept = vlines, colour = "grey50", linetype = 2) +
    facet_wrap(~pair) + scale_colour_manual(values = c(from = "#d95f02", to = "#1b9e77")) +
    labs(x = "Year", y = "Share of operations", colour = NULL) + theme_bw()
}

#' Merge species codes in a cpue_data object
#'
#' Applies a harmonisation of species codes, e.g. the `suggested` element of
#' [species_consistency()] after confirmation with the code lists, to an
#' existing `cpue_data` object: catches of each recorded species are added to
#' the species they should be counted as, and the step is recorded in the
#' processing log. Equivalent to the `synonyms` argument of [cpue_data()],
#' but for an object that already exists.
#'
#' @param x A `cpue_data` object.
#' @param synonyms Named character vector `c(recorded = "counted as")`, with
#'   column or original species names. A target that is not yet a species is
#'   created.
#' @param note Optional text for the processing log.
#' @return The `cpue_data` object with merged species.
#' @seealso [species_consistency()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' x2 <- harmonise_species(x, stats::setNames(x$species[1], x$species[2]))
#' x2$log
harmonise_species <- function(x, synonyms, note = NULL) {
  check_cpue_data(x)
  if (is.null(names(synonyms)) || any(names(synonyms) == "")) {
    stop("`synonyms` must be a named vector: c(recorded = \"counted as\").", call. = FALSE)
  }
  resolve <- function(s) {
    if (s %in% x$species) return(s)
    hit <- match(s, x$species_labels)
    if (!is.na(hit)) names(x$species_labels)[hit] else NA_character_
  }
  from <- vapply(names(synonyms), resolve, character(1))
  if (anyNA(from)) {
    message("Not found (skipped): ", paste(names(synonyms)[is.na(from)], collapse = ", "))
  }
  ok <- !is.na(from)
  d <- x$data
  n_merged <- 0
  for (i in which(ok)) {
    f <- from[[i]]
    to_name <- unname(synonyms[i])
    t_code <- resolve(to_name)
    if (!is.na(t_code) && t_code == f) next
    if (is.na(t_code)) {
      t_code <- make.names(to_name)
      if (t_code %in% names(d)) t_code <- paste0(t_code, "_sp")
      d[[t_code]] <- 0
      x$species <- c(x$species, t_code)
      x$species_labels <- c(x$species_labels, stats::setNames(to_name, t_code))
    }
    d[[t_code]] <- d[[t_code]] + d[[f]]
    d[[f]] <- NULL
    x$species <- setdiff(x$species, f)
    x$species_labels <- x$species_labels[names(x$species_labels) != f]
    n_merged <- n_merged + 1
  }
  ord <- order(colSums(as.matrix(d[, x$species, drop = FALSE])), decreasing = TRUE)
  x$species <- x$species[ord]
  x$species_labels <- x$species_labels[x$species]
  x$data <- d
  if (is.null(note)) note <- sprintf("%d species codes merged", n_merged)
  x$log <- rbind(x$log, data.frame(step = "harmonise species", operations = nrow(d), note = note))
  if (!is.null(x$targeting)) {
    message("Targeting covariates were computed before this step; rerun add_targeting().")
  }
  x
}

#' Check the consistency of the vessel set over time
#'
#' A random vessel effect separates vessel efficiency from year effects only
#' if vessels overlap in time. Abrupt changes of the vessel set, e.g. a year
#' in which most operations come from vessels never seen before (a change of
#' data source or extraction), break this link and can shift the index.
#'
#' For each year the table gives the number of active, new (first year in the
#' data) and leaving vessels (active in the previous year but not in this
#' one), the Jaccard overlap of the vessel set with the previous year, the
#' share of operations by new vessels and by vessels recorded in that year
#' only. Years are flagged when the share of operations by new vessels
#' exceeds `max_new_share`, the overlap is below `min_overlap`, or the
#' number of vessels differs from the median of all years by more than a
#' factor `max_ratio`.
#'
#' @param x A `cpue_data` object with a vessel column.
#' @param max_new_share,min_overlap,max_ratio Flagging thresholds.
#' @param cohort_width Width (years) of the cohorts of first appearance used
#'   in the `cohort` table and plot.
#' @return An object of class `vessel_consistency`: a list with `years` (the
#'   yearly table with a `flag` column), `flags` (flagged years with the
#'   reason) and `cohort` (share of each year's operations by cohort of first
#'   appearance).
#' @seealso [plot_vessel_consistency()], [plot_fleet()], [data_consistency()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' vessel_consistency(x)
vessel_consistency <- function(x, max_new_share = 0.5, min_overlap = 0.3, max_ratio = 2,
                               cohort_width = 5) {
  check_cpue_data(x)
  r <- x$roles
  if (is.null(r$vessel)) stop("No vessel column was defined in cpue_data().", call. = FALSE)
  v <- as.character(x$data[[r$vessel]])
  yr <- x$data[[r$year]]
  years <- sort(unique(yr))
  first <- tapply(yr, v, min)
  n_years_v <- tapply(yr, v, function(z) length(unique(z)))
  sets <- split(v, factor(yr, levels = years))
  sets <- lapply(sets, unique)
  out <- data.frame(year = years, operations = as.numeric(table(factor(yr, levels = years))))
  out$vessels <- vapply(sets, length, 1L)
  out$new <- vapply(seq_along(years), function(j) sum(first[sets[[j]]] == years[j]), 1L)
  out$leaving <- c(NA, vapply(2:length(years), function(j)
    length(setdiff(sets[[j - 1]], sets[[j]])), 1L))
  out$overlap <- c(NA, vapply(2:length(years), function(j) {
    a <- sets[[j - 1]]; b <- sets[[j]]
    length(intersect(a, b)) / length(union(a, b))
  }, numeric(1)))
  out$ops_new <- as.numeric(tapply(first[v] == yr, factor(yr, levels = years), mean))
  out$ops_single_year <- as.numeric(tapply(n_years_v[v] == 1, factor(yr, levels = years), mean))
  out$new[1] <- NA
  out$ops_new[1] <- NA
  med <- stats::median(out$vessels)
  reasons <- lapply(seq_along(years), function(j) {
    rs <- character(0)
    if (!is.na(out$ops_new[j]) && out$ops_new[j] > max_new_share)
      rs <- c(rs, sprintf("%.0f%% of operations by new vessels", 100 * out$ops_new[j]))
    if (!is.na(out$overlap[j]) && out$overlap[j] < min_overlap)
      rs <- c(rs, sprintf("overlap with previous year %.2f", out$overlap[j]))
    if (out$vessels[j] > max_ratio * med || out$vessels[j] < med / max_ratio)
      rs <- c(rs, sprintf("%d vessels (median %g)", out$vessels[j], med))
    paste(rs, collapse = "; ")
  })
  out$flag <- unlist(reasons) != ""
  flags <- data.frame(year = years[out$flag], reason = unlist(reasons)[out$flag])
  cstart <- min(years) + cohort_width * floor((first - min(years)) / cohort_width)
  clab <- stats::setNames(paste0(cstart, "-", pmin(cstart + cohort_width - 1, max(years))),
                          names(first))
  cohort <- prop.table(table(year = yr, cohort = clab[v]), 1)
  structure(list(years = out, flags = flags, cohort = cohort,
                 settings = list(max_new_share = max_new_share, min_overlap = min_overlap,
                                 max_ratio = max_ratio)),
            class = "vessel_consistency")
}

#' @export
print.vessel_consistency <- function(x, ...) {
  y <- x$years
  cat("<vessel_consistency>", nrow(y), "years, median", stats::median(y$vessels),
      "vessels per year,", sum(y$flag), "year(s) flagged\n")
  if (nrow(x$flags)) print(x$flags, row.names = FALSE)
  invisible(x)
}

#' Plot the vessel-set consistency check
#'
#' @param x A [vessel_consistency()] result.
#' @param type `"turnover"`: active, new and leaving vessels, the overlap with
#'   the previous year and the share of operations by new and single-year
#'   vessels (flagged years in red); `"cohort"`: share of each year's
#'   operations by cohort of first appearance.
#' @return A ggplot.
#' @export
plot_vessel_consistency <- function(x, type = c("turnover", "cohort")) {
  if (!inherits(x, "vessel_consistency")) stop("`x` must come from vessel_consistency().", call. = FALSE)
  type <- match.arg(type)
  y <- x$years
  if (type == "turnover") {
    d <- rbind(
      data.frame(year = y$year, value = y$vessels, series = "active", panel = "Vessels"),
      data.frame(year = y$year, value = y$new, series = "new", panel = "Vessels"),
      data.frame(year = y$year, value = y$leaving, series = "leaving", panel = "Vessels"),
      data.frame(year = y$year, value = y$overlap, series = "overlap (Jaccard)", panel = "Overlap with previous year"),
      data.frame(year = y$year, value = y$ops_new, series = "by new vessels", panel = "Share of operations"),
      data.frame(year = y$year, value = y$ops_single_year, series = "by single-year vessels", panel = "Share of operations"))
    d$panel <- factor(d$panel, levels = c("Vessels", "Overlap with previous year", "Share of operations"))
    d <- d[!is.na(d$value), ]
    fl <- d[d$year %in% y$year[y$flag], ]
    return(ggplot(d, aes(x = .data$year, y = .data$value, colour = .data$series)) +
             geom_line() + geom_point(size = 1) +
             geom_point(data = fl, colour = "red", size = 2.5, shape = 1) +
             facet_wrap(~panel, ncol = 1, scales = "free_y") +
             labs(x = "Year", y = NULL, colour = NULL,
                  caption = if (any(y$flag)) "Red circles: flagged years") + theme_bw())
  }
  d <- as.data.frame(x$cohort, stringsAsFactors = FALSE)
  d$year <- as.numeric(d$year)
  ggplot(d, aes(x = .data$year, y = .data$Freq, fill = .data$cohort)) +
    geom_col(width = 1) +
    labs(x = "Year", y = "Share of operations", fill = "First year\nof vessel") + theme_bw()
}

#' Check the consistency of the catch composition and recording over time
#'
#' Yearly indicators that reveal changes in how catches were extracted or
#' recorded rather than in the fishery, e.g. a year extracted for the target
#' species only (every operation positive, one species per operation), a
#' change of recording precision (whole kilograms to decimals) or of units.
#'
#' Indicators per year: number of species recorded, mean number of species
#' per operation, share of single-species operations, Shannon diversity of
#' the total catch, median total catch per operation, share of catch values
#' that are whole numbers and, with `target`, the share of operations
#' catching the target and its share of the catch. A year is flagged for an
#' indicator when its robust z-score (deviation from the median of all years
#' divided by 1.4826 x MAD) exceeds `z` **and** it differs from the median by
#' more than a fraction `min_change` of the median (at least 0.05; so that tiny deviations
#' in very stable indicators are not flagged).
#'
#' @param x A `cpue_data` object (run before [select_species()]).
#' @param target Optional target species (column or original name).
#' @param z Threshold of the robust z-score.
#' @param min_change Minimum relative difference from the median.
#' @return An object of class `composition_consistency`: a list with
#'   `years` (one row per year and indicator), `table` (one row per year,
#'   one column per indicator) and `flags`.
#' @seealso [plot_composition_consistency()], [data_consistency()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' composition_consistency(x, target = x$species[1])
composition_consistency <- function(x, target = NULL, z = 3, min_change = 0.2) {
  check_cpue_data(x)
  yr <- x$data[[x$roles$year]]
  years <- sort(unique(yr))
  m <- as.matrix(x$data[, x$species, drop = FALSE])
  pos <- m > 0
  fy <- factor(yr, levels = years)
  nsp <- rowSums(pos)
  tot <- rowSums(m)
  tab <- data.frame(year = years)
  tab$`species recorded` <- vapply(split(seq_along(yr), fy), function(i)
    sum(colSums(pos[i, , drop = FALSE]) > 0), numeric(1))
  tab$`species per operation` <- as.numeric(tapply(nsp, fy, mean))
  tab$`single-species operations` <- as.numeric(tapply(nsp == 1, fy, mean))
  tab$`Shannon diversity` <- vapply(split(seq_along(yr), fy), function(i) {
    p <- colSums(m[i, , drop = FALSE]); p <- p[p > 0] / sum(p)
    -sum(p * log(p))
  }, numeric(1))
  tab$`median catch per operation` <- as.numeric(tapply(tot, fy, stats::median))
  vals <- m[pos]; vy <- yr[row(m)[pos]]
  tab$`whole-number catches` <- as.numeric(tapply(abs(vals - round(vals)) < 1e-8,
                                                 factor(vy, levels = years), mean))
  if (!is.null(target)) {
    tc <- species_codes(x, target)
    tcatch <- rowSums(m[, tc, drop = FALSE])
    tab$`operations with target` <- as.numeric(tapply(tcatch > 0, fy, mean))
    tab$`target share of catch` <- as.numeric(tapply(tcatch, fy, sum) / tapply(tot, fy, sum))
  }
  ind <- setdiff(names(tab), "year")
  long <- do.call(rbind, lapply(ind, function(k) {
    v <- tab[[k]]
    med <- stats::median(v, na.rm = TRUE)
    s <- 1.4826 * stats::mad(v, constant = 1, na.rm = TRUE)
    zz <- if (is.na(s) || s == 0) ifelse(abs(v - med) > 1e-8, Inf * sign(v - med), 0) else (v - med) / s
    rel <- abs(v - med) / max(abs(med), 0.05)
    data.frame(year = years, indicator = k, value = v, median = med, z = zz,
               flag = abs(zz) > z & rel > min_change, stringsAsFactors = FALSE)
  }))
  long$indicator <- factor(long$indicator, levels = ind)
  flags <- long[long$flag, c("year", "indicator", "value", "median", "z")]
  rownames(flags) <- NULL
  structure(list(years = long, table = tab, flags = flags, settings = list(z = z, min_change = min_change, target = target)),
            class = "composition_consistency")
}

#' @export
print.composition_consistency <- function(x, ...) {
  cat("<composition_consistency>", nrow(x$table), "years,", nrow(x$flags),
      "flag(s) (|robust z| >", x$settings$z, ")\n")
  if (nrow(x$flags)) {
    f <- x$flags
    f$value <- signif(f$value, 3); f$median <- signif(f$median, 3); f$z <- round(f$z, 1)
    print(f, row.names = FALSE)
  }
  invisible(x)
}

#' Plot the composition and recording consistency check
#'
#' @param x A [composition_consistency()] result.
#' @param indicators Optional subset of indicators.
#' @return A ggplot with one panel per indicator; flagged years in red and the
#'   median of all years as a dashed line.
#' @export
plot_composition_consistency <- function(x, indicators = NULL) {
  if (!inherits(x, "composition_consistency")) stop("`x` must come from composition_consistency().", call. = FALSE)
  d <- x$years
  if (!is.null(indicators)) d <- d[d$indicator %in% indicators, ]
  med <- unique(d[, c("indicator", "median")])
  ggplot(d, aes(x = .data$year, y = .data$value)) +
    geom_hline(data = med, aes(yintercept = .data$median), linetype = 2, colour = "grey50") +
    geom_line() + geom_point(size = 1) +
    geom_point(data = d[d$flag, ], colour = "red", size = 2.5) +
    facet_wrap(~indicator, scales = "free_y", ncol = 2) +
    labs(x = "Year", y = NULL, caption = "Red: flagged years; dashed: median of all years") +
    theme_bw()
}

#' Run all consistency checks of a data set
#'
#' Convenience wrapper around [species_consistency()],
#' [vessel_consistency()] (if a vessel column exists) and
#' [composition_consistency()], returning one table of all flags. Run it on
#' the full data (before [select_species()]) as the first step after
#' [cpue_data()].
#'
#' @param x A `cpue_data` object.
#' @param target Optional target species.
#' @param ... Passed to [species_consistency()].
#' @return An object of class `data_consistency`: a list with `species`,
#'   `vessels`, `composition` (the three results) and `flags`, a data frame
#'   with columns `check`, `year` and `issue`.
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' data_consistency(x, target = x$species[1])
data_consistency <- function(x, target = NULL, ...) {
  check_cpue_data(x)
  sc <- species_consistency(x, ...)
  vc <- if (!is.null(x$roles$vessel)) vessel_consistency(x) else NULL
  cc <- composition_consistency(x, target = target)
  fl <- list()
  if (nrow(sc$breaks)) {
    fl$species <- data.frame(check = "species recording", year = sc$breaks$year,
                             issue = sprintf("%s: %d species (%s)", sc$breaks$type,
                                             sc$breaks$n_species, sc$breaks$species))
  }
  if (!is.null(vc) && nrow(vc$flags)) {
    fl$vessels <- data.frame(check = "vessel set", year = vc$flags$year, issue = vc$flags$reason)
  }
  if (nrow(cc$flags)) {
    fl$composition <- data.frame(check = "composition", year = cc$flags$year,
                                 issue = sprintf("%s = %s (median %s)", cc$flags$indicator,
                                                 signif(cc$flags$value, 3), signif(cc$flags$median, 3)))
  }
  flags <- if (length(fl)) do.call(rbind, fl) else
    data.frame(check = character(0), year = numeric(0), issue = character(0))
  flags <- flags[order(flags$year, flags$check), ]
  rownames(flags) <- NULL
  structure(list(species = sc, vessels = vc, composition = cc, flags = flags),
            class = "data_consistency")
}

#' @export
print.data_consistency <- function(x, ...) {
  cat("<data_consistency>", nrow(x$flags), "flag(s)\n")
  if (nrow(x$flags)) {
    cat(sprintf("  %-6g %-18s %s", x$flags$year, x$flags$check, x$flags$issue), sep = "\n")
  }
  if (length(x$species$suggested)) {
    cat("Suggested species harmonisation (confirm with the code lists):\n")
    cat(paste0("  ", names(x$species$suggested), " -> ", x$species$suggested, collapse = "\n"), "\n")
  }
  invisible(x)
}
