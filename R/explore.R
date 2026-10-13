#' Cumulative species curve
#'
#' Ranks species by how often they are caught (`by = "occurrence"`, the share
#' of all positive species records, as in the CPUE guidelines) or by total
#' catch (`by = "weight"`), with the cumulative share. Used to choose the set
#' of species that carries the targeting signal: rare species add noise to
#' the catch composition, so the PCA and clustering use only the species that
#' together make up most of the records (e.g. 95%).
#'
#' @param x A `cpue_data` object.
#' @param by `"occurrence"` (default) or `"weight"`.
#' @return A data frame with `species` (column name), `label`, `rank`,
#'   `value` (number of positive records or total catch), `share` and
#'   `cumulative`.
#' @seealso [plot_species_curve()], [select_species()]
#' @export
species_curve <- function(x, by = c("occurrence", "weight")) {
  check_cpue_data(x)
  by <- match.arg(by)
  m <- as.matrix(x$data[, x$species, drop = FALSE])
  value <- if (by == "occurrence") colSums(m > 0) else colSums(m)
  value <- sort(value, decreasing = TRUE)
  out <- data.frame(species = names(value), label = sp_label(x, names(value)),
                    rank = seq_along(value), value = unname(value))
  out$share <- out$value / sum(out$value)
  out$cumulative <- cumsum(out$share)
  out
}

#' Plot the cumulative species curve
#'
#' @inheritParams species_curve
#' @param threshold Cumulative share used to select species (drawn as a
#'   line; the selected species are coloured).
#' @param max_species Show at most this many species.
#' @return A ggplot object.
#' @export
plot_species_curve <- function(x, threshold = 0.95, by = c("occurrence", "weight"),
                               max_species = 30) {
  by <- match.arg(by)
  sc <- species_curve(x, by)
  sc$selected <- seq_len(nrow(sc)) <= n_species_for(sc, threshold)
  sc <- sc[seq_len(min(nrow(sc), max_species)), ]
  sc$label <- factor(sc$label, levels = sc$label)
  sel <- sc[sc$selected, ]
  ggplot(sc, aes(x = .data$rank, y = .data$cumulative)) +
    geom_line(colour = "grey50") +
    geom_point(colour = "grey50") +
    geom_point(data = sel, aes(colour = .data$label), size = 3) +
    geom_hline(yintercept = threshold, linetype = 2) +
    geom_vline(xintercept = sum(sc$selected) + 0.5, linetype = 3) +
    scale_x_continuous(breaks = sc$rank) +
    labs(x = "Species rank",
         y = paste("Cumulative share of", if (by == "occurrence") "records" else "catch"),
         colour = "Selected species") +
    theme_bw()
}

#' Number of species needed to reach a cumulative share
#' @noRd
n_species_for <- function(sc, threshold) {
  if (threshold >= 1) return(nrow(sc))
  min(nrow(sc), sum(sc$cumulative < threshold) + 1)
}

#' Select the main species
#'
#' Keeps the species that, ranked by [species_curve()], are needed to reach
#' the cumulative share `threshold`, plus any species listed in `keep` (e.g.
#' the species of interest if it is rare). The other species are dropped or
#' pooled into one column `Other`. Operations that caught none of the
#' remaining species are removed, because their catch composition is
#' undefined.
#'
#' Pooling keeps the information that an operation caught "something else",
#' which can help to separate tactics, but a large `Other` group can also
#' dominate the composition; try both.
#'
#' @inheritParams species_curve
#' @param threshold Cumulative share to reach (default 0.95).
#' @param keep Species to keep in any case (column names or original names).
#' @param n Alternatively, the number of species to keep (overrides
#'   `threshold`).
#' @param pool_other Logical; pool the dropped species into `Other`.
#' @return The `cpue_data` object restricted to the selected species.
#' @export
#' @examples
#' \dontrun{
#' plot_species_curve(x, threshold = 0.95)
#' x <- select_species(x, threshold = 0.95, keep = "Pagellus bogaraveo")
#' }
select_species <- function(x, threshold = 0.95, by = c("occurrence", "weight"),
                           keep = NULL, n = NULL, pool_other = FALSE) {
  check_cpue_data(x)
  by <- match.arg(by)
  sc <- species_curve(x, by)
  n_sel <- if (is.null(n)) n_species_for(sc, threshold) else min(n, nrow(sc))
  chosen <- sc$species[seq_len(n_sel)]
  if (!is.null(keep)) {
    keep_codes <- species_codes(x, keep)
    chosen <- union(chosen, keep_codes)
  }
  dropped <- setdiff(x$species, chosen)
  n_orig <- length(x$species)
  n_kept <- length(chosen)
  d <- x$data
  if (pool_other && length(dropped)) {
    d$Other <- rowSums(d[, dropped, drop = FALSE])
    x$species_labels <- c(x$species_labels, Other = "Other")
    chosen <- c(chosen, "Other")
  }
  d <- d[, setdiff(names(d), dropped), drop = FALSE]
  x$data <- d
  x$species <- chosen
  x$species_labels <- x$species_labels[chosen]
  note <- sprintf("%d of %d species kept (%s >= %s)%s", n_kept, n_orig, by,
                  if (is.null(n)) threshold else paste("top", n),
                  if (pool_other && length(dropped)) ", rest pooled as Other" else "")
  x$log <- rbind(x$log, data.frame(step = "select species", operations = nrow(d), note = note))
  empty <- rowSums(as.matrix(d[, chosen, drop = FALSE])) == 0
  if (any(empty)) x <- filter_ops(x, !empty, "remove operations without selected species",
                                  "zero catch of all selected species")
  x
}

#' Translate original species names or column names to column names
#' @noRd
species_codes <- function(x, species) {
  codes <- ifelse(species %in% x$species, species,
                  names(x$species_labels)[match(species, x$species_labels)])
  if (anyNA(codes)) {
    stop("Species not found: ", paste(species[is.na(codes)], collapse = ", "), call. = FALSE)
  }
  codes
}

#' Overview of the fleet and the records
#'
#' A set of exploratory plots of how the data are distributed over years,
#' vessels and covariates, i.e. the non-random allocation of effort that
#' standardisation has to deal with.
#'
#' @param x A `cpue_data` object.
#' @param by Optional covariate (e.g. `"month"`, `"port"`, `"HBF"`) for the
#'   `by_covariate` plot.
#' @param max_vessels Show at most this many vessels in the vessel plot (the
#'   ones with most records).
#' @return A list of ggplots: `records` (operations per year), `vessels`
#'   (vessel x year tile plot of the number of operations, showing fleet
#'   turnover), `vessel_years` (vessels active per year and new entries) and,
#'   if `by` is given, `by_covariate` (share of operations by covariate level
#'   over years).
#' @export
plot_fleet <- function(x, by = NULL, max_vessels = 60) {
  check_cpue_data(x)
  r <- x$roles
  d <- x$data
  yr <- d[[r$year]]
  out <- list()
  rec <- as.data.frame(table(year = yr), stringsAsFactors = FALSE)
  rec$year <- as.numeric(rec$year)
  out$records <- ggplot(rec, aes(x = .data$year, y = .data$Freq)) +
    geom_col(fill = "grey60", colour = "grey30") +
    labs(x = "Year", y = "Number of operations") + theme_bw()
  if (!is.null(r$vessel)) {
    v <- as.character(d[[r$vessel]])
    tab <- as.data.frame(table(vessel = v, year = yr), stringsAsFactors = FALSE)
    tab <- tab[tab$Freq > 0, ]
    tab$year <- as.numeric(tab$year)
    top <- names(sort(table(v), decreasing = TRUE))[seq_len(min(max_vessels, length(unique(v))))]
    tab_top <- tab[tab$vessel %in% top, ]
    first <- tapply(tab_top$year, tab_top$vessel, min)
    last <- tapply(tab_top$year, tab_top$vessel, max)
    lev <- names(first)[order(first, last)]
    tab_top$vessel <- factor(tab_top$vessel, levels = lev)
    out$vessels <- ggplot(tab_top, aes(x = .data$year, y = .data$vessel,
                                       fill = log10(.data$Freq))) +
      geom_tile(colour = "white") +
      scale_fill_viridis_c() +
      labs(x = "Year", y = "Vessel", fill = "log10\noperations",
           subtitle = if (length(unique(v)) > max_vessels)
             sprintf("%d vessels with most records", max_vessels)) +
      theme_bw() + theme(axis.text.y = element_text(size = 6))
    first_all <- tapply(tab$year, tab$vessel, min)
    vy <- data.frame(year = sort(unique(tab$year)))
    vy$active <- as.numeric(table(tab$year)[as.character(vy$year)])
    vy$new <- as.numeric(table(factor(first_all, levels = vy$year)))
    vy$new[1] <- NA
    vy_long <- rbind(data.frame(year = vy$year, n = vy$active, what = "Active"),
                     data.frame(year = vy$year, n = vy$new, what = "New"))
    out$vessel_years <- ggplot(vy_long[!is.na(vy_long$n), ],
                               aes(x = .data$year, y = .data$n, colour = .data$what)) +
      geom_line() + geom_point() +
      labs(x = "Year", y = "Number of vessels", colour = NULL) + theme_bw()
  }
  if (!is.null(by)) {
    if (!by %in% names(d)) stop("Column `", by, "` not found.", call. = FALSE)
    cov <- data.frame(year = yr, level = factor(d[[by]]))
    out$by_covariate <- ggplot(cov, aes(x = .data$year, fill = .data$level)) +
      geom_bar(position = "fill") +
      labs(x = "Year", y = "Share of operations", fill = by) + theme_bw()
  }
  out
}

#' Select core vessels
#'
#' Keeps vessels that fished in at least `min_years` years with at least
#' `min_records` operations in each of those years, a common way to define
#' an indicator fleet with consistent reporting.
#'
#' @param x A `cpue_data` object with a vessel column.
#' @param min_years Minimum number of qualifying years.
#' @param min_records Minimum number of operations for a year to qualify.
#' @return The subsetted `cpue_data` object.
#' @export
select_vessels <- function(x, min_years = 3, min_records = 10) {
  check_cpue_data(x)
  r <- x$roles
  if (is.null(r$vessel)) stop("No vessel column was defined in cpue_data().", call. = FALSE)
  v <- as.character(x$data[[r$vessel]])
  tab <- table(v, x$data[[r$year]])
  years_ok <- rowSums(tab >= min_records)
  core <- names(years_ok)[years_ok >= min_years]
  filter_ops(x, v %in% core, "select core vessels",
             sprintf("%d of %d vessels (>= %d years with >= %d records)",
                     length(core), nrow(tab), min_years, min_records))
}

#' Nominal CPUE by year
#'
#' @param x A `cpue_data` object.
#' @param species Species (column or original names; default: all).
#' @return A data frame with `species`, `year`, `n` (operations), `positive`
#'   (proportion of operations catching the species), `mean` (mean CPUE),
#'   `mean_positive` (mean CPUE of positive operations), `median` and
#'   `relative` (mean scaled to its average over years).
#' @export
nominal_cpue <- function(x, species = x$species) {
  check_cpue_data(x)
  species <- species_codes(x, species)
  cp <- cpue_matrix(x, species)
  yr <- x$data[[x$roles$year]]
  out <- do.call(rbind, lapply(species, function(s) {
    y <- cp[, s]
    years <- sort(unique(yr))
    data.frame(species = sp_label(x, s), year = years,
               n = as.numeric(table(yr)[as.character(years)]),
               positive = as.numeric(tapply(y > 0, yr, mean)[as.character(years)]),
               mean = as.numeric(tapply(y, yr, mean)[as.character(years)]),
               mean_positive = as.numeric(tapply(y, yr, function(z) {
                 if (any(z > 0)) mean(z[z > 0]) else NA_real_
               })[as.character(years)]),
               median = as.numeric(tapply(y, yr, stats::median)[as.character(years)]))
  }))
  out$relative <- out$mean / stats::ave(out$mean, out$species, FUN = function(z) mean(z, na.rm = TRUE))
  out$species <- factor(out$species, levels = sp_label(x, species))
  out
}

#' Plot nominal CPUE
#'
#' @inheritParams nominal_cpue
#' @return A list of ggplots: `index` (relative nominal CPUE by year),
#'   `positive` (proportion of positive operations by year), `distribution`
#'   (histogram of log CPUE of positive operations) and `boxplot` (CPUE of
#'   positive operations by year, log scale).
#' @export
plot_nominal <- function(x, species = x$species) {
  nom <- nominal_cpue(x, species)
  codes <- species_codes(x, species)
  cp <- cpue_matrix(x, codes)
  yr <- x$data[[x$roles$year]]
  long <- data.frame(species = factor(rep(sp_label(x, codes), each = nrow(cp)),
                                      levels = sp_label(x, codes)),
                     year = rep(yr, length(codes)), cpue = as.vector(cp))
  pos <- long[long$cpue > 0, ]
  list(
    index = ggplot(nom, aes(x = .data$year, y = .data$relative)) +
      geom_line() + geom_point() + geom_hline(yintercept = 1, linetype = 2) +
      facet_wrap(~species, scales = "free_y") +
      labs(x = "Year", y = "Relative nominal CPUE") + theme_bw(),
    positive = ggplot(nom, aes(x = .data$year, y = .data$positive)) +
      geom_line() + geom_point() + facet_wrap(~species) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Year", y = "Proportion of positive operations") + theme_bw(),
    distribution = ggplot(pos, aes(x = log(.data$cpue))) +
      geom_histogram(bins = 50, fill = "grey60", colour = "grey30") +
      facet_wrap(~species, scales = "free") +
      labs(x = "log CPUE (positive operations)", y = "Frequency") + theme_bw(),
    boxplot = ggplot(pos, aes(x = factor(.data$year), y = .data$cpue)) +
      geom_boxplot(outlier.size = 0.4, fill = "grey85") +
      scale_y_log10() + facet_wrap(~species, scales = "free_y") +
      labs(x = "Year", y = "CPUE (positive operations, log scale)") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
  )
}

#' Plot the catch composition over time
#'
#' Share of each species in the total catch (`type = "weight"`) or in the
#' positive records (`type = "occurrence"`) by year or by another grouping
#' column. Changes in composition are the first sign of changes in targeting
#' (or in relative abundance).
#'
#' @param x A `cpue_data` object.
#' @param type `"weight"` or `"occurrence"`.
#' @param by Grouping column (default: the year column).
#' @return A ggplot object.
#' @export
plot_composition <- function(x, type = c("weight", "occurrence"), by = NULL) {
  check_cpue_data(x)
  type <- match.arg(type)
  if (is.null(by)) by <- x$roles$year
  m <- as.matrix(x$data[, x$species, drop = FALSE])
  if (type == "occurrence") m <- (m > 0) * 1
  g <- x$data[[by]]
  tot <- rowsum(m, g)
  share <- tot / rowSums(tot)
  long <- data.frame(group = rep(rownames(share), ncol(share)),
                     species = factor(rep(sp_label(x, colnames(share)), each = nrow(share)),
                                      levels = sp_label(x, x$species)),
                     share = as.vector(share))
  if (is.numeric(g)) long$group <- as.numeric(long$group)
  ggplot(long, aes(x = .data$group, y = .data$share, fill = .data$species)) +
    geom_col(width = if (is.numeric(g)) 0.9 else 0.8) +
    scale_fill_manual(values = species_colours(length(x$species))) +
    labs(x = by, y = paste("Share of", if (type == "weight") "catch" else "records"),
         fill = "Species") +
    theme_bw()
}

#' Map of operations and CPUE
#'
#' @param x A `cpue_data` object with `lon` and `lat`.
#' @param species Optional species for the CPUE map (column or original
#'   name); without it, the number of operations per cell is mapped.
#' @param res Cell size in units of `lon`/`lat` (e.g. 0.5 degrees).
#' @param by_period Optional number of years per panel (e.g. 5) to see shifts
#'   in fishing grounds.
#' @return A ggplot object.
#' @export
plot_cpue_map <- function(x, species = NULL, res = 1, by_period = NULL) {
  check_cpue_data(x)
  r <- x$roles
  if (is.null(r$lon) || is.null(r$lat)) {
    stop("No position columns were defined in cpue_data().", call. = FALSE)
  }
  d <- data.frame(lon = floor(x$data[[r$lon]] / res) * res + res / 2,
                  lat = floor(x$data[[r$lat]] / res) * res + res / 2,
                  year = x$data[[r$year]])
  d$period <- if (is.null(by_period)) "All years" else {
    start <- min(d$year) + ((d$year - min(d$year)) %/% by_period) * by_period
    paste0(start, "-", pmin(start + by_period - 1, max(d$year)))
  }
  if (is.null(species)) {
    d$value <- 1
    agg <- stats::aggregate(value ~ lon + lat + period, d, sum)
    fill_lab <- "Operations"
  } else {
    code <- species_codes(x, species)
    d$value <- cpue_matrix(x, code)[, 1]
    agg <- stats::aggregate(value ~ lon + lat + period, d, mean)
    agg$value <- log1p(agg$value)
    fill_lab <- paste0("log(1 + mean CPUE)\n", sp_label(x, code))
  }
  ggplot(agg, aes(x = .data$lon, y = .data$lat, fill = .data$value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "turbo") +
    coord_quickmap() +
    facet_wrap(~period) +
    labs(x = r$lon, y = r$lat, fill = fill_lab) + theme_bw()
}

#' Distinct colours for species
#' @noRd
species_colours <- function(n) {
  base <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
            "#a6761d", "#1f78b4", "#b2df8a", "#fb9a99", "#cab2d6", "#666666")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}
