#' Prepare catch and effort data for analysis
#'
#' First step of the real-data workflow. Brings logbook or observer records
#' into one standard object: one row per fishing operation (trip, set or
#' haul) and one column per species with the catch. All later steps
#' ([species_curve()], [stephens_maccall()], [add_targeting()], ...) work on
#' this object, so they know which columns are species and which play the
#' roles of year, vessel, effort and position.
#'
#' Two input formats are accepted:
#' * **long**: one row per operation and species, e.g. the usual logbook
#'   extract with columns `trip`, `year`, `vessel`, `species`, `catch_kg`.
#'   Give the name of the species column in `species` and of the catch column
#'   in `catch`. Catches of the same species within an operation are summed.
#' * **wide**: one row per operation with one column per species. Give the
#'   names of the species columns in `species` and leave `catch` empty.
#'
#' Operations are identified by `id` (one or more columns, e.g.
#' `c("vessel", "trip")`). Without `id`, every unique combination of the
#' operation-level columns is treated as one operation, which can merge
#' different operations that happen to share all values; a message reports
#' this case.
#'
#' Species names are converted to syntactic column names (e.g.
#' `"Pagellus bogaraveo"` becomes `Pagellus.bogaraveo`) so that they can be
#' used in model formulas; the original names are kept in `species_labels`
#' and used in plots.
#'
#' @param data A data frame.
#' @param year Name of the year column (numeric).
#' @param species Long format: name of the species column. Wide format:
#'   names of the species catch columns.
#' @param catch Long format only: name of the catch column (weight or
#'   numbers).
#' @param effort Optional name of the effort column (e.g. hours fished,
#'   hooks/1000, days). Operations with missing or non-positive effort are
#'   removed. Without effort, the catch per operation is the CPUE.
#' @param id Optional name(s) of the column(s) identifying an operation.
#' @param vessel Optional name of the vessel column.
#' @param lon,lat Optional names of position columns (decimal degrees or
#'   grid coordinates).
#' @param covariates Optional names of other operation-level columns to keep
#'   (e.g. `"month"`, `"port"`, `"HBF"`, `"gear"`).
#' @param synonyms Optional named character vector to merge species before
#'   anything else, e.g. `c("Trichiurus lepturus" = "Lepidopus caudatus")`
#'   (names: as recorded; values: as they should be counted).
#' @param drop_species Optional species names to drop entirely (e.g. bait or
#'   discarded categories).
#'
#' @return An object of class `cpue_data`: a list with
#'   * `data`: one row per operation with the role columns, `effort` if
#'     given, and one catch column per species;
#'   * `species`: the species column names, ordered by total catch;
#'   * `species_labels`: original species names (named by column name);
#'   * `roles`: which column plays which role;
#'   * `log`: a record of the processing steps and how many operations each
#'     removed (printed by `print()`).
#' @seealso [species_curve()], [select_species()], [add_targeting()]
#' @export
#' @examples
#' logbook <- data.frame(
#'   trip = rep(1:3, each = 2), year = rep(c(2020, 2020, 2021), each = 2),
#'   vessel = rep(c("A", "B", "A"), each = 2), hours = rep(c(5, 8, 6), each = 2),
#'   species = c("Pagellus bogaraveo", "Conger conger", "Pagellus bogaraveo",
#'               "Trichiurus lepturus", "Lepidopus caudatus", "Conger conger"),
#'   kg = c(20, 5, 12, 3, 4, 9))
#' x <- cpue_data(logbook, year = "year", species = "species", catch = "kg",
#'                effort = "hours", id = "trip", vessel = "vessel",
#'                synonyms = c("Trichiurus lepturus" = "Lepidopus caudatus"))
#' x
#' head(x$data)
cpue_data <- function(data, year, species, catch = NULL, effort = NULL, id = NULL,
                      vessel = NULL, lon = NULL, lat = NULL, covariates = NULL,
                      synonyms = NULL, drop_species = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  data <- as.data.frame(data)
  format <- if (!is.null(catch)) "long" else "wide"
  if (format == "long" && length(species) != 1) {
    stop("In long format, `species` is the name of one column (the species ",
         "names) and `catch` the name of the catch column.", call. = FALSE)
  }
  roles <- list(id = id, year = year, vessel = vessel, effort = effort,
                lon = lon, lat = lat, covariates = covariates)
  op_cols <- unique(unlist(roles))
  needed <- unique(c(op_cols, species, catch))
  missing <- setdiff(needed, names(data))
  if (length(missing)) {
    stop("Columns not found in `data`: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(data[[year]])) {
    year_num <- suppressWarnings(as.numeric(as.character(data[[year]])))
    if (anyNA(year_num[!is.na(data[[year]])])) {
      stop("The year column `", year, "` must be numeric.", call. = FALSE)
    }
    data[[year]] <- year_num
  }
  log <- data.frame(step = character(0), operations = integer(0), note = character(0))
  add_log <- function(log, step, n, note = "") {
    rbind(log, data.frame(step = step, operations = n, note = note))
  }

  if (format == "long") {
    sp <- as.character(data[[species]])
    if (!is.null(synonyms)) sp <- recode_species(sp, synonyms)
    keep <- !is.na(sp) & !(sp %in% drop_species)
    value <- data[[catch]]
    if (!is.numeric(value)) stop("The catch column must be numeric.", call. = FALSE)
    if (any(value < 0, na.rm = TRUE)) stop("Negative catches found.", call. = FALSE)
    value[is.na(value)] <- 0
    wide <- long_to_wide(data[keep, op_cols, drop = FALSE], sp[keep], value[keep], id)
    sets <- wide$sets
    catch_mat <- wide$catch
    log <- add_log(log, "input (long format)", nrow(sets),
                   sprintf("%d records, %d species", sum(keep), ncol(catch_mat)))
  } else {
    sets <- data[, op_cols, drop = FALSE]
    catch_mat <- as.matrix(data[, species, drop = FALSE])
    if (!is.numeric(catch_mat)) stop("Species columns must be numeric.", call. = FALSE)
    if (any(catch_mat < 0, na.rm = TRUE)) stop("Negative catches found.", call. = FALSE)
    catch_mat[is.na(catch_mat)] <- 0
    colnames(catch_mat) <- species
    if (!is.null(synonyms)) {
      merged <- recode_species(colnames(catch_mat), synonyms)
      catch_mat <- t(rowsum(t(catch_mat), merged, reorder = FALSE))
    }
    catch_mat <- catch_mat[, !colnames(catch_mat) %in% drop_species, drop = FALSE]
    log <- add_log(log, "input (wide format)", nrow(sets),
                   sprintf("%d species", ncol(catch_mat)))
  }

  # Order species by total catch and make syntactic column names
  ord <- order(colSums(catch_mat), decreasing = TRUE)
  catch_mat <- catch_mat[, ord, drop = FALSE]
  labels <- colnames(catch_mat)
  codes <- make.names(labels, unique = TRUE)
  clash <- intersect(codes, c(op_cols, "effort"))
  if (length(clash)) codes[codes %in% clash] <- paste0(codes[codes %in% clash], "_sp")
  colnames(catch_mat) <- codes
  names(labels) <- codes

  out <- sets
  if (!is.null(effort)) {
    eff <- out[[effort]]
    bad <- is.na(eff) | eff <= 0
    if (any(bad)) {
      out <- out[!bad, , drop = FALSE]
      catch_mat <- catch_mat[!bad, , drop = FALSE]
      log <- add_log(log, "remove missing or zero effort", nrow(out),
                     sprintf("%d removed", sum(bad)))
    }
  }
  bad_year <- is.na(out[[year]])
  if (any(bad_year)) {
    out <- out[!bad_year, , drop = FALSE]
    catch_mat <- catch_mat[!bad_year, , drop = FALSE]
    log <- add_log(log, "remove missing year", nrow(out), sprintf("%d removed", sum(bad_year)))
  }
  out <- cbind(out, as.data.frame(catch_mat))
  rownames(out) <- NULL
  structure(list(data = out, species = codes, species_labels = labels,
                 roles = roles, log = log),
            class = "cpue_data")
}

#' Merge species names
#' @noRd
recode_species <- function(sp, synonyms) {
  if (is.null(names(synonyms)) || any(names(synonyms) == "")) {
    stop("`synonyms` must be a named vector: c(recorded = \"counted as\").", call. = FALSE)
  }
  hit <- sp %in% names(synonyms)
  sp[hit] <- unname(synonyms[sp[hit]])
  sp
}

#' Long records (operation x species) to one row per operation
#' @noRd
long_to_wide <- function(sets_long, sp, value, id = NULL) {
  key_cols <- if (is.null(id)) names(sets_long) else id
  key <- do.call(paste, c(lapply(sets_long[key_cols], as.character), sep = "\r"))
  ukey <- unique(key)
  first <- match(ukey, key)
  sets <- sets_long[first, , drop = FALSE]
  if (!is.null(id)) {
    other <- setdiff(names(sets_long), id)
    if (length(other)) {
      full <- do.call(paste, c(lapply(sets_long[c(id, other)], as.character), sep = "\r"))
      n_full <- length(unique(full))
      if (n_full > length(ukey)) {
        message(n_full - length(ukey), " operation(s) have more than one value of an ",
                "operation-level column within `id`; the first value is kept.")
      }
    }
  } else {
    message("No `id` given: each unique combination of the operation-level columns ",
            "is treated as one operation.")
  }
  species <- sort(unique(sp))
  tab <- tapply(value, list(factor(key, levels = ukey), factor(sp, levels = species)),
                sum)
  tab[is.na(tab)] <- 0
  dimnames(tab) <- list(NULL, species)
  rownames(sets) <- NULL
  list(sets = sets, catch = unclass(tab))
}

#' @export
print.cpue_data <- function(x, ...) {
  d <- x$data
  r <- x$roles
  yrs <- range(d[[r$year]])
  cat("<cpue_data>", nrow(d), "operations,", length(x$species), "species,",
      "years", yrs[1], "-", yrs[2], "\n")
  if (!is.null(r$vessel)) cat("  vessels:", length(unique(d[[r$vessel]])), "\n")
  cat("  effort:", if (is.null(r$effort)) "none (catch per operation)" else r$effort, "\n")
  if (!is.null(r$lon)) cat("  position:", r$lon, "/", r$lat, "\n")
  if (length(r$covariates)) cat("  covariates:", paste(r$covariates, collapse = ", "), "\n")
  top <- utils::head(x$species, 5)
  cat("  main species:", paste(x$species_labels[top], collapse = ", "),
      if (length(x$species) > 5) "..." else "", "\n")
  if (!is.null(x$targeting)) {
    cat("  targeting: ", x$targeting$n_pc, " PC(s) used, ", x$targeting$n_clusters,
        " clusters (", x$targeting$method, ")\n", sep = "")
  }
  cat("Processing log:\n")
  print(x$log, row.names = FALSE)
  invisible(x)
}

#' Keep a subset of operations
#'
#' Works like [base::subset()] on the operation table and records the step in
#' the processing log. `subset` can be a logical expression evaluated within
#' `x$data` (e.g. `year >= 2005 & month %in% 4:9`) or a logical vector, e.g.
#' the `selected` element of [stephens_maccall()].
#'
#' @param x A `cpue_data` object.
#' @param subset Logical expression or vector.
#' @param note Optional text for the processing log.
#' @param ... Not used.
#' @return The subsetted `cpue_data` object.
#' @export
#' @examples
#' \dontrun{
#' x2 <- subset(x, year >= 2010, note = "recent years")
#' }
subset.cpue_data <- function(x, subset, note = NULL, ...) {
  keep <- eval(substitute(subset), x$data, parent.frame())
  if (!is.logical(keep) || length(keep) != nrow(x$data)) {
    stop("`subset` must give one TRUE/FALSE per operation.", call. = FALSE)
  }
  keep[is.na(keep)] <- FALSE
  if (is.null(note)) note <- paste(deparse(substitute(subset)), collapse = " ")
  if (nchar(note) > 60) note <- paste0(substr(note, 1, 57), "...")
  filter_ops(x, keep, "subset", note)
}

#' Keep operations and update the log and targeting results
#' @noRd
filter_ops <- function(x, keep, step, note = "") {
  x$data <- x$data[keep, , drop = FALSE]
  rownames(x$data) <- NULL
  x$log <- rbind(x$log, data.frame(step = step, operations = nrow(x$data),
                                   note = sprintf("%s (%d removed)", note, sum(!keep))))
  if (!is.null(x$targeting)) {
    message("Targeting covariates were computed before this step; ",
            "rerun add_targeting() on the subset.")
  }
  x
}

#' CPUE of each operation
#'
#' Catch divided by effort (or the catch itself when no effort column was
#' given).
#'
#' @param x A `cpue_data` object.
#' @param species Species column names (default: all).
#' @return A matrix with one column per species.
#' @export
cpue_matrix <- function(x, species = x$species) {
  m <- as.matrix(x$data[, species, drop = FALSE])
  if (!is.null(x$roles$effort)) m <- m / x$data[[x$roles$effort]]
  m
}

#' Species label for plots
#' @noRd
sp_label <- function(x, codes) {
  out <- x$species_labels[codes]
  out[is.na(out)] <- codes[is.na(out)]
  unname(out)
}

#' Check the object class
#' @noRd
check_cpue_data <- function(x) {
  if (!inherits(x, "cpue_data")) {
    stop("`x` must be a cpue_data object; create it with cpue_data().", call. = FALSE)
  }
  invisible(x)
}
