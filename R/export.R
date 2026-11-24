# Export standardised indices to assessment-model inputs:
# JABBA / JABBA3 data frames, FLR (FLIndexBiomass / FLIndices) and SPiCT.

#' Collect standardised indices for export
#'
#' Harmonises one or several standardised indices into one long table, the
#' common input of the exporters [as_jabba()], [as_FLIndexBiomass()],
#' [as_FLIndices()] and [as_spict()].
#'
#' `x` can be
#' * a `cpue_fit` from [fit_cpue()]: the index of `model` (default: the last
#'   model fitted, i.e. the most complete); several models give several
#'   series named after the model;
#' * a named list of `cpue_fit` objects (e.g. scenarios): one series per
#'   element and model, named `<element>` (one model) or
#'   `<element>_<model>`;
#' * a data frame with columns `year`, `index` and `se_log` (log-scale SE),
#'   and optionally `name` (or `model`) for several series; or a named list of
#'   such data frames.
#'
#' @param x A `cpue_fit`, a data frame, or a named list of either.
#' @param model Model(s) to export from `cpue_fit` objects.
#' @param years Optional years to keep (e.g. to drop years with inconsistent
#'   recording).
#' @param normalise If `TRUE`, each series is divided by its mean over the
#'   years kept (the log-scale SE is unchanged).
#' @return A data frame of class `cpue_index_set` with columns `name`,
#'   `year`, `index` and `se_log`.
#' @seealso [as_jabba()], [as_FLIndices()], [as_spict()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' x <- add_targeting(select_species(x, keep = "Pagellus bogaraveo"))
#' fit <- fit_cpue(x, "Pagellus bogaraveo", models = c("Base", "Cluster"))
#' head(index_set(fit, model = c("Base", "Cluster")))
index_set <- function(x, model = NULL, years = NULL, normalise = FALSE) {
  one <- function(obj, nm) {
    if (inherits(obj, "cpue_fit")) {
      mods <- if (is.null(model)) utils::tail(names(obj$models), 1) else model
      miss <- setdiff(mods, unique(obj$indices$model))
      if (length(miss)) stop("Model(s) not in the fit: ", paste(miss, collapse = ", "), call. = FALSE)
      do.call(rbind, lapply(mods, function(m) {
        z <- obj$indices[obj$indices$model == m, ]
        lab <- if (is.null(nm)) m else if (length(mods) > 1) paste(nm, m, sep = "_") else nm
        data.frame(name = lab, year = z$year, index = z$index, se_log = z$se_log)
      }))
    } else if (is.data.frame(obj)) {
      need <- c("year", "index", "se_log")
      if (!all(need %in% names(obj))) {
        stop("Data frames need columns year, index and se_log.", call. = FALSE)
      }
      lab <- if ("name" %in% names(obj)) as.character(obj$name) else
        if ("model" %in% names(obj)) as.character(obj$model) else
          if (is.null(nm)) "index" else nm
      if (!is.null(nm) && any(c("name", "model") %in% names(obj))) lab <- paste(nm, lab, sep = "_")
      data.frame(name = lab, year = obj$year, index = obj$index, se_log = obj$se_log)
    } else {
      stop("Unsupported input: use a cpue_fit, a data frame or a named list of these.", call. = FALSE)
    }
  }
  if (inherits(x, "cpue_index_set")) {
    out <- as.data.frame(x)
  } else if (inherits(x, "cpue_fit") || is.data.frame(x)) {
    out <- one(x, NULL)
  } else if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms) || any(nms == "")) stop("A list of indices must be named.", call. = FALSE)
    out <- do.call(rbind, Map(one, x, nms))
  } else {
    stop("Unsupported input: use a cpue_fit, a data frame or a named list of these.", call. = FALSE)
  }
  out <- out[is.finite(out$index) & out$index > 0, ]
  if (!is.null(years)) out <- out[out$year %in% years, ]
  if (!nrow(out)) stop("No index values left.", call. = FALSE)
  if (any(duplicated(out[, c("name", "year")]))) {
    stop("Duplicated name-year combinations; give each series a unique name.", call. = FALSE)
  }
  if (normalise) out$index <- out$index / stats::ave(out$index, out$name)
  out$name <- factor(out$name, levels = unique(out$name))
  out <- out[order(out$name, out$year), ]
  out$name <- as.character(out$name)
  rownames(out) <- NULL
  class(out) <- c("cpue_index_set", "data.frame")
  out
}

#' Indices as JABBA / JABBA3 data frames
#'
#' Returns the `cpue` and `se` data frames used by `JABBA::build_jabba()`
#' and accepted by `JABBA3::jb3_data()` (`indices`, `se`): one row per year,
#' one column per index, `NA` where an index has no value. `se` holds the
#' log-scale standard errors.
#'
#' @inheritParams index_set
#' @param ... Passed to [index_set()] (`model`, `years`, `normalise`).
#' @return A list with data frames `cpue` and `se` (first column `year`).
#' @seealso [index_set()], [as_FLIndices()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' x <- add_targeting(select_species(x, keep = "Pagellus bogaraveo"))
#' fit <- fit_cpue(x, "Pagellus bogaraveo", models = c("Base", "Cluster"))
#' jb <- as_jabba(fit, model = c("Base", "Cluster"), normalise = TRUE)
#' head(jb$cpue); head(jb$se)
#' \dontrun{
#' # classic JABBA
#' JABBA::build_jabba(catch = catch, cpue = jb$cpue, se = jb$se, ...)
#' # JABBA3
#' d <- JABBA3::jb3_data(catch, jb$cpue, jb$se)
#' }
as_jabba <- function(x, ...) {
  s <- index_set(x, ...)
  yrs <- sort(unique(s$year))
  nms <- unique(s$name)
  cpue <- data.frame(year = yrs)
  se <- data.frame(year = yrs)
  for (n in nms) {
    z <- s[s$name == n, ]
    cpue[[n]] <- z$index[match(yrs, z$year)]
    se[[n]] <- z$se_log[match(yrs, z$year)]
  }
  list(cpue = cpue, se = se)
}

#' Indices as FLR FLIndexBiomass / FLIndices
#'
#' Converts standardised indices to `FLIndexBiomass` objects of the FLR
#' package `FLCore`, the native index input of `JABBA3` (`jb3_data()`), with
#' the index in `index` and the **log-scale standard error in `index.var`**
#' (the JABBA3 convention, also used by `ss3diags2::ss3_indices()`, with
#' `distribution = "lognormal"`). Several series (models or scenarios) are
#' returned as a named `FLIndices`. Requires the suggested package `FLCore`
#' (<https://flr-project.org>).
#'
#' The timing of the index within the year is stored in the `startf` and
#' `endf` elements of `range` (FLR convention). `JABBA3` takes the timing of
#' annual indices from the `index.time` argument of `jb3_data()` instead.
#'
#' @inheritParams index_set
#' @param ... Passed to [index_set()] (`model`, `years`, `normalise`).
#' @param units Units of the index.
#' @param startf,endf Start and end of the index period as fractions of the
#'   year (e.g. 0 and 2/12 for January-February); `NA` leaves them unset.
#' @param desc Description stored in the object.
#' @return `as_FLIndexBiomass()`: an `FLIndexBiomass` (one series only);
#'   `as_FLIndices()`: a named `FLIndices` of `FLIndexBiomass`.
#' @seealso [index_set()], [as_jabba()]
#' @export
#' @examples
#' \dontrun{
#' fi <- as_FLIndexBiomass(fit, model = "Cluster", startf = 0, endf = 2 / 12)
#' fl <- as_FLIndices(list(jf = fit_jf, winter = fit_winter), model = "Cluster")
#' d <- JABBA3::jb3_data(catch, fl, index.time = 1 / 12)
#' }
as_FLIndexBiomass <- function(x, ..., units = "NA", startf = NA, endf = NA, desc = "") {
  s <- index_set(x, ...)
  if (length(unique(s$name)) > 1) {
    stop("More than one series; use as_FLIndices().", call. = FALSE)
  }
  .flindex(s, units = units, startf = startf, endf = endf, desc = desc)
}

#' @rdname as_FLIndexBiomass
#' @export
as_FLIndices <- function(x, ..., units = "NA", startf = NA, endf = NA, desc = "") {
  s <- index_set(x, ...)
  nms <- unique(s$name)
  rec <- function(v) rep_len(v, length(nms))
  units <- rec(units); startf <- rec(startf); endf <- rec(endf); desc <- rec(desc)
  .need_flcore()
  out <- lapply(seq_along(nms), function(j)
    .flindex(s[s$name == nms[j], ], units = units[j], startf = startf[j], endf = endf[j],
             desc = desc[j]))
  FLCore::FLIndices(stats::setNames(out, nms))
}

.need_flcore <- function() {
  if (!requireNamespace("FLCore", quietly = TRUE)) {
    stop("Package 'FLCore' is needed: install.packages('FLCore', repos = 'https://flr-project.org/R')",
         call. = FALSE)
  }
}

.flindex <- function(s, units, startf, endf, desc) {
  .need_flcore()
  yrs <- seq(min(s$year), max(s$year))
  v <- s$index[match(yrs, s$year)]
  e <- s$se_log[match(yrs, s$year)]
  q <- FLCore::FLQuant(v, dimnames = list(age = "all", year = yrs), units = units)
  fi <- FLCore::FLIndexBiomass(index = q, name = s$name[1], desc = desc)
  FLCore::index.var(fi) <- FLCore::FLQuant(e, dimnames = list(age = "all", year = yrs),
                                           units = "se")
  fi@distribution <- "lognormal"
  if (!is.na(startf) || !is.na(endf)) {
    rg <- fi@range
    if (!is.na(startf)) rg["startf"] <- startf
    if (!is.na(endf)) rg["endf"] <- endf
    fi@range <- rg
  }
  fi
}

#' Indices as SPiCT input
#'
#' Returns the index elements of a SPiCT input list (`inp`, see
#' `spict::check.inp()`): `obsI`, `timeI` and `stdevfacI` (lists with one
#' vector per index) and a prior on the log observation standard deviation of
#' each index in `priors$logsdi`.
#'
#' The observation SD of each index (`sdI`) is the mean of its log-scale
#' standard errors. `stdevfacI` holds the standard errors relative to
#' `sdI`, so that SPiCT scales its estimated observation SD by year as the
#' standardisation does, and `priors$logsdi` is centred on `log(sdI)` with
#' standard deviation `sdi_sd`: `c(log(sdI), sdi_sd, 1)`, a list with one
#' prior per index when there are several. Merge the result into an `inp`
#' that holds the catch (`obsC`, `timeC`).
#'
#' @inheritParams index_set
#' @param ... Passed to [index_set()] (`model`, `years`, `normalise`).
#' @param timing Timing of each index within the year (fraction of the year,
#'   recycled), added to the year in `timeI`; e.g. `1/12` for the middle of
#'   January-February.
#' @param sdi_sd Standard deviation of the prior on `logsdi`; `NA` omits the
#'   prior.
#' @return A list with `obsI`, `timeI`, `stdevfacI`, `priors` (`logsdi`) and
#'   `sdI` (named vector, for reference).
#' @seealso [index_set()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' x <- add_targeting(select_species(x, keep = "Pagellus bogaraveo"))
#' fit <- fit_cpue(x, "Pagellus bogaraveo", models = "Cluster")
#' si <- as_spict(fit, timing = 0.5)
#' str(si)
#' \dontrun{
#' inp <- list(obsC = catch$catch, timeC = catch$year)
#' inp[names(si)[1:3]] <- si[1:3]
#' inp$priors$logsdi <- si$priors$logsdi
#' fit <- spict::fit.spict(inp)
#' }
as_spict <- function(x, ..., timing = 0, sdi_sd = 0.5) {
  s <- index_set(x, ...)
  nms <- unique(s$name)
  timing <- rep_len(timing, length(nms))
  sp <- split(s, factor(s$name, levels = nms))
  sdI <- vapply(sp, function(z) mean(z$se_log, na.rm = TRUE), numeric(1))
  out <- list(
    obsI = unname(lapply(sp, `[[`, "index")),
    timeI = unname(Map(function(z, t) z$year + t, sp, timing)),
    stdevfacI = unname(Map(function(z, m) {
      f <- z$se_log / m
      f[!is.finite(f)] <- 1
      f
    }, sp, sdI)))
  names(out$obsI) <- names(out$timeI) <- names(out$stdevfacI) <- nms
  pri <- lapply(sdI, function(m) c(log(m), sdi_sd, 1))
  out$priors <- if (is.na(sdi_sd)) list() else
    list(logsdi = if (length(nms) == 1) unname(pri[[1]]) else unname(pri))
  out$sdI <- sdI
  out
}

#' Indices as Stock Synthesis (SS3) CPUE input
#'
#' Returns the CPUE/abundance-index block of an SS3 data file (v3.30): one row
#' per observation with year, month, fleet, observation and log-scale standard
#' error, and a comment naming the index. Write it with [utils::write.csv()]
#' (e.g. `ss3idx.csv`) and paste the rows into the data file.
#'
#' Observations to keep in the file but exclude from the likelihood are
#' written with a negative year (SS3 convention), e.g. years with inconsistent
#' recording.
#'
#' @inheritParams index_set
#' @param ... Passed to [index_set()] (`model`, `years`, `normalise`).
#' @param fleet SS3 fleet/survey number of each index (default `1, 2, ...` in
#'   the order of the series; remap to the fleet numbers of the SS3 model).
#' @param month Month of each index (recycled). In SS3 v3.30 month 1.0 is the
#'   start of January and decimals are allowed: e.g. 2.0 for the middle of
#'   January-February, 1.5 for the middle of December-February (mid-January),
#'   7.0 for the middle of the year.
#' @param exclude Years written with a negative year (excluded in SS3).
#' @param names `"ss3"` (default): columns `#_yr`, `month`, `fleet`, `obs`,
#'   `stderr` as in the data file; `"r4ss"`: `year`, `month`, `index`, `obs`,
#'   `se_log` as in `r4ss::SS_readdat()` / `SS_writedat()`.
#' @param comment If `TRUE`, a last column `comment` with `#_` and the index
#'   name (and `excluded` for negative years), read as a comment by SS3.
#' @return A data frame with one row per observation.
#' @seealso [index_set()], [as_jabba()], [as_FLIndices()]
#' @export
#' @examples
#' sim <- simulate_logbook(seed = 1, plots = FALSE)
#' x <- cpue_data(sim$logbook, year = "year", species = "species",
#'                catch = "kg", id = "trip", vessel = "vessel")
#' x <- add_targeting(select_species(x, keep = "Pagellus bogaraveo"))
#' fit <- fit_cpue(x, "Pagellus bogaraveo", models = c("Base", "Cluster"))
#' ss3 <- as_ss3(fit, model = c("Base", "Cluster"), fleet = c(5, 6),
#'               month = 7, normalise = TRUE, exclude = 2000:2002)
#' head(ss3)
#' \dontrun{
#' write.csv(ss3, "ss3idx.csv", row.names = FALSE, quote = FALSE)
#' }
as_ss3 <- function(x, ..., fleet = NULL, month = 1, exclude = NULL,
                   names = c("ss3", "r4ss"), comment = TRUE) {
  names <- match.arg(names)
  s <- index_set(x, ...)
  nms <- unique(s$name)
  if (is.null(fleet)) fleet <- seq_along(nms)
  if (length(fleet) != length(nms)) {
    stop("`fleet` needs one number per series (", length(nms), ").", call. = FALSE)
  }
  month <- rep_len(month, length(nms))
  j <- match(s$name, nms)
  off <- s$year %in% exclude
  out <- data.frame(year = ifelse(off, -s$year, s$year), month = month[j],
                    fleet = as.integer(fleet[j]), obs = s$index, se_log = s$se_log)
  if (names == "ss3") {
    names(out) <- c("#_yr", "month", "fleet", "obs", "stderr")
  } else {
    names(out) <- c("year", "month", "index", "obs", "se_log")
  }
  if (comment) out$comment <- paste0("#_ ", s$name, ifelse(off, " excluded", ""))
  out
}
