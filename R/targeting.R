#' Principal components of the catch composition
#'
#' The catch of each operation is converted to proportions by weight (to
#' remove the effect of catch volume), transformed (fourth root by default,
#' which down-weights dominant species), and decomposed with a PCA on the
#' correlation matrix. The PC scores describe the catch composition of each
#' operation on a few continuous axes: the basis of both the DPC approach
#' (scores as smooth covariates) and of clustering into fishing tactics
#' (Winker et al. 2013, 2014).
#'
#' The number of non-trivial PCs is chosen with the Optimal Coordinates test
#' ([compute_oc()]).
#'
#' @param x A `cpue_data` object, ideally after [select_species()].
#' @param species Species to use (default: all in `x`).
#' @param transform Transformation of the proportions: `"fourth_root"`
#'   (default, as in the CPUE guidelines), `"sqrt"`, `"none"`, or
#'   `"presence"` (binary presence/absence instead of proportions).
#' @param scale Logical; PCA on the correlation matrix (default). The
#'   Kaiser-Guttman part of the OC test (eigenvalue >= 1) assumes this.
#' @param fit Optional logical vector (one value per operation): estimate the
#'   PCA (centring, scaling, rotation) and the OC test on these operations
#'   only and project all operations onto the axes with
#'   [stats::predict.prcomp()]. Used by [add_targeting()] with `fit_years`.
#' @return An object of class `cpue_pca` with `pca` (the [stats::prcomp()]
#'   result), `scores` (all operations), `eigen` (eigenvalues and shares of
#'   variance), `oc` (result of [compute_oc()]), `composition` (the
#'   transformed matrix, used for clustering), `species`, `labels`,
#'   `transform` and `fit` (operations used for estimation).
#' @export
composition_pca <- function(x, species = x$species,
                            transform = c("fourth_root", "sqrt", "none", "presence"),
                            scale = TRUE, fit = NULL) {
  check_cpue_data(x)
  transform <- match.arg(transform)
  species <- species_codes(x, species)
  comp <- composition_matrix(x$data, species, transform)
  if (is.null(fit)) fit <- rep(TRUE, nrow(comp))
  if (length(fit) != nrow(comp) || !is.logical(fit) || anyNA(fit)) {
    stop("`fit` must be a logical vector with one value per operation.", call. = FALSE)
  }
  if (sum(fit) < 3) stop("Fewer than three operations to fit the PCA.", call. = FALSE)
  constant <- apply(comp[fit, , drop = FALSE], 2, stats::sd) == 0
  if (any(constant)) {
    message("Dropped species with no variation in composition: ",
            paste(sp_label(x, species[constant]), collapse = ", "))
    comp <- comp[, !constant, drop = FALSE]
    species <- species[!constant]
  }
  if (ncol(comp) < 2) stop("At least two species are needed for a PCA.", call. = FALSE)
  pca <- stats::prcomp(comp[fit, , drop = FALSE], scale. = scale)
  scores <- if (all(fit)) pca$x else stats::predict(pca, comp)
  eig <- pca$sdev^2
  structure(list(pca = pca, scores = scores,
                 eigen = data.frame(PC = seq_along(eig), eigenvalue = eig,
                                    share = eig / sum(eig), cumulative = cumsum(eig) / sum(eig)),
                 oc = compute_oc(eig), composition = comp, species = species,
                 labels = sp_label(x, species), transform = transform, scale = scale,
                 fit = fit),
            class = "cpue_pca")
}

#' Transformed catch composition of operations
#' @noRd
composition_matrix <- function(data, species, transform) {
  m <- as.matrix(data[, species, drop = FALSE])
  tot <- rowSums(m)
  if (any(tot == 0)) {
    stop(sum(tot == 0), " operation(s) caught none of the species; their composition is ",
         "undefined. Use select_species() (which removes them) or subset().", call. = FALSE)
  }
  switch(transform,
         fourth_root = (m / tot)^0.25,
         sqrt = sqrt(m / tot),
         none = m / tot,
         presence = (m > 0) * 1)
}

#' @export
print.cpue_pca <- function(x, ...) {
  fitted <- if (is.null(x$fit)) nrow(x$scores) else sum(x$fit)
  cat("<cpue_pca>", nrow(x$scores), "operations", if (fitted < nrow(x$scores))
    sprintf("(PCA fitted on %d)", fitted), ",", length(x$species), "species,",
      "transform:", x$transform, "\n")
  e <- x$eigen
  e$eigenvalue <- round(e$eigenvalue, 3)
  e$share <- round(e$share, 3)
  e$cumulative <- round(e$cumulative, 3)
  print(utils::head(e, 8), row.names = FALSE)
  cat("Optimal Coordinates: ", x$oc$npc, " non-trivial PC(s); ",
      x$oc$n_clusters, " clusters by the npc + 1 rule\n", sep = "")
  invisible(x)
}

#' Optimal Coordinates test for the number of non-trivial components
#'
#' Non-graphical solution to Cattell's scree test (Raiche et al. 2013, as in
#' the R package nFactors), combined with the Kaiser-Guttman rule. For each
#' component \eqn{i}, the eigenvalue is predicted by the straight line through
#' eigenvalues \eqn{i + 1} and \eqn{k} (the last); component \eqn{i} is
#' retained if its eigenvalue is at least this prediction and at least 1. The
#' number of non-trivial components is the number retained in sequence from
#' the first. Winker et al. (2014) take the number of clusters as this number
#' plus one.
#'
#' Same rule as `computeOC()` in the CPUE guidelines, with the last
#' comparisons defined so that the test also works when every component
#' passes (at most \eqn{k - 1} components are retained).
#'
#' @param eig Eigenvalues in decreasing order.
#' @return A list with `npc`, `n_clusters` (`npc + 1`), `tests` (eigenvalue,
#'   OC prediction and whether each component passes) and `lines` (the OC
#'   lines for plotting).
#' @references Raiche, G., Walls, T.A., Magis, D., Riopel, M. and Blais, J.-G.
#'   (2013) Non-graphical solutions for Cattell's scree test. Methodology 9:
#'   23-29.
#' @export
#' @examples
#' compute_oc(c(2.1, 1.3, 0.4, 0.2))$npc
compute_oc <- function(eig) {
  nk <- length(eig)
  if (nk < 2) stop("Need at least two eigenvalues.", call. = FALSE)
  pred <- rep(NA_real_, nk)
  pass <- rep(FALSE, nk)
  lines <- list()
  for (i in seq_len(nk - 1)) {
    if (i + 1 == nk) {
      pred[i] <- eig[nk]
      slope <- 0
    } else {
      slope <- (eig[nk] - eig[i + 1]) / (nk - i - 1)
      pred[i] <- eig[i + 1] - slope
    }
    pass[i] <- eig[i] >= pred[i] && eig[i] >= 1
    lines[[i]] <- data.frame(tested = i, PC = i:nk, eigenvalue = pred[i] + slope * (0:(nk - i)))
  }
  npc <- which(!pass)[1] - 1
  if (npc < 1) {
    warning("No component passed the OC test; using 1.", call. = FALSE)
    npc <- 1
  }
  shown <- seq_len(min(npc + 1, nk - 1))
  list(npc = npc, n_clusters = npc + 1,
       tests = data.frame(PC = seq_len(nk), eigenvalue = eig, oc_prediction = pred,
                          retained = seq_len(nk) <= npc),
       lines = do.call(rbind, lines[shown]))
}

#' Plot the scree and the Optimal Coordinates test
#'
#' @param x A `cpue_data` object after [add_targeting()], or a `cpue_pca`.
#' @return A ggplot object.
#' @export
plot_oc <- function(x) {
  p <- get_pca(x)
  t <- p$oc$tests
  t$status <- ifelse(t$retained, "retained", "trivial")
  ggplot(t, aes(x = .data$PC, y = .data$eigenvalue)) +
    geom_line(data = p$oc$lines,
              aes(group = factor(.data$tested), colour = factor(.data$tested)),
              linetype = 2) +
    geom_hline(yintercept = 1, linetype = 3) +
    geom_line() +
    geom_point(aes(shape = .data$status), size = 3) +
    scale_shape_manual(values = c(retained = 16, trivial = 1)) +
    scale_x_continuous(breaks = t$PC) +
    labs(x = "Principal component", y = "Eigenvalue", shape = NULL,
         colour = "OC line for PC",
         title = sprintf("Optimal Coordinates: %d non-trivial PC(s)", p$oc$npc)) +
    theme_bw()
}

#' Cluster operations into fishing tactics
#'
#' @param pca Output of [composition_pca()].
#' @param n_clusters Number of clusters (default: OC result + 1, Winker et
#'   al. 2014).
#' @param method `"clara"` (default, as in the CPUE guidelines; needs the
#'   `cluster` package), `"kmeans"`, or `"hclust"` (Ward's method; only for
#'   up to 20,000 operations).
#' @param on Cluster the transformed composition (`"composition"`, default,
#'   as in the guidelines) or the first `n_pc` PC scores (`"pcs"`).
#' @param n_pc Number of PCs when `on = "pcs"`.
#' @param seed Random seed (CLARA and k-means are stochastic).
#' @param samples,sampsize CLARA settings: number of samples and sample size
#'   (default 250, as in the guidelines).
#' @param fit Optional logical vector: cluster these operations only and
#'   assign the others to the nearest cluster centre (default: `pca$fit`, i.e.
#'   the operations used to fit the PCA).
#' @return An integer vector of clusters, numbered by decreasing size (in the
#'   fitted operations), with attribute `centres`: the cluster centres in the
#'   clustering space (CLARA medoids, k-means centres or, for hclust, cluster
#'   means; row names are the cluster numbers), used to assign new
#'   operations to the nearest centre.
#' @export
cluster_tactics <- function(pca, n_clusters = pca$oc$n_clusters,
                            method = c("clara", "kmeans", "hclust"),
                            on = c("composition", "pcs"), n_pc = pca$oc$npc,
                            seed = 1, samples = 100, sampsize = 250, fit = pca$fit) {
  if (!inherits(pca, "cpue_pca")) stop("`pca` must come from composition_pca().", call. = FALSE)
  method <- match.arg(method)
  on <- match.arg(on)
  z_all <- if (on == "composition") pca$composition else pca$scores[, seq_len(n_pc), drop = FALSE]
  if (is.null(fit)) fit <- rep(TRUE, nrow(z_all))
  z <- z_all[fit, , drop = FALSE]
  n <- nrow(z)
  if (n_clusters < 2 || n_clusters >= n) stop("`n_clusters` must be between 2 and n - 1.",
                                              call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  res <- switch(method,
    clara = {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Install the `cluster` package for CLARA, or use method = \"kmeans\".",
             call. = FALSE)
      }
      r <- cluster::clara(z, n_clusters, samples = samples,
                          sampsize = min(n, max(sampsize, 40 + 2 * n_clusters)),
                          pamLike = TRUE, rngR = TRUE)
      list(cl = r$clustering, centres = r$medoids)
    },
    kmeans = {
      r <- stats::kmeans(z, centers = n_clusters, nstart = 20, iter.max = 50)
      list(cl = r$cluster, centres = r$centers)
    },
    hclust = {
      if (n > 20000) stop("hclust needs all pairwise distances; use \"clara\" or ",
                          "\"kmeans\" for more than 20,000 operations.", call. = FALSE)
      cl <- stats::cutree(stats::hclust(stats::dist(z), method = "ward.D2"), k = n_clusters)
      list(cl = cl, centres = rowsum(z, cl) / as.numeric(table(cl)))
    })
  # number clusters by decreasing size so that results are comparable
  tab <- table(res$cl)
  ord <- as.integer(names(tab))[order(tab, decreasing = TRUE)]
  # centres are kept in the clustering algorithm's order (which decides
  # ties between equidistant centres); row names give the renumbered cluster
  centres <- res$centres
  rownames(centres) <- match(seq_len(n_clusters), ord)
  out <- integer(nrow(z_all))
  out[fit] <- match(res$cl, ord)
  if (any(!fit)) out[!fit] <- nearest_centre(z_all[!fit, , drop = FALSE], centres)
  attr(out, "centres") <- centres
  out
}

#' Assign rows to the nearest centre (Euclidean); ties go to the first centre
#' in the stored order, as in CLARA. Returns the cluster labels (row names).
#' @noRd
nearest_centre <- function(z, centres) {
  d <- vapply(seq_len(nrow(centres)), function(k)
    rowSums((z - matrix(centres[k, ], nrow(z), ncol(z), byrow = TRUE))^2), numeric(nrow(z)))
  if (!is.matrix(d)) d <- matrix(d, nrow = nrow(z))
  d <- round(d, 12)
  lab <- if (is.null(rownames(centres))) seq_len(nrow(centres)) else as.integer(rownames(centres))
  lab[max.col(-d, ties.method = "first")]
}

#' Add targeting covariates: PC scores and fishing tactics
#'
#' One step that runs [composition_pca()], the Optimal Coordinates test and
#' [cluster_tactics()], and adds the results to the data:
#' * `PC1`, `PC2`, ...: the PC scores of the catch composition (all PCs are
#'   added; the DPC approach uses the first `n_pc`);
#' * `FT`: the fishing tactic (cluster) as a factor.
#'
#' These are the targeting covariates of the DPC and Cluster models. Check
#' them with [plot_oc()], [plot_biplot()], [tactic_table()] and
#' [plot_targeting()] before using them.
#'
#' **Fitting years and freezing the targeting.** By default the PCA, the OC
#' test and the clusters are estimated from all operations, so adding a year
#' of data changes the PC scores and clusters of all earlier years (and can
#' shift the whole index). Two options keep them fixed:
#' * `fit_years`: estimate the PCA (centring, scaling and rotation), the OC
#'   test and the cluster centres from the operations of these years only;
#'   all operations are then projected onto the PC axes
#'   ([stats::predict.prcomp()]) and assigned to the nearest cluster centre.
#'   Useful to exclude years with inconsistent recording from the estimation,
#'   or to fit on a reference period and predict new years.
#' * `reference`: reuse the targeting of an earlier run (a `cpue_data` object
#'   after `add_targeting()`, or its `targeting` element) without
#'   re-estimating anything: the composition of every operation is projected
#'   onto the reference PCA and assigned to the reference cluster centres.
#'   This freezes the targeting covariates when an index is updated with a
#'   new year. See also [predict_targeting()].
#'
#' Operations outside the fitting years are only as well described as their
#' composition resembles that of the fitting years; check the PC scores and
#' clusters of predicted years with [plot_targeting()].
#'
#' @param x A `cpue_data` object, ideally after [select_species()].
#' @param n_pc Number of PCs to use (default: Optimal Coordinates test).
#' @param n_clusters Number of clusters (default: `n_pc + 1`).
#' @param method,on,seed Passed to [cluster_tactics()].
#' @param transform,species Passed to [composition_pca()].
#' @param fit_years Optional years whose operations are used to estimate the
#'   PCA, the OC test and the clusters (default: all years).
#' @param reference Optional earlier targeting to reuse unchanged (see
#'   Details); other settings are then taken from it.
#' @param ... Further arguments for [cluster_tactics()].
#' @return The `cpue_data` object with the new columns and a `targeting`
#'   element (PCA, settings, cluster centres, fitting years and cluster
#'   composition).
#' @export
#' @examples
#' \dontrun{
#' x <- add_targeting(x)              # OC test decides n_pc and n_clusters
#' plot_oc(x)
#' plot_biplot(x)
#' tactic_table(x)
#' plot_targeting(x)$time
#' # estimate on 2005-2024, predict the PCs and tactics of all years
#' x <- add_targeting(x, fit_years = 2005:2024)
#' # freeze the targeting of last year's assessment for the update
#' x_new <- add_targeting(x_new, reference = x_last_year)
#' }
add_targeting <- function(x, n_pc = NULL, n_clusters = NULL,
                          method = c("clara", "kmeans", "hclust"),
                          on = c("composition", "pcs"), transform = "fourth_root",
                          species = x$species, seed = 1, fit_years = NULL,
                          reference = NULL, ...) {
  check_cpue_data(x)
  if (!is.null(reference)) return(predict_targeting(x, reference))
  method <- match.arg(method)
  on <- match.arg(on)
  yr <- x$data[[x$roles$year]]
  fit <- if (is.null(fit_years)) rep(TRUE, nrow(x$data)) else yr %in% fit_years
  if (!any(fit)) stop("No operations in `fit_years`.", call. = FALSE)
  pca <- composition_pca(x, species = species, transform = transform, fit = fit)
  if (is.null(n_pc)) n_pc <- pca$oc$npc
  if (n_pc < 1 || n_pc > ncol(pca$scores)) {
    stop("`n_pc` must be between 1 and ", ncol(pca$scores), call. = FALSE)
  }
  if (is.null(n_clusters)) n_clusters <- n_pc + 1
  ft <- cluster_tactics(pca, n_clusters = n_clusters, method = method, on = on,
                        n_pc = n_pc, seed = seed, fit = fit, ...)
  x <- set_targeting(x, pca, ft, n_pc, n_clusters)
  fy <- sort(unique(yr[fit]))
  x$targeting <- list(pca = pca, n_pc = n_pc, n_clusters = n_clusters, method = method,
                      on = on, transform = transform, oc_npc = pca$oc$npc,
                      centres = attr(ft, "centres"), fit_years = fy, frozen = FALSE,
                      assignment = targeting_assignment(x, pca$composition, ft))
  x$targeting$composition <- tactic_table(x)
  x$log <- rbind(x$log, data.frame(
    step = "add targeting", operations = nrow(x$data),
    note = sprintf("%d PC(s) used (OC test: %d), %d clusters (%s on %s)%s",
                   n_pc, pca$oc$npc, n_clusters, method, on,
                   if (all(fit)) "" else sprintf("; fitted on %d operations of %s",
                                                 sum(fit), year_range_text(fy)))))
  x
}

#' Apply an earlier targeting to new data
#'
#' Projects the catch composition of every operation of `x` onto the PCA of
#' an earlier [add_targeting()] run and assigns it to the nearest of the
#' earlier cluster centres, without re-estimating anything. Use it to freeze
#' the targeting covariates when an index is updated with new years, so that
#' the PC scores and tactics of earlier years do not change (no retrospective
#' shift from the targeting step). Equivalent to
#' `add_targeting(x, reference = reference)`.
#'
#' Operations of the reference (same id) whose composition is unchanged keep
#' their reference tactic exactly; new or revised operations are assigned to
#' the nearest reference cluster centre (ties broken as in the clustering
#' algorithm). Species of the reference that are missing in `x` are treated
#' as zero catches; species of `x` not in the reference are ignored. Operations that
#' caught none of the reference species cannot be projected and cause an
#' error (remove them with [subset()]).
#'
#' @param x A `cpue_data` object.
#' @param reference A `cpue_data` object after [add_targeting()], or its
#'   `targeting` element.
#' @return `x` with the PC scores, the tactics `FT` and a `targeting` element
#'   marked `frozen = TRUE`.
#' @seealso [add_targeting()]
#' @export
predict_targeting <- function(x, reference) {
  check_cpue_data(x)
  ref <- if (inherits(reference, "cpue_data")) reference$targeting else reference
  if (is.null(ref$pca) || is.null(ref$centres)) {
    stop("`reference` must come from add_targeting() of this package version (rerun it).",
         call. = FALSE)
  }
  sp <- ref$pca$species
  d <- x$data
  miss <- setdiff(sp, names(d))
  if (length(miss)) {
    message("Reference species missing in the data (treated as zero): ",
            paste(ref$pca$labels[match(miss, sp)], collapse = ", "))
    for (m in miss) d[[m]] <- 0
    x$data <- d
  }
  comp <- composition_matrix(d, sp, ref$transform)
  scores <- stats::predict(ref$pca$pca, comp)
  pca <- ref$pca
  pca$scores <- scores
  pca$composition <- comp
  pca$fit <- rep(FALSE, nrow(comp))
  z <- if (ref$on == "composition") comp else scores[, seq_len(ref$n_pc), drop = FALSE]
  ft <- nearest_centre(z, ref$centres)
  # operations of the reference with unchanged composition keep their tactic
  # (exact reproduction, including ties resolved inside the clustering)
  kept <- 0L
  a <- ref$assignment
  if (!is.null(a) && !is.null(x$roles$id)) {
    m <- match(as.character(x$data[[x$roles$id]]), a$id)
    same <- !is.na(m) & abs(composition_key(comp) - a$key[m]) < 1e-9
    ft[same] <- a$FT[m[same]]
    kept <- sum(same)
  }
  x <- set_targeting(x, pca, ft, ref$n_pc, ref$n_clusters)
  x$targeting <- ref[setdiff(names(ref), c("composition", "previous", "assignment"))]
  x$targeting$pca <- pca
  x$targeting$frozen <- TRUE
  x$targeting$assignment <- targeting_assignment(x, comp, ft)
  x$targeting$composition <- tactic_table(x)
  x$log <- rbind(x$log, data.frame(
    step = "predict targeting", operations = nrow(x$data),
    note = sprintf("frozen targeting from %s: %d PC(s), %d clusters; %d operations as in the reference, %d assigned to the nearest centre",
                   year_range_text(ref$fit_years), ref$n_pc, ref$n_clusters, kept,
                   nrow(x$data) - kept)))
  x
}

#' Operation ids, tactics and composition keys (to reproduce tactics exactly)
#' @noRd
targeting_assignment <- function(x, comp, ft) {
  if (is.null(x$roles$id)) return(NULL)
  data.frame(id = as.character(x$data[[x$roles$id]]), FT = as.integer(ft),
             key = composition_key(comp), stringsAsFactors = FALSE)
}

#' Numeric fingerprint of each composition row
#' @noRd
composition_key <- function(comp) {
  w <- 1 / sqrt(seq_len(ncol(comp)) + 0.5)
  as.vector(comp %*% w)
}

#' Add PC scores and tactics to the data
#' @noRd
set_targeting <- function(x, pca, ft, n_pc, n_clusters) {
  old_pc <- grep("^PC[0-9]+$", names(x$data), value = TRUE)
  x$data <- x$data[, setdiff(names(x$data), c(old_pc, "FT")), drop = FALSE]
  scores <- pca$scores
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  x$data <- cbind(x$data, as.data.frame(scores))
  x$data$FT <- factor(as.integer(ft), levels = seq_len(n_clusters))
  x
}

#' Text for a set of years, e.g. "2005-2024"
#' @noRd
year_range_text <- function(y) {
  if (is.null(y) || !length(y)) return("unknown years")
  y <- sort(unique(y))
  if (all(diff(y) == 1)) paste0(min(y), "-", max(y)) else paste(y, collapse = ", ")
}

#' Keep or drop fishing tactics and recompute the targeting covariates
#'
#' An intermediate cleaning step: remove operations of fishing tactics
#' (clusters) that are irrelevant for the species to be standardised, e.g.
#' the pelagic tactic when analysing a demersal species, and then recompute
#' the PCA and clusters on the remaining operations. Without the other
#' tactics the PCA describes the variation that matters for the target (e.g.
#' reef versus sediment, shallow versus deep) instead of spending its first
#' axes on the contrast with the removed tactics.
#'
#' Clusters can be chosen by number (`keep`, `drop`) or by their species
#' content: with `species` and `min_share`, clusters are kept if the given
#' species make up at least `min_share` of their catch (by weight), e.g. all
#' demersal species with `min_share = 0.5`.
#'
#' **Caution.** The clusters are derived from the catch composition, so this
#' selection is based on the data being standardised. Select on groups of
#' *other* species (or on the tactic that clearly does not catch the target)
#' rather than on the target alone, and compare the index with and without
#' the selection.
#'
#' @param x A `cpue_data` object after [add_targeting()].
#' @param keep,drop Cluster numbers to keep or drop.
#' @param species Species (column or original names) whose combined share of
#'   a cluster's catch decides whether it is kept.
#' @param min_share Minimum combined share of `species` (default 0.5).
#' @param retarget Logical; rerun [add_targeting()] on the remaining
#'   operations (default). If `FALSE`, the old PC scores and clusters are
#'   removed.
#' @param n_pc,n_clusters Passed to [add_targeting()] (default: the Optimal
#'   Coordinates test is run again on the subset).
#' @param pca_species Species used for the new PCA (default: all species of
#'   `x`; species without catches in the subset are dropped automatically).
#' @param ... Further arguments for [add_targeting()]; by default the
#'   clustering method, `on`, `transform` and `fit_years` of the previous call
#'   are reused.
#' @return The subsetted `cpue_data` object with new targeting covariates.
#'   `x$targeting$previous` holds the composition of the previous clusters
#'   and which were kept.
#' @seealso [tactic_table()], [select_vessels()], [stephens_maccall()]
#' @export
#' @examples
#' \dontrun{
#' tactic_table(x)                                   # which cluster is pelagic?
#' xd <- select_tactics(x, drop = 3)                 # by number
#' demersal <- c("Dentex dentex", "Pagellus bogaraveo", "Conger conger")
#' xd <- select_tactics(x, species = demersal, min_share = 0.5)
#' tactic_table(xd)
#' }
select_tactics <- function(x, keep = NULL, drop = NULL, species = NULL, min_share = 0.5,
                           retarget = TRUE, n_pc = NULL, n_clusters = NULL,
                           pca_species = x$species, ...) {
  check_targeting(x)
  levs <- levels(x$data$FT)
  ft <- as.character(x$data$FT)
  sel <- rep(TRUE, length(levs))
  why <- character(0)
  if (!is.null(keep)) {
    sel <- sel & levs %in% as.character(keep)
    why <- c(why, paste("keep", paste(keep, collapse = ", ")))
  }
  if (!is.null(drop)) {
    sel <- sel & !levs %in% as.character(drop)
    why <- c(why, paste("drop", paste(drop, collapse = ", ")))
  }
  share <- NULL
  if (!is.null(species)) {
    codes <- species_codes(x, species)
    if (anyNA(codes)) {
      stop("Species not found: ", paste(species[is.na(codes)], collapse = ", "), call. = FALSE)
    }
    m <- as.matrix(x$data[, x$species, drop = FALSE])
    agg <- rowsum(m, ft)[levs, , drop = FALSE]
    share <- rowSums(agg[, codes, drop = FALSE]) / rowSums(agg)
    sel <- sel & share >= min_share
    why <- c(why, sprintf("clusters with >= %.0f%% of catch from %d species",
                          100 * min_share, length(codes)))
  }
  if (!length(why)) stop("Give `keep`, `drop` or `species`.", call. = FALSE)
  if (!any(sel)) stop("No cluster is kept.", call. = FALSE)
  kept <- levs[sel]
  prev <- x$targeting
  previous <- list(composition = prev$composition, kept = kept, share = share)
  x$targeting <- NULL
  x <- filter_ops(x, ft %in% kept, "select tactics",
                  sprintf("clusters %s kept (%s)", paste(kept, collapse = ", "),
                          paste(why, collapse = "; ")))
  if (!retarget) {
    old <- c(grep("^PC[0-9]+$", names(x$data), value = TRUE), "FT")
    x$data <- x$data[, setdiff(names(x$data), old), drop = FALSE]
    return(x)
  }
  args <- list(...)
  for (a in c("method", "on", "transform", "fit_years")) {
    if (is.null(args[[a]])) args[[a]] <- prev[[a]]
  }
  x <- do.call(add_targeting, c(list(x, n_pc = n_pc, n_clusters = n_clusters,
                                     species = pca_species), args))
  x$targeting$previous <- previous
  x
}

#' Species composition of each fishing tactic
#'
#' @param x A `cpue_data` object after [add_targeting()].
#' @param type `"weight"` (percentage of the cluster's total catch, as in the
#'   CPUE guidelines) or `"occurrence"` (percentage of operations catching
#'   each species).
#' @return A data frame with one row per cluster: `FT`, number and share of
#'   operations, the dominant species, and one column per species.
#' @export
tactic_table <- function(x, type = c("weight", "occurrence")) {
  check_targeting(x)
  type <- match.arg(type)
  sp <- x$targeting$pca$species
  m <- as.matrix(x$data[, sp, drop = FALSE])
  if (type == "occurrence") m <- (m > 0) * 1
  agg <- rowsum(m, x$data$FT)
  pct <- if (type == "weight") 100 * agg / rowSums(agg) else 100 * agg / as.numeric(table(x$data$FT))
  n <- as.numeric(table(x$data$FT))
  out <- data.frame(FT = rownames(agg), operations = n, share = round(n / sum(n), 3),
                    dominant = sp_label(x, sp[apply(agg, 1, which.max)]))
  pct <- round(pct, 1)
  colnames(pct) <- sp_label(x, sp)
  cbind(out, as.data.frame(pct, check.names = FALSE), row.names = NULL)
}

#' PCA biplot of the catch composition
#'
#' Operations as points (coloured by fishing tactic, by default) and species
#' loadings as arrows, on two chosen PCs.
#'
#' @param x A `cpue_data` object after [add_targeting()] (or a `cpue_pca`,
#'   then points are not coloured).
#' @param pcs The two PCs to show.
#' @param colour_by Column of `x$data` used to colour points (default `"FT"`;
#'   `NULL` for none).
#' @param max_points Plot a random sample of at most this many operations.
#' @return A ggplot object.
#' @export
plot_biplot <- function(x, pcs = c(1, 2), colour_by = "FT", max_points = 5000) {
  p <- get_pca(x)
  sc <- as.data.frame(p$scores[, pcs, drop = FALSE])
  names(sc) <- c("a", "b")
  if (inherits(x, "cpue_data") && !is.null(colour_by)) {
    sc$col <- factor(x$data[[colour_by]])
  } else {
    sc$col <- factor("all")
    colour_by <- NULL
  }
  if (nrow(sc) > max_points) sc <- sc[sample.int(nrow(sc), max_points), ]
  ld <- as.data.frame(p$pca$rotation[, pcs, drop = FALSE])
  names(ld) <- c("a", "b")
  stretch <- 0.8 * max(abs(c(sc$a, sc$b))) / max(abs(c(ld$a, ld$b)))
  ld$a <- ld$a * stretch
  ld$b <- ld$b * stretch
  ld$label <- p$labels
  share <- round(100 * p$eigen$share[pcs], 1)
  g <- ggplot(sc, aes(x = .data$a, y = .data$b)) +
    geom_point(aes(colour = .data$col), alpha = 0.4, size = 0.9) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = .data$a, yend = .data$b),
                 arrow = grid::arrow(length = grid::unit(0.2, "cm")), colour = "grey20") +
    geom_label(data = ld, aes(label = .data$label), size = 3, alpha = 0.8) +
    labs(x = sprintf("PC%d (%s%%)", pcs[1], share[1]),
         y = sprintf("PC%d (%s%%)", pcs[2], share[2]), colour = colour_by) +
    theme_bw()
  if (is.null(colour_by)) g <- g + guides(colour = "none")
  g
}

#' What the targeting covariates capture
#'
#' @param x A `cpue_data` object after [add_targeting()].
#' @param by Optional column (e.g. `"month"`, `"HBF"`, `"port"`, `"vessel"`)
#'   to compare PC scores and tactics across its levels.
#' @param res Cell size for the maps (if positions are available).
#' @return A list of ggplots:
#'   * `loadings`: species loadings on the PCs used;
#'   * `composition`: species composition of each tactic;
#'   * `time`: share of tactics by year;
#'   * `pcs_time`: median and 10-90% range of the PC scores by year;
#'   * `by`: PC1 scores and tactic shares by the levels of `by` (if given);
#'   * `map_tactic`, `map_pc`: most frequent tactic and mean PC1 per cell (if
#'     positions are available).
#' @export
plot_targeting <- function(x, by = NULL, res = 1) {
  check_targeting(x)
  tg <- x$targeting
  r <- x$roles
  d <- x$data
  yr <- d[[r$year]]
  used <- paste0("PC", seq_len(tg$n_pc))
  cols <- cluster_colours(tg$n_clusters)
  out <- list()

  ld <- tg$pca$pca$rotation[, seq_len(tg$n_pc), drop = FALSE]
  ld <- data.frame(species = factor(rep(tg$pca$labels, ncol(ld)), levels = rev(tg$pca$labels)),
                   PC = rep(used, each = nrow(ld)), loading = as.vector(ld))
  out$loadings <- ggplot(ld, aes(x = .data$loading, y = .data$species,
                                 fill = .data$loading > 0)) +
    geom_col(show.legend = FALSE) + geom_vline(xintercept = 0) +
    facet_wrap(~PC) + scale_fill_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "#d95f0e")) +
    labs(x = "Loading", y = NULL) + theme_bw()

  tt <- tactic_table(x)
  sp_cols <- names(tt)[-(1:4)]
  comp <- data.frame(FT = rep(paste0(tt$FT, " (", tt$dominant, ")"), length(sp_cols)),
                     species = factor(rep(sp_cols, each = nrow(tt)), levels = sp_cols),
                     pct = unlist(tt[sp_cols], use.names = FALSE))
  out$composition <- ggplot(comp, aes(x = .data$FT, y = .data$pct, fill = .data$species)) +
    geom_col() + scale_fill_manual(values = species_colours(length(sp_cols))) +
    labs(x = "Fishing tactic (dominant species)", y = "% of catch", fill = "Species") +
    theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  out$time <- ggplot(data.frame(year = yr, FT = d$FT), aes(x = .data$year, fill = .data$FT)) +
    geom_bar(position = "fill") + scale_fill_manual(values = cols) +
    labs(x = "Year", y = "Share of operations", fill = "Tactic") + theme_bw()

  q <- do.call(rbind, lapply(used, function(pc) {
    s <- tapply(d[[pc]], yr, stats::quantile, probs = c(0.1, 0.5, 0.9))
    data.frame(PC = pc, year = as.numeric(names(s)), do.call(rbind, s))
  }))
  names(q)[3:5] <- c("lo", "mid", "hi")
  out$pcs_time <- ggplot(q, aes(x = .data$year, y = .data$mid)) +
    geom_ribbon(aes(ymin = .data$lo, ymax = .data$hi), fill = "grey80") +
    geom_line() + geom_point(size = 1) + geom_hline(yintercept = 0, linetype = 2) +
    facet_wrap(~PC, scales = "free_y") +
    labs(x = "Year", y = "PC score (median, 10-90%)") + theme_bw()

  if (!is.null(by)) {
    if (!by %in% names(d)) stop("Column `", by, "` not found.", call. = FALSE)
    bd <- data.frame(level = factor(d[[by]]), PC1 = d$PC1, FT = d$FT)
    out$by <- list(
      pc = ggplot(bd, aes(x = .data$level, y = .data$PC1)) +
        geom_boxplot(outlier.size = 0.3, fill = "grey85") +
        labs(x = by, y = "PC1 score") + theme_bw(),
      tactic = ggplot(bd, aes(x = .data$level, fill = .data$FT)) +
        geom_bar(position = "fill") + scale_fill_manual(values = cols) +
        labs(x = by, y = "Share of operations", fill = "Tactic") + theme_bw())
  }

  if (!is.null(r$lon) && !is.null(r$lat)) {
    cell <- data.frame(lon = floor(d[[r$lon]] / res) * res + res / 2,
                       lat = floor(d[[r$lat]] / res) * res + res / 2,
                       FT = d$FT, PC1 = d$PC1)
    key <- paste(cell$lon, cell$lat)
    tab <- table(key, cell$FT)
    mode_ft <- data.frame(key = rownames(tab),
                          FT = factor(colnames(tab)[apply(tab, 1, which.max)],
                                      levels = levels(d$FT)))
    xy <- cell[!duplicated(key), c("lon", "lat")]
    xy$key <- key[!duplicated(key)]
    mode_ft <- merge(xy, mode_ft, by = "key")
    pcm <- stats::aggregate(PC1 ~ lon + lat, cell, mean)
    out$map_tactic <- ggplot(mode_ft, aes(x = .data$lon, y = .data$lat, fill = .data$FT)) +
      geom_tile() + scale_fill_manual(values = cols) + coord_quickmap() +
      labs(x = r$lon, y = r$lat, fill = "Most frequent\ntactic") + theme_bw()
    out$map_pc <- ggplot(pcm, aes(x = .data$lon, y = .data$lat, fill = .data$PC1)) +
      geom_tile() + scale_fill_gradient2(low = "#d95f0e", mid = "white", high = "#2c7fb8") +
      coord_quickmap() + labs(x = r$lon, y = r$lat, fill = "Mean PC1") + theme_bw()
  }
  out
}

#' PCA from a cpue_data object or a cpue_pca
#' @noRd
get_pca <- function(x) {
  if (inherits(x, "cpue_pca")) return(x)
  if (inherits(x, "cpue_data") && !is.null(x$targeting)) return(x$targeting$pca)
  stop("Give a cpue_pca object or a cpue_data object after add_targeting().", call. = FALSE)
}

#' Check that add_targeting() has been run
#' @noRd
check_targeting <- function(x) {
  check_cpue_data(x)
  if (is.null(x$targeting) || is.null(x$data$FT)) {
    stop("Run add_targeting() first.", call. = FALSE)
  }
  if (nrow(x$targeting$pca$scores) != nrow(x$data)) {
    stop("The data were subset after add_targeting(); run add_targeting() again.",
         call. = FALSE)
  }
  invisible(x)
}
