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
#' @return An object of class `cpue_pca` with `pca` (the [stats::prcomp()]
#'   result), `scores`, `eigen` (eigenvalues and shares of variance), `oc`
#'   (result of [compute_oc()]), `composition` (the transformed matrix, used
#'   for clustering), `species`, `labels` and `transform`.
#' @export
composition_pca <- function(x, species = x$species,
                            transform = c("fourth_root", "sqrt", "none", "presence"),
                            scale = TRUE) {
  check_cpue_data(x)
  transform <- match.arg(transform)
  species <- species_codes(x, species)
  m <- as.matrix(x$data[, species, drop = FALSE])
  tot <- rowSums(m)
  if (any(tot == 0)) {
    stop(sum(tot == 0), " operation(s) caught none of the species; their composition is ",
         "undefined. Use select_species() (which removes them) or subset().", call. = FALSE)
  }
  comp <- switch(transform,
                 fourth_root = (m / tot)^0.25,
                 sqrt = sqrt(m / tot),
                 none = m / tot,
                 presence = (m > 0) * 1)
  constant <- apply(comp, 2, stats::sd) == 0
  if (any(constant)) {
    message("Dropped species with no variation in composition: ",
            paste(sp_label(x, species[constant]), collapse = ", "))
    comp <- comp[, !constant, drop = FALSE]
    species <- species[!constant]
  }
  if (ncol(comp) < 2) stop("At least two species are needed for a PCA.", call. = FALSE)
  pca <- stats::prcomp(comp, scale. = scale)
  eig <- pca$sdev^2
  structure(list(pca = pca, scores = pca$x,
                 eigen = data.frame(PC = seq_along(eig), eigenvalue = eig,
                                    share = eig / sum(eig), cumulative = cumsum(eig) / sum(eig)),
                 oc = compute_oc(eig), composition = comp, species = species,
                 labels = sp_label(x, species), transform = transform, scale = scale),
            class = "cpue_pca")
}

#' @export
print.cpue_pca <- function(x, ...) {
  cat("<cpue_pca>", nrow(x$scores), "operations,", length(x$species), "species,",
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
#' @return An integer vector of clusters, numbered by decreasing size.
#' @export
cluster_tactics <- function(pca, n_clusters = pca$oc$n_clusters,
                            method = c("clara", "kmeans", "hclust"),
                            on = c("composition", "pcs"), n_pc = pca$oc$npc,
                            seed = 1, samples = 100, sampsize = 250) {
  if (!inherits(pca, "cpue_pca")) stop("`pca` must come from composition_pca().", call. = FALSE)
  method <- match.arg(method)
  on <- match.arg(on)
  z <- if (on == "composition") pca$composition else pca$scores[, seq_len(n_pc), drop = FALSE]
  n <- nrow(z)
  if (n_clusters < 2 || n_clusters >= n) stop("`n_clusters` must be between 2 and n - 1.",
                                              call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  cl <- switch(method,
    clara = {
      if (!requireNamespace("cluster", quietly = TRUE)) {
        stop("Install the `cluster` package for CLARA, or use method = \"kmeans\".",
             call. = FALSE)
      }
      cluster::clara(z, n_clusters, samples = samples,
                     sampsize = min(n, max(sampsize, 40 + 2 * n_clusters)),
                     pamLike = TRUE, rngR = TRUE)$clustering
    },
    kmeans = stats::kmeans(z, centers = n_clusters, nstart = 20, iter.max = 50)$cluster,
    hclust = {
      if (n > 20000) stop("hclust needs all pairwise distances; use \"clara\" or ",
                          "\"kmeans\" for more than 20,000 operations.", call. = FALSE)
      stats::cutree(stats::hclust(stats::dist(z), method = "ward.D2"), k = n_clusters)
    })
  # number clusters by decreasing size so that results are comparable
  ord <- order(table(cl), decreasing = TRUE)
  match(cl, as.integer(names(table(cl)))[ord])
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
#' @param x A `cpue_data` object, ideally after [select_species()].
#' @param n_pc Number of PCs to use (default: Optimal Coordinates test).
#' @param n_clusters Number of clusters (default: `n_pc + 1`).
#' @param method,on,seed Passed to [cluster_tactics()].
#' @param transform,species Passed to [composition_pca()].
#' @param ... Further arguments for [cluster_tactics()].
#' @return The `cpue_data` object with the new columns and a `targeting`
#'   element (PCA, settings and cluster composition).
#' @export
#' @examples
#' \dontrun{
#' x <- add_targeting(x)              # OC test decides n_pc and n_clusters
#' plot_oc(x)
#' plot_biplot(x)
#' tactic_table(x)
#' plot_targeting(x)$time
#' }
add_targeting <- function(x, n_pc = NULL, n_clusters = NULL,
                          method = c("clara", "kmeans", "hclust"),
                          on = c("composition", "pcs"), transform = "fourth_root",
                          species = x$species, seed = 1, ...) {
  check_cpue_data(x)
  method <- match.arg(method)
  on <- match.arg(on)
  pca <- composition_pca(x, species = species, transform = transform)
  if (is.null(n_pc)) n_pc <- pca$oc$npc
  if (n_pc < 1 || n_pc > ncol(pca$scores)) {
    stop("`n_pc` must be between 1 and ", ncol(pca$scores), call. = FALSE)
  }
  if (is.null(n_clusters)) n_clusters <- n_pc + 1
  ft <- cluster_tactics(pca, n_clusters = n_clusters, method = method, on = on,
                        n_pc = n_pc, seed = seed, ...)
  old_pc <- grep("^PC[0-9]+$", names(x$data), value = TRUE)
  x$data <- x$data[, setdiff(names(x$data), c(old_pc, "FT")), drop = FALSE]
  scores <- pca$scores
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  x$data <- cbind(x$data, as.data.frame(scores))
  x$data$FT <- factor(ft, levels = seq_len(n_clusters))
  x$targeting <- list(pca = pca, n_pc = n_pc, n_clusters = n_clusters, method = method,
                      on = on, transform = transform, oc_npc = pca$oc$npc)
  x$targeting$composition <- tactic_table(x)
  x$log <- rbind(x$log, data.frame(
    step = "add targeting", operations = nrow(x$data),
    note = sprintf("%d PC(s) used (OC test: %d), %d clusters (%s on %s)",
                   n_pc, pca$oc$npc, n_clusters, method, on)))
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
