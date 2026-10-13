# Synthetic logbook in long format: three tactics with distinct compositions,
# a rare species, a synonym, and tactic shares that change over time.
make_logbook <- function(n = 600, seed = 1) {
  set.seed(seed)
  year <- sort(sample(2001:2010, n, replace = TRUE))
  p_tactic <- cbind(seq(0.7, 0.1, length.out = 10), 0.2, seq(0.1, 0.7, length.out = 10))
  tactic <- vapply(year - 2000, function(y) sample(1:3, 1, prob = p_tactic[y, ]), 1L)
  mu <- rbind(c(50, 5, 0.3, 1), c(10, 30, 10, 1), c(0.3, 5, 40, 1))  # tactic x species
  species <- c("Pagellus bogaraveo", "Conger conger", "Lepidopus caudatus", "Rara avis")
  catch <- t(vapply(tactic, function(k) stats::rpois(4, mu[k, ] * stats::rgamma(4, 2, 2)),
                    numeric(4)))
  catch[, 4] <- catch[, 4] * (stats::runif(n) < 0.03)   # rare species
  ops <- data.frame(trip = seq_len(n), year = year,
                    vessel = sample(sprintf("V%02d", 1:15), n, replace = TRUE),
                    hours = round(stats::runif(n, 2, 10), 1),
                    month = sample(1:12, n, replace = TRUE),
                    lon = stats::runif(n, -6, -2), lat = stats::runif(n, 35, 37),
                    tactic = tactic)
  long <- ops[rep(seq_len(n), each = 4), ]
  long$species <- rep(species, n)
  long$kg <- as.vector(t(catch))
  long <- long[long$kg > 0, ]
  # a synonym recorded in some trips
  syn <- long$species == "Lepidopus caudatus" & stats::runif(nrow(long)) < 0.2
  long$species[syn] <- "Trichiurus lepturus"
  rownames(long) <- NULL
  list(long = long, ops = ops, catch = catch, species = species)
}

test_that("cpue_data() builds one row per operation from long data", {
  lb <- make_logbook()
  x <- cpue_data(lb$long, year = "year", species = "species", catch = "kg",
                 effort = "hours", id = "trip", vessel = "vessel",
                 lon = "lon", lat = "lat", covariates = "month",
                 synonyms = c("Trichiurus lepturus" = "Lepidopus caudatus"))
  expect_s3_class(x, "cpue_data")
  expect_equal(anyDuplicated(x$data$trip), 0L)
  expect_false("Trichiurus.lepturus" %in% x$species)
  # catches match the generating matrix (trips without any catch are absent)
  ops_kept <- match(x$data$trip, lb$ops$trip)
  expect_equal(x$data$Pagellus.bogaraveo, lb$catch[ops_kept, 1])
  expect_equal(x$data$Lepidopus.caudatus, lb$catch[ops_kept, 3])
  expect_equal(unname(x$species_labels["Pagellus.bogaraveo"]), "Pagellus bogaraveo")
  # species ordered by total catch
  expect_true(all(diff(colSums(x$data[, x$species])) <= 0))
  expect_equal(cpue_matrix(x, "Conger.conger")[, 1], x$data$Conger.conger / x$data$hours)
})

test_that("wide input gives the same result as long input", {
  lb <- make_logbook()
  long <- cpue_data(lb$long, year = "year", species = "species", catch = "kg",
                    id = "trip", synonyms = c("Trichiurus lepturus" = "Lepidopus caudatus"))
  wide_df <- data.frame(lb$ops[, c("trip", "year")], lb$catch)
  names(wide_df)[3:6] <- lb$species
  wide_df <- wide_df[rowSums(lb$catch) > 0, ]
  wide <- cpue_data(wide_df, year = "year", species = lb$species, id = "trip")
  expect_equal(sort(wide$species), sort(long$species))
  o <- match(long$data$trip, wide$data$trip)
  for (s in long$species) expect_equal(wide$data[[s]][o], long$data[[s]])
})

test_that("missing effort is removed and logged", {
  d <- data.frame(id = 1:4, year = 2001, A = c(1, 0, 2, 3), B = c(0, 1, 1, 0),
                  effort = c(1, NA, 0, 2))
  x <- cpue_data(d, year = "year", species = c("A", "B"), effort = "effort", id = "id")
  expect_equal(nrow(x$data), 2)
  expect_match(paste(x$log$step, collapse = " "), "effort")
})

test_that("species curve and selection follow the threshold", {
  lb <- make_logbook()
  x <- cpue_data(lb$long, year = "year", species = "species", catch = "kg", id = "trip",
                 synonyms = c("Trichiurus lepturus" = "Lepidopus caudatus"))
  sc <- species_curve(x)
  expect_equal(sum(sc$share), 1)
  expect_equal(sc$cumulative[nrow(sc)], 1)
  expect_equal(sc$label[nrow(sc)], "Rara avis")
  xs <- select_species(x, threshold = 0.95)
  expect_false("Rara.avis" %in% xs$species)
  xk <- select_species(x, threshold = 0.95, keep = "Rara avis")
  expect_true("Rara.avis" %in% xk$species)
  xo <- select_species(x, n = 2, pool_other = TRUE)
  expect_true("Other" %in% xo$species)
  expect_equal(length(xo$species), 3)
  expect_s3_class(plot_species_curve(x), "ggplot")
})

test_that("subset() and select_vessels() filter and log", {
  lb <- make_logbook()
  x <- cpue_data(lb$long, year = "year", species = "species", catch = "kg",
                 id = "trip", vessel = "vessel")
  x2 <- subset(x, year >= 2005)
  expect_true(all(x2$data$year >= 2005))
  expect_equal(nrow(x2$log), nrow(x$log) + 1)
  x3 <- select_vessels(x, min_years = 5, min_records = 3)
  expect_lte(nrow(x3$data), nrow(x$data))
  expect_error(subset(x, c(TRUE, FALSE)), "one TRUE/FALSE")
})

test_that("exploration plots are produced", {
  lb <- make_logbook()
  x <- cpue_data(lb$long, year = "year", species = "species", catch = "kg",
                 effort = "hours", id = "trip", vessel = "vessel", lon = "lon",
                 lat = "lat", covariates = "month")
  fl <- plot_fleet(x, by = "month")
  expect_true(all(c("records", "vessels", "vessel_years", "by_covariate") %in% names(fl)))
  nom <- nominal_cpue(x, "Pagellus bogaraveo")
  expect_equal(nrow(nom), 10)
  expect_true(all(nom$positive >= 0 & nom$positive <= 1))
  expect_s3_class(plot_nominal(x, c("Pagellus bogaraveo", "Conger conger"))$index, "ggplot")
  expect_s3_class(plot_composition(x), "ggplot")
  expect_s3_class(plot_cpue_map(x, species = "Conger conger", res = 0.5, by_period = 5),
                  "ggplot")
})

test_that("compute_oc() reproduces the guideline function", {
  # computeOC() of the CPUE guidelines, valid when not all components pass
  computeOC_guidelines <- function(eig) {
    nk <- length(eig); k <- 1:nk; aparallel <- rep(1, nk)
    cond1 <- TRUE; cond2 <- TRUE; i <- 0
    while (cond1 && cond2 && i < nk) {
      i <- i + 1
      ind <- k[c(i + 1, nk)]
      vp.p <- lm(eig[c(i + 1, nk)] ~ ind)
      vp.prec <- sum(c(1, i) * coef(vp.p))
      cond1 <- eig[i] >= vp.prec
      cond2 <- eig[i] >= aparallel[i]
      nc <- i - 1
    }
    nc
  }
  for (eig in list(c(2.1, 1.3, 0.4, 0.2), c(3.0, 0.6, 0.3, 0.1),
                   c(2.5, 1.6, 1.2, 0.4, 0.2, 0.1), c(1.8, 1.1, 0.9, 0.2))) {
    expect_equal(compute_oc(eig)$npc, computeOC_guidelines(eig))
  }
  # edge case where the guideline loop fails: all components pass
  expect_equal(compute_oc(c(2, 1.5))$npc, 1)
  expect_equal(compute_oc(c(2, 1.5))$n_clusters, 2)
})

test_that("Stephens & MacCall keeps operations in the target habitat", {
  lb <- make_logbook(n = 1500)
  x <- cpue_data(lb$long, year = "year", species = "species", catch = "kg", id = "trip",
                 covariates = "tactic",
                 synonyms = c("Trichiurus lepturus" = "Lepidopus caudatus"))
  x <- select_species(x, threshold = 0.99)
  # the target must be absent from some operations for the method to apply
  x$data$Rare <- x$data$Lepidopus.caudatus * (x$data$tactic == 3)
  x$species <- c(x$species, "Rare")
  x$species_labels <- c(x$species_labels, Rare = "Rare")
  sm <- stephens_maccall(x, target = "Rare", species = c("Pagellus.bogaraveo", "Conger.conger"))
  expect_s3_class(sm, "stephens_maccall")
  expect_true(sm$critical > 0 && sm$critical < 1)
  # kept operations are mostly from tactic 3
  expect_gt(mean(x$data$tactic[sm$selected] == 3), mean(x$data$tactic == 3))
  pl <- plot_stephens_maccall(sm)
  expect_true(all(vapply(pl, inherits, logical(1), "ggplot")))
  xs <- subset(x, sm$selected, note = "Stephens & MacCall")
  expect_equal(nrow(xs$data), sum(sm$selected))
})

test_that("add_targeting() recovers the tactics", {
  skip_if_not_installed("cluster")
  lb <- make_logbook(n = 1500)
  x <- cpue_data(lb$long, year = "year", species = "species", catch = "kg", id = "trip",
                 vessel = "vessel", lon = "lon", lat = "lat", covariates = c("month", "tactic"),
                 synonyms = c("Trichiurus lepturus" = "Lepidopus caudatus"))
  x <- select_species(x, threshold = 0.99)
  x <- add_targeting(x, n_clusters = 3)
  expect_true(all(c("PC1", "PC2", "FT") %in% names(x$data)))
  expect_equal(nlevels(x$data$FT), 3)
  # clusters agree with the true tactics (up to labels)
  purity <- sum(apply(table(x$data$FT, x$data$tactic), 1, max)) / nrow(x$data)
  expect_gt(purity, 0.8)
  # clusters are numbered by size
  expect_true(all(diff(as.numeric(table(x$data$FT))) <= 0))
  tt <- tactic_table(x)
  expect_equal(nrow(tt), 3)
  expect_equal(unname(rowSums(tt[, -(1:4)])), rep(100, 3), tolerance = 0.2)
  for (m in c("kmeans", "hclust")) {
    expect_equal(nlevels(add_targeting(x, n_clusters = 3, method = m)$data$FT), 3)
  }
  expect_s3_class(plot_oc(x), "ggplot")
  expect_s3_class(plot_biplot(x), "ggplot")
  pt <- plot_targeting(x, by = "month", res = 0.5)
  expect_true(all(c("loadings", "composition", "time", "pcs_time", "by",
                    "map_tactic", "map_pc") %in% names(pt)))
  # subsetting after add_targeting() is flagged
  expect_message(xs <- subset(x, year > 2003), "rerun add_targeting")
  expect_error(tactic_table(xs), "again")
})

test_that("the workflow runs on simulated data in long format", {
  skip_on_cran()
  skip_if_not_installed("cluster")
  sim <- simulate_bathymetric(seed = 1)
  logbook <- sim$Data[sim$Data$CPUE > 0, ]
  x <- cpue_data(logbook, year = "year", species = "Sp", catch = "CPUE", id = "SetID",
                 vessel = "vessel", lon = "Xobs", lat = "Yobs")
  expect_equal(length(x$species), 4)
  x <- add_targeting(x)
  expect_gte(x$targeting$n_pc, 1)
  expect_equal(x$targeting$n_clusters, x$targeting$n_pc + 1)
})
