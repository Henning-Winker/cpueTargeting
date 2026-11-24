sim_x <- function() {
  sim <- simulate_logbook(seed = 1, plots = FALSE)
  x <- cpue_data(sim$logbook, year = "year", species = "species", catch = "kg",
                 id = "trip", vessel = "vessel")
  select_species(x, keep = "Pagellus bogaraveo")
}

test_that("fit_years = all years reproduces the default", {
  x <- sim_x()
  a <- add_targeting(x)
  b <- add_targeting(x, fit_years = unique(x$data$year))
  expect_equal(a$data$PC1, b$data$PC1)
  expect_identical(a$data$FT, b$data$FT)
  expect_false(a$targeting$frozen)
  expect_equal(a$targeting$fit_years, sort(unique(x$data$year)))
})

test_that("fit_years estimates on a subset and predicts all operations", {
  x <- sim_x()
  last <- max(x$data$year)
  k <- x$data$year < last
  old <- add_targeting(subset(x, year < last))
  xf <- add_targeting(x, fit_years = sort(unique(x$data$year))[-length(unique(x$data$year))])
  expect_equal(xf$targeting$pca$pca$rotation, old$targeting$pca$pca$rotation)
  expect_equal(xf$data$PC1[k], old$data$PC1)
  expect_identical(as.integer(xf$data$FT[k]), as.integer(old$data$FT))
  # new year: scores from the old axes, tactics from the nearest centre
  comp_new <- xf$targeting$pca$composition[!k, , drop = FALSE]
  expect_equal(unname(xf$targeting$pca$scores[!k, 1]),
               unname(stats::predict(old$targeting$pca$pca, comp_new)[, 1]))
  expect_true(all(!is.na(xf$data$FT[!k])))
  expect_equal(sum(xf$targeting$pca$fit), sum(k))
  expect_match(tail(xf$log$note, 1), "fitted on")
  expect_error(add_targeting(x, fit_years = 1800), "No operations")
})

test_that("predict_targeting() freezes the targeting exactly for earlier years", {
  x <- sim_x()
  last <- max(x$data$year)
  for (m in c("clara", "kmeans")) {
    old <- add_targeting(subset(x, year < last), method = m)
    xn <- predict_targeting(x, old)
    k <- x$data$year < last
    expect_equal(xn$data$PC1[k], old$data$PC1)
    expect_identical(as.integer(xn$data$FT[k]), as.integer(old$data$FT))
    expect_true(xn$targeting$frozen)
    expect_equal(xn$targeting$n_clusters, old$targeting$n_clusters)
  }
  # same via add_targeting(reference = )
  xr <- add_targeting(x, reference = old$targeting)
  expect_identical(xr$data$FT, xn$data$FT)
  expect_match(tail(xn$log$step, 1), "predict targeting")
  expect_match(tail(xn$log$note, 1), "as in the reference")
  # re-estimating instead changes earlier years
  xa <- add_targeting(x, method = "kmeans")
  expect_false(isTRUE(all.equal(xa$data$PC1[x$data$year < last], old$data$PC1)))
})

test_that("revised operations are reassigned and missing species are zero", {
  x <- sim_x()
  old <- add_targeting(x)
  x2 <- x
  sp <- x2$species
  # make one operation catch only the target: its tactic is recomputed
  j <- 1
  x2$data[j, sp] <- 0
  x2$data[j, "Pagellus.bogaraveo"] <- 10
  xn <- predict_targeting(x2, old)
  expect_equal(unname(xn$targeting$pca$composition[j, "Pagellus.bogaraveo"]), 1)
  expect_identical(xn$data$FT[-j], old$data$FT[-j])
  # a reference species missing in the new data is treated as zero
  drop_sp <- setdiff(sp, "Pagellus.bogaraveo")[1]
  x3 <- x
  x3$data[[drop_sp]] <- NULL
  x3$species <- setdiff(x3$species, drop_sp)
  x3 <- subset(x3, rowSums(as.matrix(x3$data[, x3$species])) > 0)
  expect_message(predict_targeting(x3, old), "treated as zero")
  expect_error(predict_targeting(x, list(a = 1)), "must come from add_targeting")
})

test_that("frozen targeting works downstream", {
  x <- sim_x()
  last <- max(x$data$year)
  old <- add_targeting(subset(x, year < last))
  xn <- predict_targeting(x, old)
  f <- fit_cpue(xn, "Pagellus bogaraveo", models = c("DPC", "Cluster"))
  expect_true(last %in% f$indices$year)
  expect_s3_class(plot_targeting(xn)$time, "ggplot")
  expect_true(is.data.frame(tactic_table(xn)))
  # select_tactics() keeps the fitting years
  xf <- add_targeting(x, fit_years = 2000:2010)
  xs <- select_tactics(xf, drop = 4)
  expect_equal(xs$targeting$fit_years, sort(unique(xs$data$year[xs$data$year %in% 2000:2010])))
})
