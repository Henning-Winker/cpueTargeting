small_bathymetric <- function() {
  modifyList(bathymetric_settings(), list(n_years = 5, Range_X = 1:15, Range_Y = 1:15,
                                          Tot_effort = 100))
}

small_tuna <- function() {
  s <- tuna_settings()
  modifyList(s, list(n_years = 6, nx = 15, ny = 15, effort = 200,
                     price = s$price[c(1, 1, 11, 11, 21, 21), ]))
}

test_that("bathymetric model returns consistent output", {
  sim <- simulate_bathymetric(small_bathymetric(), seed = 1, plots = FALSE)
  expect_named(sim, c("Data", "Biomass", "Biomass.XY", "Env", "settings"))
  expect_true(all(c("year", "SetID", "X", "Y", "Xobs", "Yobs", "Sp", "CPUE",
                    "depth", "vessel") %in% names(sim$Data)))
  expect_equal(nrow(sim$Data), 4 * length(unique(sim$Data$SetID)))
  expect_false(anyNA(sim$Data$CPUE))
  expect_true(all(sim$Data$CPUE >= 0))
  expect_equal(nrow(sim$Biomass), 5 * 4)
  expect_equal(as.vector(tapply(sim$Biomass$True, sim$Biomass$Sp, mean)), rep(1, 4))
})

test_that("simulations are reproducible with a seed", {
  a <- simulate_tuna(small_tuna(), seed = 9, plots = FALSE)
  b <- simulate_tuna(small_tuna(), seed = 9, plots = FALSE)
  expect_identical(a$Data, b$Data)
})

test_that("tuna model returns HBF, latent overlap and sensible values", {
  sim <- simulate_tuna(small_tuna(), seed = 2, plots = FALSE)
  d <- sim$Data
  expect_true(all(c("HBF", "sst", "thermocline", "hook_depth",
                    "vertical_avail", "hbs_avail") %in% names(d)))
  expect_true(all(d$HBF %in% tuna_settings()$HBF_levels))
  expect_true(all(d$vertical_avail >= 0))
  expect_false(anyNA(d$CPUE))
  expect_equal(levels(sim$Biomass$Sp), c("YFT", "BET", "BFT", "ALB"))
  # deeper gear sets fish deeper
  depth_by_hbf <- tapply(d$hook_depth, d$HBF, mean)
  expect_true(all(diff(depth_by_hbf) > 0))
  sets <- as_set_data(d)
  expect_equal(nrow(sets), length(unique(d$SetID)))
  expect_true(all(c("YFT", "BET", "BFT", "ALB") %in% names(sets)))
})

test_that("plots are built", {
  sim <- simulate_tuna(small_tuna(), seed = 3)
  expect_true(all(c("Env", "Thermocline", "SpDens", "Hist", "CPUE", "HBF")
                  %in% names(sim$Plots)))
  expect_s3_class(sim$Plots$CPUE, "ggplot")
  expect_s3_class(plot_tuna_vertical_distributions(), "ggplot")
})

test_that("standardisation runs and returns indices for all models", {
  skip_if_not_installed("mgcv")
  sim <- simulate_tuna(small_tuna(), seed = 4, plots = FALSE)
  res <- standardize_tuna_indices(sim, models = c("Base+S", "HBF+S", "HBS+S"),
    add_models = c(HBF_period = "year + s(Xobs, Yobs, k = 8) + HBF_period"))
  expect_equal(nrow(res$indices), 4 * 4 * 6)
  expect_false(anyNA(res$indices$index))
  expect_s3_class(plot_indices(res), "ggplot")
  expect_s3_class(plot_indices(list(A = res, B = res), ci = TRUE), "ggplot")
  expect_s3_class(plot_indices(res, type = "smooth"), "ggplot")
  expect_s3_class(plot_indices(res, type = "points+smooth", span = 0.8), "ggplot")
  expect_s3_class(plot_diagnostics(list(A = res, B = res)), "ggplot")
  expect_true(all(c("Base+S", "HBF+S", "HBS+S", "HBF_period") %in% levels(res$indices$model)))
  fast <- standardize_tuna_indices(sim, models = c("Cluster", "DPC_HBF",
                                                  "Cluster+S", "DPC_HBF+S"),
                                   spatial_res = 1, engine = "bam")
  expect_false(anyNA(fast$indices$index))
  expect_equal(levels(fast$indices$model),
               c("Cluster", "DPC_HBF", "Cluster+S", "DPC_HBF+S"))
  expect_false(grepl("Xs", fast$formulas[["Cluster"]]))
  expect_true(grepl("Xs", fast$formulas[["Cluster+S"]]))
  expect_error(standardize_tuna_indices(sim, models = "Foo+S"), "Unknown")
  expect_equal(length(unique(paste(fast$data$Xs, fast$data$Ys))),
               length(unique(paste(fast$data$X, fast$data$Y))))
})

test_that("thermocline is fixed by default and moves with annual settings", {
  fixed <- simulate_tuna(small_tuna(), seed = 5, plots = FALSE)
  by_year <- tapply(fixed$Thermocline$anomaly, fixed$Thermocline$year, mean)
  expect_equal(as.vector(by_year), rep(by_year[1], 6), ignore_attr = TRUE)

  moving_cfg <- modifyList(small_tuna(), list(thermocline_annual_sd = 30,
                                              thermocline_annual_field_sd = 10))
  moving <- simulate_tuna(moving_cfg, seed = 5, plots = FALSE)
  expect_gt(sd(tapply(moving$Thermocline$anomaly, moving$Thermocline$year, mean)), 1)
  # same seed: same SST field and same fixed thermocline pattern
  expect_equal(moving$Env$SST, fixed$Env$SST)
  expect_equal(moving$Env$thermocline, fixed$Env$thermocline)
  # thermocline in Data is the anomaly of that year
  d <- moving$Data[moving$Data$Species == 1, ]
  th <- moving$Thermocline
  key <- match(paste(d$year, d$X, d$Y), paste(th$year, th$x, th$y))
  expect_equal(d$thermocline, th$anomaly[key])

  constant <- simulate_tuna(modifyList(small_tuna(), list(thermocline_sd = 0)),
                            seed = 5, plots = FALSE)
  expect_true(all(constant$Thermocline$anomaly == 0))
  expect_equal(constant$Env$SST, fixed$Env$SST)
})

test_that("effect plots are built", {
  moving_cfg <- modifyList(small_tuna(), list(thermocline_annual_sd = 30))
  sim <- simulate_tuna(moving_cfg, seed = 6, plots = FALSE)
  th <- plot_thermocline(sim)
  expect_s3_class(th$series, "ggplot")
  expect_s3_class(th$maps, "ggplot")
  expect_s3_class(plot_thermocline_effect(), "ggplot")
  expect_s3_class(plot_availability(sim), "ggplot")
})

test_that("fits at different resolutions can be combined", {
  skip_if_not_installed("mgcv")
  sim <- simulate_tuna(small_tuna(), seed = 7, plots = FALSE)
  a <- standardize_tuna_indices(sim, models = c("Base", "HBF+S"), engine = "bam")
  b <- standardize_tuna_indices(sim, models = "Base+S", spatial_res = 1,
                                engine = "bam", label = " (1x1)")
  ab <- bind_comparisons(a, b)
  expect_equal(levels(ab$indices$model), c("Base", "HBF+S", "Base+S (1x1)"))
  expect_equal(nrow(ab$diagnostics), 4 * 3)
  expect_error(bind_comparisons(a, a), "not unique")
  pc3 <- standardize_tuna_indices(sim, models = "DPC+S", n_pc = 3, engine = "bam")
  expect_true(grepl("s(PC3, k = 3)", pc3$formulas[["DPC+S"]], fixed = TRUE))
  expect_error(standardize_tuna_indices(sim, models = "DPC", n_pc = 5), "n_pc")
  pcp <- plot_pcs(pc3, pcs = c("PC1", "PC2", "PC3"), min_sets = 1)
  expect_named(pcp, c("loadings", "space", "depth", "time"))
  for (p in pcp) expect_s3_class(p, "ggplot")
  expect_equal(dim(pc3$pca$loadings), c(4, 4))
  expect_s3_class(plot_diagnostics(ab), "ggplot")
})

test_that("select_models keeps a subset that can be re-attached", {
  skip_if_not_installed("mgcv")
  sim <- simulate_tuna(small_tuna(), seed = 8, plots = FALSE)
  a <- standardize_tuna_indices(sim, models = c("Base", "Base+S"), engine = "bam")
  b <- standardize_tuna_indices(sim, models = "Base+S", spatial_res = 1,
                                engine = "bam", label = " (1x1)")
  control <- select_models(bind_comparisons(a, b), "Base+S (1x1)")
  expect_equal(levels(control$indices$model), "Base+S (1x1)")
  expect_equal(nrow(control$diagnostics), 4)
  expect_equal(levels(bind_comparisons(a, control)$indices$model),
               c("Base", "Base+S", "Base+S (1x1)"))
  expect_error(select_models(a, "HBF+S"), "not found")
})

test_that("3D distribution plots are built", {
  sim <- simulate_tuna(small_tuna(), seed = 9, plots = FALSE)
  p <- plot_sections(sim, year = 2)
  expect_s3_class(p$section, "ggplot")
  expect_s3_class(p$layers, "ggplot")
})

test_that("replicates run, resume from file and plot", {
  skip_if_not_installed("mgcv")
  f <- tempfile(fileext = ".rds")
  sc <- list(A = small_tuna())
  r1 <- run_replicates(sc, seeds = 1, models = c("Base", "DPC+S"), engine = "bam",
                       file = f,
                       extra = list(list(models = "Base+S", spatial_res = 1,
                                         engine = "bam", label = " (1x1)")))
  expect_equal(sort(unique(r1$diagnostics$seed)), 1)
  expect_true("Base+S (1x1)" %in% levels(r1$diagnostics$model))
  r2 <- run_replicates(sc, seeds = 1:2, models = c("Base", "DPC+S"), engine = "bam",
                       file = f,
                       extra = list(list(models = "Base+S", spatial_res = 1,
                                         engine = "bam", label = " (1x1)")))
  expect_equal(sort(unique(r2$diagnostics$seed)), 1:2)
  expect_equal(nrow(r2$diagnostics), 2 * 3 * 4)
  expect_s3_class(plot_replicates(r2), "ggplot")
})

test_that("vertical distribution plot works with and without overlay", {
  expect_s3_class(plot_tuna_vertical_distributions(overlay = FALSE), "ggplot")
  expect_s3_class(plot_tuna_vertical_distributions(overlay = TRUE), "ggplot")
})

test_that("settings presets differ only in targeting strength", {
  d <- tuna_settings()
  x <- tuna_settings("stress_test")
  expect_equal(max(d$price) / min(d$price), 2)
  expect_equal(max(x$price) / min(x$price), 4)
  expect_equal(x$preference, 1.5)
  expect_equal(x$HBF_choice_temp, 0.15)
  same <- setdiff(names(d), c("price", "preference", "HBF_choice_temp"))
  expect_identical(d[same], x[same])
  expect_error(tuna_settings("foo"))
})

test_that("B_start sets initial depletion", {
  s <- tuna_settings()
  expect_length(s$B_start, length(s$species))
  s$B_start <- c(1, 1, 0, 1)
  expect_error(simulate_tuna(s, seed = 1, plots = FALSE), "B_start")
  s$B_start <- c(1, 1)
  expect_error(simulate_tuna(s, seed = 1, plots = FALSE), "B_start")
})
