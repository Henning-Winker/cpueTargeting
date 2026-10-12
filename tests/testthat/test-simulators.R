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
  expect_true(all(c("Env", "SpDens", "Hist", "CPUE", "HBF") %in% names(sim$Plots)))
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
  expect_s3_class(plot_diagnostics(ab), "ggplot")
})
