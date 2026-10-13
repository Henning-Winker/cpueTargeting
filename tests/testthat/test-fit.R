# Two tactics: tactic 1 targets species A, tactic 2 targets B. A declines,
# B is flat, and the share of tactic 1 increases, so nominal CPUE of A is
# biased upwards. Vessels differ in catchability.
make_fit_data <- function(n = 2500, seed = 2) {
  set.seed(seed)
  years <- 2001:2010
  year <- sample(years, n, replace = TRUE)
  iy <- year - 2000
  true_a <- exp(seq(0, -0.8, length.out = 10))
  tactic <- ifelse(stats::runif(n) < seq(0.15, 0.85, length.out = 10)[iy], 1, 2)
  vessels <- sprintf("V%02d", 1:20)
  veff <- stats::setNames(stats::rnorm(20, 0, 0.3), vessels)
  vessel <- sample(vessels, n, replace = TRUE)
  hours <- round(stats::runif(n, 3, 10), 1)
  q_a <- c(1, 0.12)[tactic]
  q_b <- c(0.12, 1)[tactic]
  v <- exp(veff[vessel])
  cpue_a <- mgcv::rTweedie(5 * true_a[iy] * q_a * v, p = 1.4, phi = 1)
  cpue_b <- mgcv::rTweedie(5 * q_b * v, p = 1.4, phi = 1)
  cpue_c <- mgcv::rTweedie(rep(2, n) * v, p = 1.4, phi = 1)
  d <- data.frame(op = seq_len(n), year = year, vessel = vessel, hours = hours,
                  month = sample(1:12, n, replace = TRUE),
                  lon = stats::runif(n, 0, 10) + 3 * (tactic == 1),
                  lat = stats::runif(n, 0, 10), tactic = tactic,
                  A = cpue_a * hours, B = cpue_b * hours, C = cpue_c * hours)
  list(data = d, true = data.frame(year = years, true = true_a))
}

prepare <- function(fd) {
  x <- cpue_data(fd$data, year = "year", species = c("A", "B", "C"), effort = "hours",
                 id = "op", vessel = "vessel", lon = "lon", lat = "lat",
                 covariates = c("month", "tactic"))
  x <- select_species(x, threshold = 1)
  add_targeting(x, n_pc = 1, n_clusters = 2, method = "kmeans")
}

get_gam_for_test <- function(fit, model) fit$models[[model]]$fits$main

rmse_log <- function(fit, model, true) {
  idx <- cpue_index(fit, model)
  t <- true$true / mean(true$true)
  sqrt(mean((log(idx$relative) - log(t))^2))
}

test_that("DPC and Cluster models correct the targeting bias", {
  skip_if_not_installed("mgcv")
  fd <- make_fit_data()
  x <- prepare(fd)
  fit <- fit_cpue(x, "A", models = c("Base", "DPC", "Cluster"), vessel = TRUE)
  expect_s3_class(fit, "cpue_fit")
  expect_equal(levels(fit$indices$model), c("Base", "DPC", "Cluster"))
  base <- rmse_log(fit, "Base", fd$true)
  expect_lt(rmse_log(fit, "DPC", fd$true), 0.6 * base)
  expect_lt(rmse_log(fit, "Cluster", fd$true), 0.6 * base)
  idx <- cpue_index(fit, "DPC")
  expect_equal(nrow(idx), 10)
  expect_true(all(is.finite(idx$se_log) & idx$se_log > 0))
  expect_equal(mean(idx$relative), 1)
  expect_true(all(idx$lower < idx$relative & idx$upper > idx$relative))
  # the reference set averages over both tactics
  expect_setequal(as.character(unique(fit$reference$FT)), c("1", "2"))
})

test_that("diagnostics: table, influence, residuals, effects, vessels", {
  skip_if_not_installed("mgcv")
  fd <- make_fit_data()
  x <- prepare(fd)
  fit <- fit_cpue(x, "A", models = c("Base", "DPC", "Cluster+S"), vessel = TRUE,
                  spatial_k = 10)
  dg <- cpue_diagnostics(fit)
  expect_equal(nrow(dg), 3)
  expect_true(all(dg$tweedie_p > 1 & dg$tweedie_p < 2))
  expect_true(all(dg$vessel_sd > 0.1 & dg$vessel_sd < 0.6))
  expect_equal(min(dg$dAIC), 0)
  inf <- influence_cpue(fit, "DPC")
  expect_true(all(c("s(PC1)", "s(vessel_re)") %in% levels(inf$term)))
  pc_inf <- inf[inf$term == "s(PC1)", ]
  expect_gt(stats::cor(pc_inf$year, pc_inf$influence), 0.8)
  expect_s3_class(plot_influence(fit, "DPC"), "ggplot")
  expect_s3_class(plot_index(fit, true = fd$true), "ggplot")
  expect_s3_class(plot_step(fit), "ggplot")
  res <- plot_residuals(fit, "DPC", by = "FT")
  expect_named(res, c("qq", "fitted", "year", "by"))
  r <- cpueTargeting:::cpue_residuals(get_gam_for_test(fit, "DPC"))
  expect_equal(attr(r, "type"), "quantile")
  expect_lt(abs(mean(r)), 0.1)
  expect_lt(abs(stats::sd(r) - 1), 0.1)
  eff <- plot_effects(fit, "DPC")
  expect_true("numeric" %in% names(eff))
  eff2 <- plot_effects(fit, "Cluster+S")
  expect_true(all(c("factor", "space") %in% names(eff2)))
  expect_s3_class(plot_vessel_effects(fit, "DPC"), "ggplot")
  f <- tempfile(fileext = ".csv")
  tab <- export_index(fit, f, model = "DPC")
  back <- utils::read.csv(f)
  expect_equal(names(back), names(tab))
  expect_equal(nrow(back), 10)
})

test_that("delta models, offset, covariate terms and bam work", {
  skip_if_not_installed("mgcv")
  fd <- make_fit_data(n = 1500)
  x <- prepare(fd)
  for (fam in c("delta_lognormal", "delta_gamma")) {
    fit <- fit_cpue(x, "A", models = c("Base", "DPC"), family = fam)
    idx <- cpue_index(fit, "DPC")
    expect_true(all(is.finite(idx$index) & idx$se_log > 0))
    expect_equal(nrow(cpue_diagnostics(fit)), 4)
    expect_s3_class(plot_influence(fit, "DPC", component = "binomial"), "ggplot")
    expect_named(plot_residuals(fit, "DPC"), c("qq", "fitted", "year"))
  }
  fo <- fit_cpue(x, "A", models = "DPC", offset = TRUE)
  fr <- fit_cpue(x, "A", models = "DPC", offset = FALSE)
  expect_gt(stats::cor(cpue_index(fo)$relative, cpue_index(fr)$relative), 0.9)
  fm <- fit_cpue(x, "A", terms = 's(month, bs = "cc", k = 6)',
                 knots = list(month = c(0.5, 12.5)))
  expect_equal(names(fm$models), c("Base", "Covariates", "DPC", "Cluster"))
  expect_true("month" %in% fm$settings$average_over)
  fb <- fit_cpue(x, "A", models = c("Base", "DPC"), engine = "bam", vessel = TRUE)
  expect_true(all(is.finite(cpue_index(fb, "DPC")$se_log)))
  fa <- fit_cpue(x, "A", models = "Base",
                 add_models = c(Tactic_PC = "fyear + FT + s(PC1, k = 4)"))
  expect_equal(names(fa$models), c("Base", "Tactic_PC"))
})

test_that("informative errors", {
  skip_if_not_installed("mgcv")
  fd <- make_fit_data(n = 600)
  x <- cpue_data(fd$data, year = "year", species = c("A", "B", "C"), id = "op")
  expect_error(fit_cpue(x, "A", models = "DPC"), "add_targeting")
  expect_error(fit_cpue(x, "A", models = "Base+S"), "positions")
  expect_error(fit_cpue(x, "A", vessel = TRUE), "vessel")
  expect_error(fit_cpue(x, "A", add_models = c(M = "s(month)")), "fyear")
})
