test_that("Tweedie generator has the right mean and variance", {
  set.seed(1)
  x <- cpueTargeting:::rtweedie_cp(rep(0.1, 2e5), power = 1.2, phi = 0.1)
  expect_equal(mean(x), 0.1, tolerance = 0.02)
  expect_equal(var(x), 0.1 * 0.1^1.2, tolerance = 0.05)
  expect_true(any(x == 0))
})

test_that("Gaussian random field has roughly the requested SD", {
  set.seed(2)
  sds <- replicate(20, sd(cpueTargeting:::simulate_grf(40, 40, sd = 2, scale = 4)))
  expect_equal(mean(sds), 2, tolerance = 0.15)
})

test_that("stable distribution is invariant under movement", {
  set.seed(3)
  x <- rep(1:10, 10) + 0.5
  y <- rep(1:10, each = 10) + 0.5
  habitat <- runif(100)
  k <- cpueTargeting:::dispersal_kernels(x, y, 1)[[1]]
  op <- cpueTargeting:::movement_operator(k, habitat)
  moved <- cpueTargeting:::move_biomass(op, op$stable)
  expect_equal(moved, op$stable, tolerance = 1e-10)
  # movement conserves biomass and matches the explicit movement matrix
  b <- runif(100)
  m <- k * rep(habitat, each = 100)
  m <- m / rowSums(m)
  expect_equal(cpueTargeting:::move_biomass(op, b), drop(t(m) %*% b))
})

test_that("catch never exceeds biomass and is shared among sets", {
  biomass <- matrix(c(10, 0, 5, 2), 2, 2)
  cell <- c(1, 1, 2)
  q <- matrix(c(0.1, 0.3, 0.2, 0.5, 0, 0.1), 3, 2)
  out <- cpueTargeting:::allocate_catch(biomass, cell, q)
  expect_true(all(out$cell <= biomass))
  expect_equal(rowsum(out$set, cell), out$cell, ignore_attr = TRUE)
  expect_equal(out$cell[1, 1], 10 * (1 - exp(-0.4)))
})

test_that("Schaefer update applies r by species", {
  b <- matrix(5, 4, 2)
  out <- cpueTargeting:::schaefer_update(b, b * 2, r = c(0.1, 0.5), catch = 0 * b)
  expect_equal(out[, 1], rep(5 + 0.1 * 5 * 0.5, 4))
  expect_equal(out[, 2], rep(5 + 0.5 * 5 * 0.5, 4))
})

test_that("coarse coordinates are block centres", {
  expect_equal(cpueTargeting:::coarse_coord(c(1.5, 5.5, 6.5, 50.5), 5),
               c(3.5, 3.5, 8.5, 48.5))
})

test_that("hook depth increases with HBF and matches Bigelow's catenary", {
  s <- tuna_settings()
  mid <- function(h) max(cpueTargeting:::catenary_depth(seq_len(h), h, 60, s))
  depths <- vapply(c(5, 10, 15, 20), mid, numeric(1))
  expect_true(all(diff(depths) > 0))
  expect_equal(round(depths), c(132, 203, 277, 348))
  edges <- seq(0, s$max_depth, by = s$depth_bin)
  prof <- cpueTargeting:::nominal_hook_profiles(s, edges)
  expect_equal(rowSums(prof), rep(1, 4))
})
