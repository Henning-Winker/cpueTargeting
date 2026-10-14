small_logbook <- function() {
  s <- logbook_settings()
  s$n_years <- 9
  s$trips <- rep(400, 9)
  s$price_periods$years <- c(1, 4, 7)
  s$closure$year <- 7
  s
}

test_that("simulate_logbook() returns a messy logbook and the truth", {
  sim <- simulate_logbook(small_logbook(), seed = 1)
  lb <- sim$logbook
  expect_named(lb, c("trip", "year", "month", "vessel", "hours", "lon", "lat", "species", "kg"))
  expect_true(all(lb$kg > 0))
  expect_true(all(abs(lb$kg / 0.5 - round(lb$kg / 0.5)) < 1e-8))
  expect_true("Blackspot seabream" %in% lb$species)
  expect_true(any(is.na(lb$lon)))
  expect_lt(mean(is.na(lb$lon)), 0.1)
  expect_equal(range(lb$year), c(2000, 2008))
  tr <- sim$trips
  expect_true(all(tr$tactic %in% c("Reef", "Sediment", "Pelagic")))
  expect_equal(mean(tr$n_sets == 2), 0.3, tolerance = 0.05)
  expect_equal(nrow(sim$truth), 9 * 10)
  expect_equal(as.numeric(tapply(sim$truth$relative, sim$truth$species, mean)), rep(1, 10))
  expect_equal(mean(sim$env$substrate == "reef"), 0.25, tolerance = 0.02)
  expect_named(logbook_truth(sim, "Pagellus bogaraveo"), c("year", "true"))
})

test_that("habitats, vessel types, closure and prices act as specified", {
  sim <- simulate_logbook(small_logbook(), seed = 2)
  tr <- sim$trips
  # reef species are denser on reef, sediment species on sediment
  dens <- sim$density
  reef_cell <- sim$env$substrate[match(paste(dens$x, dens$y), paste(sim$env$x, sim$env$y))] == "reef"
  ratio <- function(sp) {
    d <- dens$density[dens$species == sp]
    r <- reef_cell[dens$species == sp]
    mean(d[r]) / mean(d[!r])
  }
  expect_gt(ratio("Dentex dentex"), 3)
  expect_lt(ratio("Merluccius merluccius"), 1 / 3)
  # specialists
  expect_true(all(tr$tactic[tr$vessel_type == "pelagic"] == "Pelagic"))
  expect_false(any(tr$tactic[tr$vessel_type == "demersal"] == "Pelagic"))
  expect_false(any(tr$tactic2[tr$vessel_type == "pelagic"] %in% c("Reef", "Sediment")))
  # closure: no sediment sets deeper than 300 m from year 7
  closed <- tr$year >= 2006 & tr$depth >= 300 & tr$tactic == "Sediment"
  expect_equal(sum(closed), 0)
  # the pelagic boom (years 4-6) raises the share of pelagic trips
  share <- tapply(tr$tactic == "Pelagic", tr$year >= 2003 & tr$year <= 2005, mean)
  expect_gt(share[["TRUE"]], share[["FALSE"]] + 0.05)
  pl <- sim$Plots
  expect_true(all(vapply(pl, inherits, logical(1), "ggplot")))
})

test_that("custom settings are validated", {
  s <- logbook_settings()
  s$n_years <- 3
  expect_error(simulate_logbook(s), "n_years")
  s <- small_logbook()
  s$vessel_types <- c(demersal = 5, other = 5)
  expect_error(simulate_logbook(s), "home")
  s <- small_logbook()
  s$loyalty[["pelagic"]] <- 1.5
  expect_error(simulate_logbook(s), "loyalty")
})

test_that("specialists stay loyal to their home tactics most of the time", {
  sim <- simulate_logbook(small_logbook(), seed = 2, plots = FALSE)
  tr <- sim$trips
  pel <- tr[tr$vessel_type == "pelagic", ]
  expect_gt(mean(pel$loyal), 0.7)
  expect_lt(mean(pel$loyal), 0.9)
  expect_true(any(pel$tactic != "Pelagic"))
  expect_true(all(pel$tactic[pel$loyal] == "Pelagic"))
  dem <- tr[tr$vessel_type == "demersal" & tr$loyal, ]
  expect_true(all(dem$tactic %in% c("Reef", "Sediment")))
  expect_false(any(tr$loyal[tr$vessel_type == "mixed"]))
  # earlier `allowed` settings: specialists never leave their tactics
  s <- small_logbook()
  s$home <- NULL
  s$loyalty <- NULL
  s$allowed <- list(demersal = c("Reef", "Sediment"), pelagic = "Pelagic",
                    mixed = c("Reef", "Sediment", "Pelagic"))
  tr1 <- simulate_logbook(s, seed = 2, plots = FALSE)$trips
  expect_true(all(tr1$tactic[tr1$vessel_type == "pelagic"] == "Pelagic"))
})

test_that("the toolbox recovers the tactics from the logbook", {
  skip_if_not_installed("cluster")
  sim <- simulate_logbook(small_logbook(), seed = 3, plots = FALSE)
  x <- cpue_data(sim$logbook, year = "year", species = "species", catch = "kg",
                 effort = "hours", id = "trip", vessel = "vessel", lon = "lon", lat = "lat",
                 covariates = "month",
                 synonyms = c("Blackspot seabream" = "Pagellus bogaraveo"))
  expect_equal(length(x$species), 10)
  x <- select_species(x, threshold = 0.95)
  x <- add_targeting(x, n_clusters = 3)
  true_tac <- sim$trips$tactic[match(x$data$trip, sim$trips$trip)]
  # there is a clearly pelagic and a clearly demersal cluster; coastal
  # clusters may mix (trolling for coastal pelagics vs shallow handlining)
  pelagic <- as.numeric(true_tac == "Pelagic")
  pcs <- grep("^PC[0-9]+$", names(x$data), value = TRUE)
  r <- vapply(pcs, function(p) abs(stats::cor(x$data[[p]], pelagic)), numeric(1))
  expect_gt(max(r), 0.5)
  share_pelagic <- tapply(pelagic, x$data$FT, mean)
  expect_gt(max(share_pelagic) - min(share_pelagic), 0.4)
})
