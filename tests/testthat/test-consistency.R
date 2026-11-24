# Simulated logbook with injected recording problems:
# * two species codes replaced by new codes from 2012 (coding break);
# * the last year extracted for the target only, by new vessels.
make_inconsistent <- function() {
  sim <- simulate_logbook(seed = 1, plots = FALSE)
  lb <- sim$logbook
  late <- lb$year >= 2012
  lb$species[late & lb$species == "Conger conger"] <- "Conger sp."
  lb$species[late & lb$species == "Dentex dentex"] <- "Dentex sp."
  last <- max(lb$year)
  lb <- lb[lb$year < last | lb$species == "Pagellus bogaraveo", ]
  lb$vessel[lb$year == last] <- paste0("N", lb$vessel[lb$year == last])
  list(clean = sim$logbook, dirty = lb, last = last)
}
to_x <- function(lb) cpue_data(lb, year = "year", species = "species", catch = "kg",
                               id = "trip", vessel = "vessel")

test_that("species_consistency() finds no break in consistent data", {
  sim <- simulate_logbook(seed = 1, plots = FALSE)
  sc <- species_consistency(to_x(sim$logbook))
  expect_s3_class(sc, "species_consistency")
  expect_equal(nrow(sc$breaks), 0)
  expect_equal(nrow(sc$blocks), 1)
  expect_true(all(sc$occurrence >= 0 & sc$occurrence <= 1))
})

test_that("species_consistency() detects a coding break and suggests replacements", {
  d <- make_inconsistent()
  sc <- species_consistency(to_x(d$dirty))
  expect_true(2012 %in% sc$breaks$year)
  expect_true(all(c("Conger conger", "Dentex dentex") %in% names(sc$suggested)))
  expect_equal(unname(sc$suggested["Conger conger"]), "Conger sp.")
  expect_equal(unname(sc$suggested["Dentex dentex"]), "Dentex sp.")
  expect_true(any(sc$blocks$to == 2011))
  # the target-only last year is a loss of species, with no suggestion to the target
  expect_match(sc$breaks$type[sc$breaks$year == d$last], "species lost")
  expect_false("Pagellus bogaraveo" %in% sc$suggested)
  expect_match(sc$breaks$type[sc$breaks$year == 2012], "coding change")
  expect_s3_class(plot_species_consistency(sc), "ggplot")
  expect_s3_class(plot_species_consistency(sc, "years"), "ggplot")
  expect_s3_class(plot_species_consistency(sc, "pairs"), "ggplot")
  expect_output(print(sc), "Candidate coding breaks")
})

test_that("harmonise_species() merges codes, keeps totals and removes the break", {
  d <- make_inconsistent()
  x <- to_x(d$dirty)
  sc <- species_consistency(x)
  x2 <- harmonise_species(x, sc$suggested)
  expect_false(any(c("Conger.conger", "Dentex.dentex") %in% x2$species))
  expect_equal(sum(x2$data[, x2$species]), sum(x$data[, x$species]))
  expect_equal(x2$data$Conger.sp., x$data$Conger.sp. + x$data$Conger.conger)
  expect_equal(tail(x2$log$step, 1), "harmonise species")
  sc2 <- species_consistency(subset(x2, year < d$last))
  expect_false(2012 %in% sc2$breaks$year)
  # a new target name is created
  x3 <- harmonise_species(x, c("Conger conger" = "Conger group", "Conger sp." = "Conger group"))
  expect_true("Conger.group" %in% x3$species)
  expect_equal(unname(x3$species_labels["Conger.group"]), "Conger group")
  expect_message(harmonise_species(x, c("Not a fish" = "Conger sp.")), "Not found")
})

test_that("vessel_consistency() flags a year with a new vessel set", {
  d <- make_inconsistent()
  vc <- vessel_consistency(to_x(d$dirty))
  expect_s3_class(vc, "vessel_consistency")
  expect_true(d$last %in% vc$flags$year)
  expect_equal(vc$years$ops_new[vc$years$year == d$last], 1)
  expect_equal(vc$years$overlap[vc$years$year == d$last], 0)
  vc0 <- vessel_consistency(to_x(d$clean))
  expect_false(d$last %in% vc0$flags$year)
  expect_equal(sum(vc$cohort[1, ]), 1)
  expect_s3_class(plot_vessel_consistency(vc), "ggplot")
  expect_s3_class(plot_vessel_consistency(vc, "cohort"), "ggplot")
})

test_that("composition_consistency() flags a target-only year", {
  d <- make_inconsistent()
  cc <- composition_consistency(to_x(d$dirty), target = "Pagellus bogaraveo")
  f <- cc$flags[cc$flags$year == d$last, ]
  expect_true(all(c("single-species operations", "operations with target") %in% f$indicator))
  expect_equal(cc$table$`operations with target`[cc$table$year == d$last], 1)
  cc0 <- composition_consistency(to_x(d$clean), target = "Pagellus bogaraveo")
  expect_equal(nrow(cc0$flags), 0)
  expect_s3_class(plot_composition_consistency(cc), "ggplot")
})

test_that("data_consistency() combines the flags of all checks", {
  d <- make_inconsistent()
  dc <- data_consistency(to_x(d$dirty), target = "Pagellus bogaraveo")
  expect_s3_class(dc, "data_consistency")
  expect_true(all(c("species recording", "vessel set", "composition") %in% dc$flags$check))
  expect_output(print(dc), "Suggested species harmonisation")
  dc0 <- data_consistency(to_x(d$clean), target = "Pagellus bogaraveo")
  expect_equal(nrow(dc0$flags), 0)
})
