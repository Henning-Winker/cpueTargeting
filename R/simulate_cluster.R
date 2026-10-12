#' Simulate the 2014 multi-habitat cluster scenario
#'
#' This callable data generator follows the multi-habitat, multi-species
#' targeting example in Winker et al. (2014). Ten species follow exponential
#' biomass trajectories; each trip samples two habitats; effort is split
#' between them; and species catches follow habitat-specific catchability.
#' The original script is retained in `inst/extdata` for reference.
#'
#' @param years Number of annual time steps.
#' @param trips_per_year Number of trips sampled in each year.
#' @param scenario Either `"targeting"` for changing habitat allocation or
#'   `"constant"` for approximately constant allocation.
#' @param habitat_switch Probability that the second habitat uses a separate
#'   catch-rate draw. The default `0.25` matches the source script.
#' @param tweedie_power Tweedie power parameter, between 1 and 2.
#' @param tweedie_dispersion Positive Tweedie dispersion parameter.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with `data` (one row per trip, including habitat, effort
#'   fractions, and species catches), `biomass` (the true annual biomass
#'   trajectories), `effort` (annual habitat allocation probabilities), and
#'   `settings` (the arguments used).
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- simulate_cluster(years = 20, trips_per_year = 100, seed = 1)
#' head(sim$data)
#' }
simulate_cluster <- function(years = 20, trips_per_year = 500,
                           scenario = c("targeting", "constant"),
                           habitat_switch = 0.25,
                           tweedie_power = 1.3,
                           tweedie_dispersion = 10,
                           seed = NULL) {
  scenario <- match.arg(scenario)
  if (length(years) != 1L || years < 1 || years != as.integer(years)) {
    stop("`years` must be a positive integer.", call. = FALSE)
  }
  if (length(trips_per_year) != 1L || trips_per_year < 1 ||
      trips_per_year != as.integer(trips_per_year)) {
    stop("`trips_per_year` must be a positive integer.", call. = FALSE)
  }
  if (length(habitat_switch) != 1L || habitat_switch < 0 ||
      habitat_switch > 1) {
    stop("`habitat_switch` must be between 0 and 1.", call. = FALSE)
  }
  if (length(tweedie_power) != 1L || tweedie_power <= 1 ||
      tweedie_power >= 2 || tweedie_dispersion <= 0) {
    stop("Use 1 < `tweedie_power` < 2 and positive `tweedie_dispersion`.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  species <- c("KOB", "GLBK", "HAKE", "PANG", "CRPN", "SNTR",
               "ROMN", "DRGD", "RSTM", "SHRK")
  n_species <- length(species)
  year <- seq_len(as.integer(years))
  n <- years * trips_per_year

  growth <- stats::runif(n_species) * 0.1 *
    (-1)^round(stats::runif(n_species) * 4)
  initial <- c(200 * exp(stats::rnorm(3, 0, 0.5)),
               50 * exp(stats::rnorm(1, 0, 0.5)),
               200 * exp(stats::rnorm(1, 0, 0.5)),
               50 * exp(stats::rnorm(5, 0, 0.5)))
  biomass <- outer(year - 1, growth, function(t, r) exp(t * r))
  biomass <- sweep(biomass, 2, initial, `*`)
  colnames(biomass) <- species

  effort <- matrix(NA_real_, nrow = 4, ncol = years,
                   dimnames = list(paste0("habitat", 1:4), year))
  for (y in year) {
    error <- stats::rnorm(3, 0, 0.2)
    if (scenario == "targeting") {
      first <- if (y <= ceiling(years / 2)) 0.45 else 0.11
      second <- if (y <= ceiling(years / 2)) 0.10 else 0.45
      third <- 0.3 - 0.008 * (y - 1)
    } else {
      first <- 0.2
      second <- 0.25
      third <- 0.3
    }
    shares <- c(first * exp(error[1]), second * exp(error[2]),
                third * exp(error[3]), 0)
    shares[4] <- max(0.001, 1 - sum(shares[1:3]))
    effort[, y] <- shares / sum(shares)
  }

  # The original workflow randomly re-labelled habitat rows once per run.
  effort <- effort[sample.int(4), , drop = FALSE]
  trip_year <- rep(year, each = trips_per_year)
  idx <- seq_len(n)
  probabilities <- effort[, trip_year, drop = FALSE]
  draw_habitat <- function(prob) {
    u <- stats::runif(n)
    cut1 <- prob[1, ]
    cut2 <- cut1 + prob[2, ]
    cut3 <- cut2 + prob[3, ]
    ifelse(u < cut1, 1L, ifelse(u < cut2, 2L,
                                ifelse(u < cut3, 3L, 4L)))
  }
  habitat_a <- draw_habitat(probabilities)
  habitat_b <- draw_habitat(probabilities)
  effort_a <- stats::runif(n)
  effort_b <- 1 - effort_a

  q <- rbind(
    c(1.00, 0.90, 0.05, 0.10, 0.01, 0.01, 0.01, 0.01, 0.01, 0.50),
    c(0.02, 0.10, 1.00, 0.80, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05),
    c(0.01, 0.01, 0.01, 0.01, 1.00, 0.05, 0.20, 0.20, 0.30, 0.05),
    c(0.01, 0.01, 0.01, 0.01, 0.10, 1.00, 0.80, 0.70, 0.90, 0.30)
  )
  catch_a <- catch_b <- matrix(0, nrow = n, ncol = n_species)
  for (h in 1:4) {
    rows_a <- which(habitat_a == h)
    rows_b <- which(habitat_b == h)
    if (length(rows_a)) for (s in seq_len(n_species)) {
      mu <- q[h, s] * biomass[trip_year[rows_a], s]
      catch_a[rows_a, s] <- tweedie::rtweedie(
        length(rows_a), xi = tweedie_power, mu = mu,
        phi = tweedie_dispersion
      )
    }
    if (length(rows_b)) for (s in seq_len(n_species)) {
      mu <- q[h, s] * biomass[trip_year[rows_b], s]
      catch_b[rows_b, s] <- tweedie::rtweedie(
        length(rows_b), xi = tweedie_power, mu = mu,
        phi = tweedie_dispersion
      )
    }
  }
  shared <- stats::runif(n) >= habitat_switch
  catch_b[shared, ] <- catch_a[shared, , drop = FALSE]
  catches <- catch_a * effort_a + catch_b * effort_b

  data <- data.frame(
    trip = idx, year = trip_year, habitatA = habitat_a,
    habitatB = habitat_b, propA = effort_a, propB = effort_b,
    catches, check.names = FALSE
  )
  names(data)[seq.int(7L, 6L + n_species)] <- species
  rownames(biomass) <- year
  list(data = data, biomass = data.frame(year, biomass, check.names = FALSE),
       effort = data.frame(year, t(effort), check.names = FALSE),
       settings = list(years = years, trips_per_year = trips_per_year,
                       scenario = scenario, habitat_switch = habitat_switch,
                       tweedie_power = tweedie_power,
                       tweedie_dispersion = tweedie_dispersion, seed = seed))
}
