# Internal building blocks shared by the bathymetric and tuna operating models.
# None of these are exported; they are kept small so each step of the OM can be
# read and checked on its own.

#' Simulate a stationary Gaussian random field on a regular grid
#'
#' Circulant embedding with a Gaussian covariance `exp(-(d / scale)^2)`, i.e.
#' the same covariance model as `RandomFields::RMgauss(var = sd^2, scale)`
#' used in the original Run 22 code. The grid is padded to twice its size to
#' avoid wrap-around correlation; small negative eigenvalues are set to zero.
#'
#' @param nx,ny Number of cells in x and y.
#' @param sd Marginal standard deviation.
#' @param scale Correlation scale in cells.
#' @return An `nx` by `ny` matrix.
#' @noRd
simulate_grf <- function(nx, ny, sd = 1, scale = 4) {
  if (sd == 0) return(matrix(0, nx, ny))
  mx <- 2 * nx
  my <- 2 * ny
  dx <- pmin(0:(mx - 1), mx - 0:(mx - 1))
  dy <- pmin(0:(my - 1), my - 0:(my - 1))
  cov_base <- exp(-outer(dx^2, dy^2, "+") / scale^2)
  lambda <- Re(stats::fft(cov_base))
  lambda[lambda < 0] <- 0
  z <- matrix(complex(real = stats::rnorm(mx * my),
                      imaginary = stats::rnorm(mx * my)), mx, my)
  field <- Re(stats::fft(sqrt(lambda / (mx * my)) * z))
  sd * field[seq_len(nx), seq_len(ny)]
}

#' Random draws from a compound Poisson-gamma (Tweedie, 1 < p < 2) distribution
#'
#' Vectorised over `mu`; mean `mu` and variance `phi * mu^power`. Replaces
#' `tweedie::rtweedie()`, which loops internally and was the main reason for
#' the `tweedie` dependency.
#' @noRd
rtweedie_cp <- function(mu, power = 1.2, phi = 1) {
  if (power <= 1 || power >= 2) stop("Tweedie power must be between 1 and 2.",
                                     call. = FALSE)
  lambda <- mu^(2 - power) / (phi * (2 - power))
  shape <- (2 - power) / (power - 1)
  scale <- phi * (power - 1) * mu^(power - 1)
  n_events <- stats::rpois(length(mu), lambda)
  out <- numeric(length(mu))
  pos <- n_events > 0
  out[pos] <- stats::rgamma(sum(pos), shape = n_events[pos] * shape,
                            scale = scale[pos])
  out
}

#' Habitat preference functions of the original Run 22 code
#'
#' `Norm` uses `par2` as a coefficient of variation (sd = par1 * par2) and
#' `Lognorm` uses `par1` as the mean and `par2` as the CV, as in the original.
#' @noRd
habitat_density <- function(z, func, par1, par2) {
  switch(func,
    Norm = stats::dnorm(z, par1, par2 * par1),
    Exp = stats::dexp(z, rate = 1 / par1),
    Lognorm = {
      sdlog <- sqrt(log(par2^2 + 1))
      stats::dlnorm(z, log(par1) - sdlog^2 / 2, sdlog)
    },
    Unif = stats::dunif(z, par1, par2),
    stop("Unknown habitat function: ", func, call. = FALSE)
  )
}

#' Movement operator for one species
#'
#' Fish in cell i move to cell j with probability
#' `M[i, j] = K[i, j] * h[j] / c[i]`, where `K` is an exponential distance
#' kernel, `h` is habitat suitability and `c = K %*% h` normalises rows. This
#' is the movement matrix of the original code, but it is never formed: because
#' `K` is symmetric, redistribution `t(M) %*% b` equals `h * K %*% (b / c)`,
#' and the stable distribution has the closed form `h * c` (detailed balance),
#' replacing the original 1000-step iteration.
#'
#' @param kernel Symmetric dispersal kernel (cells x cells).
#' @param habitat Habitat suitability by cell.
#' @noRd
movement_operator <- function(kernel, habitat) {
  norm <- drop(kernel %*% habitat)
  stable <- habitat * norm
  list(kernel = kernel, habitat = habitat, norm = norm,
       stable = stable / sum(stable))
}

#' Redistribute biomass with a movement operator
#' @noRd
move_biomass <- function(op, biomass) {
  ratio <- ifelse(op$norm > 0, biomass / op$norm, 0)
  op$habitat * drop(op$kernel %*% ratio)
}

#' Exponential dispersal kernels, computed once per distinct dispersal distance
#' @noRd
dispersal_kernels <- function(x, y, dispersal) {
  dist <- as.matrix(stats::dist(cbind(x, y)))
  distinct <- unique(dispersal)
  kernels <- lapply(distinct, function(d) stats::dexp(dist, rate = 1 / d))
  kernels[match(dispersal, distinct)]
}

#' Annual number of sets: lognormal error around a logistic increase
#' @noRd
annual_effort <- function(year, total, cv) {
  sig <- sqrt(log(cv^2 + 1))
  mean_effort <- total / (1 + exp(-0.1 * year))
  trunc(stats::rlnorm(1, log(mean_effort) - sig^2 / 2, sig))
}

#' Share the catch in each cell among the sets fishing there
#'
#' Within a cell, total catch of a species is `B * (1 - exp(-sum(q)))` and is
#' split among sets in proportion to their effective catchability `q`, as in
#' the original code.
#'
#' @param biomass Cells x species biomass matrix.
#' @param cell Cell index of each set.
#' @param q Sets x species effective catchability.
#' @return A list with `set` (sets x species catch) and `cell` (cells x species).
#' @noRd
allocate_catch <- function(biomass, cell, q) {
  q_cell <- matrix(0, nrow(biomass), ncol(biomass))
  summed <- rowsum(q, cell)
  q_cell[as.integer(rownames(summed)), ] <- summed
  catch_cell <- biomass * (1 - exp(-q_cell))
  q_total <- q_cell[cell, , drop = FALSE]
  share <- ifelse(q_total > 0, q / q_total, 0)
  list(set = catch_cell[cell, , drop = FALSE] * share, cell = catch_cell)
}

#' Schaefer surplus production in each cell (vectorised over species)
#'
#' `r` has one value per species and is applied by column. (The original code
#' used `r * B`, which recycles `r` over cells rather than species.)
#' @noRd
schaefer_update <- function(biomass, carrying, r, catch) {
  r_mat <- matrix(r, nrow(biomass), ncol(biomass), byrow = TRUE)
  growth <- ifelse(carrying > 0, r_mat * biomass * (1 - biomass / carrying), 0)
  out <- biomass + growth - catch
  out[!is.finite(out) | out < 0] <- 0
  out
}

#' Mid-points of coarse reporting blocks (e.g. 5 x 5 cells)
#'
#' Cell centres are at `index + 0.5`; a block of `res` cells is reported at its
#' centre, so cells 1-5 map to 3.5, cells 6-10 to 8.5, etc.
#' @noRd
coarse_coord <- function(coord, res) {
  index <- coord - 0.5
  (ceiling(index / res) - 1) * res + (res + 1) / 2 + 0.5
}

#' Annual nominal (mean CPUE) and true biomass indices, both scaled to mean 1
#' @noRd
index_table <- function(data, biomass, species) {
  years <- seq_len(nrow(biomass))
  nominal <- tapply(data$CPUE, list(data$year, data$Sp), mean)
  nominal <- nominal[as.character(years), species, drop = FALSE]
  out <- lapply(seq_along(species), function(s) {
    data.frame(year = years, Sp = species[s],
               Nominal = nominal[, s] / mean(nominal[, s], na.rm = TRUE),
               True = biomass[, s] / mean(biomass[, s]))
  })
  out <- do.call(rbind, out)
  out$Sp <- factor(out$Sp, levels = species)
  rownames(out) <- NULL
  out
}

#' Reshape set-by-species records (long) to one row per set
#'
#' @param data Long data with one row per set and species, e.g. `sim$Data`.
#' @param value Column to spread across species (default `"CPUE"`).
#' @return A data frame with set-level columns and one column per species.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_tuna(seed = 1)
#' sets <- as_set_data(sim$Data)
#' head(sets)
#' }
as_set_data <- function(data, value = "CPUE") {
  species <- unique(as.character(data$Sp))
  set_cols <- setdiff(names(data),
    c("Species", "Sp", "CPUE", "vertical_avail", "hbs_avail"))
  sets <- data[!duplicated(data$SetID), set_cols, drop = FALSE]
  sets <- sets[order(sets$SetID), , drop = FALSE]
  key <- match(data$SetID, sets$SetID)
  for (s in species) {
    rows <- data$Sp == s
    col <- numeric(nrow(sets))
    col[key[rows]] <- data[[value]][rows]
    sets[[s]] <- col
  }
  rownames(sets) <- NULL
  sets
}
