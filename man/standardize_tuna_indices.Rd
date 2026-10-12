\name{standardize_tuna_indices}
\alias{standardize_tuna_indices}
\title{Compare tuna CPUE indices with and without HBF adjustment}
\usage{
standardize_tuna_indices(
  simulation,
  species_names = c("BET", "YFT", "ALB", "SWO"),
  seed = 22062017
)
}
\arguments{
\item{simulation}{Output from \code{simulate_tuna()}.}
\item{species_names}{Labels for the four species, in simulator order.}
\item{seed}{Optional seed for the fishing-tactic k-means classification.}
}
\value{A list with \code{indices} (annual relative indices and the true
biomass index) and \code{diagnostics} (mean squared log error by species and
model).}
\description{Fits base, spatial, direct-principal-component, cluster,
spatial-plus-HBF, and spatial-plus-DPC-plus-HBF models to a tuna simulation.
Latent actual hook depth and vertical overlap are not used as predictors.}
\examples{
\dontrun{
om <- simulate_tuna(seed = 42)
comparison <- standardize_tuna_indices(om)
}
}
