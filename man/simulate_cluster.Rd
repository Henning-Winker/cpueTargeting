\name{simulate_cluster}
\alias{simulate_cluster}
\title{Simulate the 2014 multi-habitat cluster scenario}
\usage{
simulate_cluster(years = 20, trips_per_year = 500,
  scenario = c("targeting", "constant"), habitat_switch = 0.25,
  tweedie_power = 1.3, tweedie_dispersion = 10, seed = NULL)
}
\arguments{
\item{years}{Number of annual time steps.}
\item{trips_per_year}{Number of trips sampled in each year.}
\item{scenario}{Either \code{"targeting"} for changing habitat allocation or
\code{"constant"} for approximately constant allocation.}
\item{habitat_switch}{Probability that the second habitat uses a separate
catch-rate draw.}
\item{tweedie_power}{Tweedie power parameter, between 1 and 2.}
\item{tweedie_dispersion}{Positive Tweedie dispersion parameter.}
\item{seed}{Optional integer seed for reproducibility.}
}
\value{A list with \code{data}, \code{biomass}, \code{effort}, and
\code{settings}.}
\description{A callable version of the multi-habitat, multi-species cluster
targeting example from Winker et al. (2014), based on the original
\code{SimulationCPUE_3d.R} script.}
\examples{
\dontrun{
sim <- simulate_cluster(years = 20, trips_per_year = 100, seed = 1)
head(sim$data)
}
}
