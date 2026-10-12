\name{simulate_bathymetric}
\alias{simulate_bathymetric}
\title{Simulate the four-species bathymetric Run 22 scenario}
\usage{simulate_bathymetric(settings = bathymetric_settings(), seed = NULL)}
\arguments{
\item{settings}{A named settings list, usually from
\code{bathymetric_settings()}.}
\item{seed}{Optional integer seed for reproducibility.}
}
\value{A list with \code{Data}, \code{Biomass}, \code{Biomass.XY}, and
\code{Plots} (ggplot diagnostics).}
\description{Runs the four-species Run 22 spatial operating model. Species
are distributed across a bathymetric landscape, and fishing effort responds
to expected multispecies catch. This model represents horizontal bathymetric
habitat only and does not include the vertical tuna availability process.}
\examples{
\dontrun{
om <- simulate_bathymetric(seed = 42)
head(om$Data)
om$Plots$Depth
}
}
