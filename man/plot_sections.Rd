\name{plot_sections}
\alias{plot_sections}
\title{Species distributions in three dimensions: sections and depth layers}
\usage{
plot_sections(sim, year = 1, layers = c(0, 100, 200, 300, 500))
}
\arguments{
\item{sim}{Output of \code{\link{simulate_tuna}}.}
\item{year}{Year whose thermocline anomaly is used (matters only when the
thermocline moves between years).}
\item{layers}{Depth-layer boundaries (m) for the layer maps.}
}
\value{A list of two ggplots, \code{section} and \code{layers}.}
\description{Combines the horizontal density of each species with its
vertical profile in each cell (including the thermocline anomaly of the
chosen year) and shows the result as a latitude x depth section (density
averaged over x) and as maps of the density in depth layers, species by
layer. Horizontal densities are the initial (unfished) distribution.}
\examples{
\dontrun{
sim <- simulate_tuna(seed = 1)
p <- plot_sections(sim)
p$section
p$layers
}
}
