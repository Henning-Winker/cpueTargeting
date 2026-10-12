\name{plot_cluster_2014}
\alias{plot_cluster_2014}
\title{Plot the 2014 cluster simulation}
\usage{plot_cluster_2014(simulation)}
\arguments{\item{simulation}{Output from \code{simulate_cluster()}.}}
\value{A named list of three ggplot objects: \code{biomass},
\code{effort}, and \code{catches}.}
\description{Creates diagnostic ggplot figures for annual biomass,
annual habitat effort allocation, and trip-level catches by species.}
\examples{
\dontrun{
sim <- simulate_cluster(seed = 1)
plots <- plot_cluster_2014(sim)
plots$biomass
}
}
