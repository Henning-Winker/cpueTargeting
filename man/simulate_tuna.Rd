\name{simulate_tuna}
\alias{simulate_tuna}
\title{Simulate spatially distributed tuna with longline hook-depth overlap}
\usage{simulate_tuna(settings = tuna_settings(), seed = NULL)}
\arguments{
\item{settings}{A named settings list, usually from \code{tuna_settings()}.}
\item{seed}{Optional integer seed for reproducibility.}
}
\value{A list with \code{Data} (set-by-species observations, including HBF
and latent hook-depth/overlap truth), \code{Biomass}, spatial biomass
information, and diagnostic \code{Plots}.}
\description{Generates the spatial operating model using species vertical
profiles and a noisy HBF-to-hook-depth relationship. Bathymetry governs
horizontal habitat association; separate vertical profiles govern overlap
between fish and longline hooks.}
\examples{
\dontrun{
om <- simulate_tuna(seed = 42)
head(om$Data)
}
}
