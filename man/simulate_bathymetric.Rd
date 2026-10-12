\name{simulate_bathymetric}
\alias{simulate_bathymetric}
\title{Simulate the four-species bathymetric Run 22 scenario}
\usage{
simulate_bathymetric(settings = bathymetric_settings(), seed = NULL, plots = TRUE)
}
\arguments{
\item{settings}{A named settings list, usually from \code{\link{bathymetric_settings}}.}
\item{seed}{Optional integer seed for reproducibility.}
\item{plots}{Logical; build the diagnostic ggplots in \code{Plots}.}
}
\value{A list with
\itemize{
\item \code{Data}: one row per set and species (\code{year}, \code{SetID},
\code{X}, \code{Y}, \code{Xobs}, \code{Yobs}, \code{Species}, \code{Sp},
\code{CPUE}, \code{depth}, \code{vessel});
\item \code{Biomass}: annual nominal and true indices, scaled to mean 1;
\item \code{Biomass.XY}: initial density by cell and species;
\item \code{Env}: the bathymetry;
\item \code{Plots}: diagnostic ggplots (if \code{plots = TRUE});
\item \code{settings}: the settings used.
}}
\description{Four species are distributed over a simulated bathymetry
according to depth preferences and short-distance movement. Fishing effort is
allocated in proportion to expected catch value, which changes with species
prices, so targeting shifts over time.

This is a re-implementation of the spatial multispecies simulator of Thorson
et al. (2017, Appendix A; function \code{Generate_scenario_data} in the
\code{spatial_DFA} package, \url{https://github.com/James-Thorson/spatial_DFA}),
which generalises the simulator of Ono et al. (2015), with the "Run 22"
settings used for the WKBMSYSPICT CPUE guidelines. A full description is in
\code{system.file("appendix", "cpueTargeting_simulator_appendix.html",
package = "cpueTargeting")}.}
\details{Compared with the original script, the dynamics are identical except
for three corrections: species-specific \code{r} is now applied by species
(the original recycled \code{r} over cells), the initial distribution is the
exact stable distribution of the movement process, and reported block centres
are aligned with the grid.}
\references{Thorson, J.T., Fonner, R., Haltuch, M.A., Ono, K. and Winker, H.
(2017) Accounting for spatiotemporal variation and fisher targeting when
estimating abundance from multispecies fishery data. Canadian Journal of
Fisheries and Aquatic Sciences 74, 1794-1807.

Ono, K., Punt, A.E. and Hilborn, R. (2015) Think outside the grids: an
objective approach to define spatial strata for catch and effort analysis.
Fisheries Research 170, 89-101.}
\examples{
\dontrun{
om <- simulate_bathymetric(seed = 42)
head(om$Data)
om$Plots$CPUE
}
}
