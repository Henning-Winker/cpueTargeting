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
prices, so targeting shifts over time.}
\details{Compared with the original script, the dynamics are identical except
for three corrections: species-specific \code{r} is now applied by species
(the original recycled \code{r} over cells), the initial distribution is the
exact stable distribution of the movement process, and reported block centres
are aligned with the grid.}
\examples{
\dontrun{
om <- simulate_bathymetric(seed = 42)
head(om$Data)
om$Plots$CPUE
}
}
