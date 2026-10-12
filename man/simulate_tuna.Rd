\name{simulate_tuna}
\alias{simulate_tuna}
\title{Simulate a longline fishery for four tuna species with vertical targeting}
\usage{
simulate_tuna(settings = tuna_settings(), seed = NULL, plots = TRUE)
}
\arguments{
\item{settings}{A named list, usually a modified \code{\link{tuna_settings}}.}
\item{seed}{Optional integer seed for reproducibility.}
\item{plots}{Logical; build the diagnostic ggplots in \code{Plots}.}
}
\value{A list with
\itemize{
\item \code{Data}: one row per set and species with \code{year},
\code{SetID}, \code{X}, \code{Y}, \code{Xobs}, \code{Yobs}, \code{vessel},
\code{HBF}, \code{sst}, \code{thermocline} (anomaly in that year, m), \code{Species},
\code{Sp}, \code{CPUE}, and the latent truth \code{hook_depth} (median hook
depth of the set), \code{vertical_avail} (realised vertical overlap) and
\code{hbs_avail} (overlap predicted by a habitat model from HBF, the expected
hook profile and the average thermocline, as in a habitat-based
standardisation);
\item \code{Biomass}: annual nominal and true indices, scaled to mean 1;
\item \code{Biomass.XY}: initial density by cell and species;
\item \code{Env}: SST and the fixed thermocline anomaly by cell;
\item \code{Thermocline}: thermocline anomaly by cell and year;
\item \code{Plots}: diagnostic ggplots (if \code{plots = TRUE});
\item \code{settings}: the settings used.
}}
\description{Extends the spatial multispecies simulator of Thorson et al.
(2017) and Ono et al. (2015), as implemented in
\code{\link{simulate_bathymetric}}, with a pelagic horizontal habitat (SST),
species vertical distributions, longline hook depth by HBF and vertical
targeting, following the habitat-based approach of Bigelow et al. (2002). A
full description of both models is in
\code{system.file("appendix", "cpueTargeting_simulator_appendix.html",
package = "cpueTargeting")}. The
catchability of a set for a species is the vessel's catchability times the
overlap between the set's hook-depth distribution and the species' vertical
distribution in that cell (Bigelow's effective effort). See
\code{\link{tuna_settings}} for the processes and parameters.}
\references{Bigelow, K.A., Hampton, J. and Miyabe, N. (2002) Application of a
habitat-based model to estimate effective longline fishing effort and relative
abundance of Pacific bigeye tuna (\emph{Thunnus obesus}). Fisheries
Oceanography 11, 143-155.

Thorson, J.T., Fonner, R., Haltuch, M.A., Ono, K. and Winker, H. (2017)
Accounting for spatiotemporal variation and fisher targeting when estimating
abundance from multispecies fishery data. Canadian Journal of Fisheries and
Aquatic Sciences 74, 1794-1807.

Ono, K., Punt, A.E. and Hilborn, R. (2015) Think outside the grids: an
objective approach to define spatial strata for catch and effort analysis.
Fisheries Research 170, 89-101.}
\examples{
\dontrun{
om <- simulate_tuna(seed = 42)
head(om$Data)
om$Plots$CPUE
}
}
