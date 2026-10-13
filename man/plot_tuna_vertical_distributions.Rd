\name{plot_tuna_vertical_distributions}
\alias{plot_tuna_vertical_distributions}
\title{Plot species vertical profiles and hook-depth profiles by HBF}
\usage{
plot_tuna_vertical_distributions(settings = tuna_settings(), overlay = TRUE,
  file = NULL)
}
\arguments{
\item{settings}{Settings list, usually from \code{\link{tuna_settings}}.}
\item{overlay}{Logical; \code{TRUE} (default) overlays the species on each HBF
panel; \code{FALSE} shows the species in a separate first panel.}
\item{file}{Optional file path; if supplied, a 300-dpi PNG is saved.}
}
\value{A ggplot object.}
\description{One panel per HBF level. Grey bars show where the hooks of a set
with that HBF fish (the expected hook-depth distribution, drawn like Fig. 2 of
Bigelow et al. 2002); the coloured translucent areas show where each species
lives (thermocline anomaly of zero). Where a species' area overlaps the grey
bars, that set can catch it. The number next to each species is its vertical
availability to that HBF, i.e. the overlap that multiplies catchability in the
operating model (see \code{\link{plot_availability}}).}
\examples{
plot_tuna_vertical_distributions()
}
