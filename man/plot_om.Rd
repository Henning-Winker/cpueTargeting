\name{plot_om}
\alias{plot_om}
\title{Diagnostic plots for an operating-model run}
\usage{
plot_om(sim)
}
\arguments{
\item{sim}{Output of \code{\link{simulate_bathymetric}} or \code{\link{simulate_tuna}}.}
}
\value{A named list of ggplots: \code{Env} (bathymetry, or SST for the tuna
model), \code{Thermocline} (fixed thermocline anomaly; tuna model only), \code{SpDens} (initial densities),
\code{Hist} (CPUE frequencies), \code{CPUE} (nominal versus true index) and,
for the tuna model, \code{HBF} (annual use of each HBF level) and
\code{Avail} (see \code{\link{plot_availability}}). \code{Depth} is
kept as an alias of \code{Env} for the bathymetric model.}
\description{Diagnostic plots for an operating-model run.}
