\name{plot_thermocline}
\alias{plot_thermocline}
\title{Plot the thermocline anomaly}
\usage{
plot_thermocline(sim, years = NULL)
}
\arguments{
\item{sim}{Output of \code{\link{simulate_tuna}}.}
\item{years}{Years to map (default: six years spread over the series).}
}
\value{A list of two ggplots: \code{series} (mean anomaly by year with the
10-90\% range across cells) and \code{maps} (anomaly fields in selected years).}
\description{Plot the thermocline anomaly.}
