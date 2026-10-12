\name{plot_thermocline_effect}
\alias{plot_thermocline_effect}
\title{Show how a thermocline shift moves species relative to the hooks}
\usage{
plot_thermocline_effect(settings = tuna_settings(), shifts = c(-40, 0, 40))
}
\arguments{
\item{settings}{Settings list, usually from \code{\link{tuna_settings}}.}
\item{shifts}{Thermocline anomalies to show (m, positive = deeper).}
}
\value{A ggplot object.}
\description{Species depth profiles for a shallow, average and deep
thermocline, with the expected hook-depth profiles of each HBF level in grey.
When the thermocline is shallow, deep-set hooks fish below the deep species;
when it is deep, shallow sets miss them.}
\examples{
plot_thermocline_effect()
}
