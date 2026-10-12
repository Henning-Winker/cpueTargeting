\name{plot_availability}
\alias{plot_availability}
\title{Vertical availability by HBF over time}
\usage{
plot_availability(sim)
}
\arguments{
\item{sim}{Output of \code{\link{simulate_tuna}}.}
}
\value{A ggplot object.}
\description{Mean realised vertical availability (overlap between hooks and
fish) of the sets fishing with each HBF level, by year and species. With a
fixed thermocline the lines are flat apart from changes in where the fleet
fishes; with a moving thermocline the catchability of a given HBF changes from
year to year, which an HBF factor cannot absorb.}
