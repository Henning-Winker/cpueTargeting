\name{plot_tuna_vertical_distributions}
\alias{plot_tuna_vertical_distributions}
\title{Plot fish vertical profiles and HBF hook-depth profiles}
\usage{
plot_tuna_vertical_distributions(
  settings = tuna_settings(),
  max_depth = 500,
  file = NULL
)
}
\arguments{
\item{settings}{Settings list, usually returned by \code{tuna_settings()}.}
\item{max_depth}{Maximum depth, in metres, shown in the plot.}
\item{file}{Optional file path. If supplied, saves a 300-dpi PNG and still
returns the plot object.}
}
\value{A \code{ggplot} object.}
\description{Plots species-specific vertical availability and marginal
hook-depth distributions for each HBF level. Hook curves combine between-set
depth uncertainty and within-set hook spread.}
\examples{plot_tuna_vertical_distributions()}
