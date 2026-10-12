\name{plot_tuna_vertical_distributions}
\alias{plot_tuna_vertical_distributions}
\title{Plot species vertical profiles and hook-depth profiles by HBF}
\usage{
plot_tuna_vertical_distributions(settings = tuna_settings(), file = NULL)
}
\arguments{
\item{settings}{Settings list, usually from \code{\link{tuna_settings}}.}
\item{file}{Optional file path; if supplied, a 300-dpi PNG is saved.}
}
\value{A ggplot object.}
\description{Fish profiles are shown for a cell with zero thermocline anomaly;
hook profiles are the expected distributions for each HBF level at the mean
mainline angle without current shoaling.}
\examples{
plot_tuna_vertical_distributions()
}
