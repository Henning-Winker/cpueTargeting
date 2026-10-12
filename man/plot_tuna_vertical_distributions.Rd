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
\description{The first panel shows the species depth profiles (for a
thermocline anomaly of zero). The other panels show the expected hook-depth
distribution of each HBF level, combining between-set error and within-set
spread, drawn like Fig. 2 of Bigelow et al. (2002).}
\examples{
plot_tuna_vertical_distributions()
}
