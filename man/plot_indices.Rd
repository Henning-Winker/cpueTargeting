\name{plot_indices}
\alias{plot_indices}
\title{Plot standardised indices against the true biomass trend}
\usage{
plot_indices(comparison, models = NULL, ci = FALSE)
}
\arguments{
\item{comparison}{Output of \code{\link{standardize_tuna_indices}}, or a named
list of such outputs (one per scenario), which are shown as rows.}
\item{models}{Optional subset of model names to show.}
\item{ci}{Logical; show 95\% confidence intervals.}
}
\value{A ggplot object.}
\description{Plot standardised indices against the true biomass trend.}
