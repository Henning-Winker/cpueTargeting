\name{plot_indices}
\alias{plot_indices}
\title{Plot standardised indices against the true biomass trend}
\usage{
plot_indices(comparison, models = NULL, ci = FALSE,
  type = c("lines", "smooth", "points+smooth"), span = 0.5)
}
\arguments{
\item{comparison}{Output of \code{\link{standardize_tuna_indices}}, or a named
list of such outputs (one per scenario), which are shown as rows.}
\item{models}{Optional subset of model names to show.}
\item{ci}{Logical; show 95\% confidence intervals.}
\item{type}{How to draw the indices: \code{"lines"} (default) joins the
annual values; \code{"smooth"} shows only a loess smoother of each index,
which makes trends easier to compare; \code{"points+smooth"} shows the annual
values as points with the smoother on top.}
\item{span}{Smoothing span of the loess smoother (larger = smoother).}
}
\value{A ggplot object.}
\description{Plot standardised indices against the true biomass trend.}
