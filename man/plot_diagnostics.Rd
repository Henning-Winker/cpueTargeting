\name{plot_diagnostics}
\alias{plot_diagnostics}
\title{Compare the accuracy of standardised indices}
\usage{
plot_diagnostics(comparison, models = NULL)
}
\arguments{
\item{comparison}{Output of \code{\link{standardize_tuna_indices}}, or a named
list of such outputs (one per scenario), shown as differently coloured bars.}
\item{models}{Optional subset of model names to show.}
}
\value{A ggplot object.}
\description{Bar chart of the root mean squared error of the log index
(relative to the true biomass index) by model and species.}
