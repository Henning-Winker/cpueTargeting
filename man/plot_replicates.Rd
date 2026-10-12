\name{plot_replicates}
\alias{plot_replicates}
\title{Accuracy of standardised indices across replicates}
\usage{
plot_replicates(replicates, models = NULL)
}
\arguments{
\item{replicates}{Output of \code{\link{run_replicates}}.}
\item{models}{Optional subset of model names to show.}
}
\value{A ggplot object.}
\description{Boxplots of the RMSE of the log index over seeds, by model,
scenario and species.}
