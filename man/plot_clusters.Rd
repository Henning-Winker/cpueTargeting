\name{plot_clusters}
\alias{plot_clusters}
\title{What do the catch-composition clusters capture?}
\usage{
plot_clusters(comparison, switch_years = NULL, max_points = 5000)
}
\arguments{
\item{comparison}{Output of \code{\link{standardize_cpue}} (or of
\code{\link{bind_comparisons}}).}
\item{switch_years}{Optional years to mark with vertical lines in the time
plot, e.g. the years in which prices change (\code{c(10.5, 20.5)}).}
\item{max_points}{Maximum number of sets drawn in the PC plot.}
}
\value{A list of ggplots: \code{composition}, \code{time}, \code{space},
\code{depth} (tuna model only) and \code{pca}.}
\description{
Diagnostic plots for the fishing tactics (clusters) used by the
\code{Cluster} models, following the CPUE guidelines: the species composition
of each cluster, the share of sets in each cluster by year, the most common
cluster in each cell, the realised hook depth by cluster (tuna model only)
and the sets in PC space coloured by cluster. Clusters are labelled by their
number and main species, e.g. \code{"2: ALB"}.
}
\details{
\code{hook_depth} is latent truth from the operating model; it is used here
only to show what the clusters pick up, not in any standardisation model.
}
\examples{
\dontrun{
sim <- simulate_tuna(seed = 123)
fit <- standardize_cpue(sim, models = c("Base", "Cluster+S"), engine = "bam")
p <- plot_clusters(fit, switch_years = c(10.5, 20.5))
p$composition
p$time
}
}
