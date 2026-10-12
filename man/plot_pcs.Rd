\name{plot_pcs}
\alias{plot_pcs}
\title{What do the principal components of the catch composition capture?}
\usage{
plot_pcs(comparison, pcs = c("PC1", "PC2"), depth_bin = 20, min_sets = 30)
}
\arguments{
\item{comparison}{Output of \code{\link{standardize_tuna_indices}} (or of
\code{\link{bind_comparisons}}).}
\item{pcs}{Principal components to show.}
\item{depth_bin}{Width of the hook-depth bins (m).}
\item{min_sets}{Minimum number of sets for a depth bin to be shown.}
}
\value{A list of ggplots: \code{loadings}, \code{space}, \code{depth} (tuna
model only) and \code{time}.}
\description{Plots the PC scores used by the DPC models in four ways: the
species loadings (what each PC means), mean scores by cell (where each tactic
is used), mean scores against the realised hook depth of the set for each HBF
(how the PCs track vertical targeting; tuna model only), and mean scores by
year (how targeting changes over time).}
\details{\code{hook_depth} is latent truth from the operating model; it is
used here only to show what the PCs pick up, not in any standardisation
model.}
