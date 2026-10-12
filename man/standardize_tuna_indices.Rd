\name{standardize_tuna_indices}
\alias{standardize_tuna_indices}
\title{Standardise simulated CPUE with alternative targeting covariates}
\usage{
standardize_tuna_indices(simulation, models = NULL, vessel_re = FALSE,
  n_clusters = 3, seed = 1)
}
\arguments{
\item{simulation}{Output of \code{\link{simulate_tuna}} or \code{\link{simulate_bathymetric}}.}
\item{models}{Names of the models to fit (default: all that apply).}
\item{vessel_re}{Logical; add a random vessel effect \code{s(vessel, bs = "re")}.}
\item{n_clusters}{Number of fishing tactics for the \code{Cluster} model.}
\item{seed}{Seed for the k-means clustering.}
}
\value{A list with \code{indices} (species, year, model, index with 95\% CI
and the true index), \code{diagnostics} (root mean squared error of the log
index by species and model), and \code{data} (the set-level data used, with
PCs and tactics).}
\description{Fits Tweedie GAMs (\code{mgcv::gam} with \code{tw()}) to each
species and compares the annual indices with the true biomass trend.}
\details{Candidate models (all include \code{year}):
\itemize{
\item \code{Base}: year only.
\item \code{Spatial}: + \code{s(Xobs, Yobs)} on the coarse reporting grid.
\item \code{DPC}: Spatial + \code{s(PC1)} + \code{s(PC2)} from the catch
composition (direct principal components).
\item \code{Cluster}: Spatial + fishing tactic from k-means clustering of the
catch composition.
\item \code{HBF}: Spatial + HBF as a factor (tuna model only).
\item \code{DPC_HBF}: DPC + HBF (tuna model only).
\item \code{HBS}: Spatial + \code{offset(log(hbs_avail))}, i.e. habitat-based
standardisation with the HBF-predicted vertical overlap as effective effort
(tuna model only; assumes the habitat model is known).
}
The latent \code{hook_depth} and realised \code{vertical_avail} are never used
as predictors.}
\examples{
\dontrun{
om <- simulate_tuna(seed = 42)
comparison <- standardize_tuna_indices(om)
comparison$diagnostics
plot_indices(comparison)
}
}
