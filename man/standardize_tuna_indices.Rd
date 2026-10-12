\name{standardize_tuna_indices}
\alias{standardize_tuna_indices}
\title{Standardise simulated CPUE with alternative targeting covariates}
\usage{
standardize_tuna_indices(simulation, models = NULL, spatial_res = 5,
  spatial_k = NULL, add_models = NULL, vessel_re = FALSE, n_pc = 2,
  n_clusters = 3, seed = 1, engine = c("gam", "bam"), label = NULL)
}
\arguments{
\item{simulation}{Output of \code{\link{simulate_tuna}} or \code{\link{simulate_bathymetric}}.}
\item{models}{Models to fit, e.g. \code{c("Base", "Base+S", "DPC+S", "HBF+S")}.
Default: \code{Base}, plus \code{Base}, \code{Cluster}, \code{DPC}, \code{HBF},
\code{Cluster_HBF} and \code{DPC_HBF} with \code{+S} (the HBF models only for
the tuna model).}
\item{spatial_res}{Size of the reporting blocks in cells (5 = 5 x 5 blocks as
in the CPUE guidelines; 1 = full resolution).}
\item{spatial_k}{Basis dimension of the spatial smoother. By default it grows
with the number of blocks (8 for 100 blocks, at most 60).}
\item{add_models}{Optional named character vector of extra right-hand sides,
e.g. \code{c(HBF_period = "year + s(Xs, Ys, k = 8) + HBF_period")}. They can
use any column of the returned \code{data}, including \code{Xs}, \code{Ys},
\code{period} (10-year block as a factor) and \code{HBF_period} (HBF x period).}
\item{vessel_re}{Logical; add a random vessel effect \code{s(vessel, bs = "re")}.}
\item{n_pc}{Number of principal components used by the DPC models (default
2, as in the CPUE guidelines; at most the number of species).}
\item{n_clusters}{Number of fishing tactics for the \code{Cluster} model.}
\item{seed}{Seed for the k-means clustering.}
\item{engine}{\code{"gam"} (default) or \code{"bam"} with \code{discrete = TRUE},
which is much faster for large data sets.}
\item{label}{Optional text appended to every model name, e.g. \code{" (1x1)"},
to tell fits at different resolutions apart after \code{\link{bind_comparisons}}.}
}
\value{A list with \code{indices} (species, year, model, index with 95\% CI and
the true index), \code{diagnostics} (root mean squared error of the log index
by species and model), \code{formulas} (the right-hand side of each model) and
\code{data} (the set-level data used, with PCs, tactics and block centres).}
\description{Fits Tweedie GAMs (\code{mgcv::gam} with \code{tw()}) to each
species and compares the annual indices with the true biomass trend.}
\details{Each model is \code{year} plus a targeting term:
\itemize{
\item \code{Base}: no targeting term.
\item \code{DPC}: \code{s(PC1) + ... + s(PCn)} from the catch composition
(direct principal components), with \code{n = n_pc}.
\item \code{Cluster}: fishing tactic from k-means clustering of the catch
composition (as a factor).
\item \code{HBF}: HBF as a factor (tuna model only).
\item \code{Cluster_HBF}: Cluster + HBF (tuna model only).
\item \code{DPC_HBF}: DPC + HBF (tuna model only).
\item \code{HBS}: \code{offset(log(hbs_avail))}, a specialised habitat-based
standardisation: effort is weighted by the overlap between the expected hook
profile of the HBF and the species' depth profile under the average
thermocline (tuna model only; not fitted by default).
}
Append \code{+S} to any model name to add a spatial smoother
\code{s(Xs, Ys)}, e.g. \code{"DPC+S"}; \code{"Base+S"} is the spatial-only
model. \code{Xs}, \code{Ys} are the centres of
\code{spatial_res} x \code{spatial_res} blocks of cells. The latent
\code{hook_depth} and realised \code{vertical_avail} are never used as
predictors.}
\examples{
\dontrun{
om <- simulate_tuna(seed = 42)
comparison <- standardize_tuna_indices(om, models = c("Base", "Base+S",
                                        "DPC+S", "Cluster+S", "HBF+S"))
comparison$diagnostics
plot_indices(comparison)
}
}
