\name{tuna_settings}
\alias{tuna_settings}
\title{Default settings for the longline tuna scenario}
\usage{
tuna_settings()
}
\value{A named list accepted by \code{\link{simulate_tuna}}.}
\description{Four tuna species in a conceptual Southern Hemisphere ocean
(y points north, towards the equator, so the north is warm).}
\details{
\itemize{
\item Horizontal habitat is a normal preference for sea surface temperature
(\code{sst_opt}, \code{sst_tol}). Yellowfin (YFT) and bigeye (BET) prefer warm
water; southern bluefin (BFT) and albacore (ALB) prefer cooler water. Wide
tolerances keep the species mixed rather than separated.
\item Vertical habitat is a normal depth profile (\code{depth_mu},
\code{depth_sd}, in m). Within each pair one species is shallow (YFT, ALB) and
one deep (BET, BFT). A thermocline-depth anomaly field (\code{thermocline_sd})
shifts all profiles up or down together in each cell, so the catch of a given
gear configuration varies in space. Set \code{thermocline_sd = 0} for the same
vertical profile everywhere.
\item Hook depth follows the catenary model of Bigelow et al. (2002) for each
HBF (hooks between floats), with between-set variation in the mainline angle
\code{phi}, within-set hook-depth variation and shoaling by currents.
\item Targeting: species prices change in three 10-year blocks (YFT, then
BET, then BFT). Skippers choose where to fish and which HBF to use from the
expected catch value.
}
All values are conceptual and meant to be tuned.}
