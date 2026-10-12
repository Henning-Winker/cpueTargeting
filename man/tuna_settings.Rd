\name{tuna_settings}
\alias{tuna_settings}
\title{Default settings for the longline tuna scenario}
\usage{
tuna_settings(preset = c("default", "stress_test"))
}
\arguments{
\item{preset}{\code{"default"} or \code{"stress_test"}.}
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
one deep (BET, BFT). A thermocline-depth anomaly shifts all profiles up or
down together in each cell, so the catch of a given gear configuration varies
in space and, optionally, between years. The anomaly in year t and cell i is
\code{fixed[i] + shift[t] + annual[i, t]}: \code{fixed} is a spatial field that
stays the same every year (\code{thermocline_sd}; 0 = same vertical profile
everywhere); \code{shift} is a basin-wide AR(1) shift
(\code{thermocline_annual_sd}, \code{thermocline_annual_rho}), e.g. ENSO-like
years; \code{annual} is a new spatial field each year
(\code{thermocline_annual_field_sd}).
\item Hook depth: each HBF (hooks between floats) level has a typical
(median) hook depth \code{HBF_depth}, roughly as in Fig. 2 of Bigelow et al.
(2002). Each set's hook profile is shifted by a lognormal set-level error
(\code{set_cv}; e.g. line angle, currents), so HBF is an imperfect proxy for
depth. Hooks within a set are lognormally spread around the set's median
(\code{hook_cv}).
\item Targeting: species prices change in three 10-year blocks: ALB
(shallow, cool water), then BET (deep, warm water), then BFT (deep, cool
water), loosely following the historical shift from albacore to bigeye
targeting and the later interest in recovering southern bluefin. Skippers
choose where to fish and which HBF to use from the expected catch value.
\item Effort increases towards 800 sets per year and is cut by 25\% in
years 21-30, so pressure on some species is relieved.
}
All values are conceptual and meant to be tuned. With the same seed,
scenarios that differ only in the thermocline settings share the same SST
field, the same fixed thermocline pattern (scaled by \code{thermocline_sd})
and the same vessels; the fishery itself then evolves differently.}
\section{Presets}{
\itemize{
\item \code{"default"}: moderate targeting (prices differ by a factor of 2
between the target species and the others, \code{preference = 1.25},
\code{HBF_choice_temp = 0.2}); nominal CPUE is biased by roughly 10-20\% by
decade.
\item \code{"stress_test"}: extreme targeting (prices differ by a factor of 4,
\code{preference = 1.5}, \code{HBF_choice_temp = 0.15}); abrupt, fleet-wide
target switches with nominal CPUE biased by up to a factor of two. These are
the settings of the first runs of the vignette.
}
}
\examples{
s <- tuna_settings()
extreme <- tuna_settings("stress_test")
extreme$price[c(1, 11, 21), ]
}
