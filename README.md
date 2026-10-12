# cpueTargeting

`cpueTargeting` is an R package for simulation experiments on fishing effort,
spatial targeting, and CPUE standardization. It keeps three related examples
separate:

1. `simulate_cluster()` is a callable version of the multi-habitat
   cluster simulator from Winker et al. (2014), based on
   `SimulationCPUE_3d.R`.
2. `simulate_bathymetric()` is the four-species Run 22 spatial operating model.
   Species distributions are associated with bathymetry; it does not include
   vertical availability.
3. `simulate_tuna()` extends the spatial operating model with vertical fish
   availability and a noisy HBF-to-hook-depth relationship.

The original Run 22 and 2014 scripts are preserved in `inst/extdata` for
reference. All plots produced by the exported package functions are ggplot2
objects, and the package code does not use or depend on `plyr`.

## Installation

```r
remotes::install_github("henning-winker/cpueTargeting")
```

Required packages are declared in `DESCRIPTION` and installed with the
package. Install `roxygen2` and run `roxygen2::roxygenise()` when changing
function documentation.

## Quick start

```r
library(cpueTargeting)

# Winker et al. (2014) multi-habitat cluster example
cluster <- simulate_cluster(seed = 1)
cluster_plots <- plot_cluster_2014(cluster)
cluster_plots$effort

# Four-species Run 22 model with bathymetric spatial distributions
bathymetric <- simulate_bathymetric(seed = 2)
head(bathymetric$Data)
bathymetric$Plots$SpDens

# Tuna extension with explicit vertical availability and HBF proxy
tuna <- simulate_tuna(seed = 3)
head(tuna$Data)
comparison <- standardize_tuna_indices(tuna)
plot_tuna_vertical_distributions()
```

Use `bathymetric_settings()` and `tuna_settings()` to inspect or edit the
respective scenarios. Tuna vertical profiles and HBF-to-depth settings are
illustrative placeholders that need fishery-specific calibration before
scientific interpretation. The tuna operating model retains latent
`hook_depth_m` and vertical overlap for validation; the analyst-facing gear
predictor is HBF.

## Main functions

- `simulate_cluster()` and `plot_cluster_2014()`: simulate and plot the
  Winker et al. (2014) multi-habitat cluster example.
- `bathymetric_settings()` and `simulate_bathymetric()`: configure and run the
  four-species bathymetric Run 22 model.
- `tuna_settings()` and `simulate_tuna()`: configure and run the spatial and
  vertical tuna extension.
- `standardize_tuna_indices()`: compare base, spatial, DPC, cluster and
  HBF-adjusted CPUE indices against simulated tuna biomass.
- `plot_tuna_vertical_distributions()`: plot species availability and HBF
  hook-depth profiles on a common depth axis.
