# cpueTargeting

`cpueTargeting` is a small R package for simulation experiments in multispecies CPUE
standardization. It brings together two related operating models:

1. a habitat-targeting trawl example based on the original 2014
   `SimulationCPUE_3d.R` script; and
2. a spatial longline tuna extension with separate horizontal habitat and
   vertical fish–hook overlap, including HBF as an imperfect hook-depth proxy.

The historical script is preserved unchanged in `inst/extdata`. The package
function `simulate_trawl()` provides a reproducible, callable version of its
data-generating scenario. `simulate_tuna()` runs the tuna operating model.
Figures returned by the package functions use ggplot2. The archived original
script is retained as-is and includes its original base-graphics code.

## Installation

After creating a GitHub repository for this folder, install with:

```r
remotes::install_github("YOUR-ACCOUNT/cpueTargeting")
```

Required packages are declared in `DESCRIPTION` and installed with the
package. Install `roxygen2` and run `roxygen2::roxygenise()` when changing
function documentation.

## Quick start

```r
library(cpueTargeting)

trawl <- simulate_trawl(seed = 1)
head(trawl$data)
trawl_plots <- plot_trawl_simulation(trawl)
trawl_plots$effort

tuna <- simulate_tuna(seed = 2)
head(tuna$Data)
comparison <- standardize_tuna_indices(tuna)

plot_tuna_vertical_distributions()
```

Use `tuna_settings()` to inspect or edit the illustrative tuna scenario.
Vertical profiles, HBF-to-depth settings, prices, and the spatial settings are
demonstration values that need fleet-specific calibration before scientific
interpretation. The operating model keeps latent `hook_depth_m` and vertical
overlap for validation; the analyst-facing gear predictor is HBF.

## Main functions

- `simulate_trawl()`: returns trip-level catches, true biomass, and effort
  allocation for the four-habitat, ten-species trawl scenario.
- `plot_trawl_simulation()`: returns ggplot figures for trawl biomass, effort,
  and catch distributions.
- `tuna_settings()`: returns default settings for the tuna operating model.
- `simulate_tuna()`: returns spatial and vertical tuna simulation output.
- `standardize_tuna_indices()`: compares the base, spatial, DPC, cluster and
  HBF-adjusted CPUE indices against simulated biomass.
- `plot_tuna_vertical_distributions()`: plots species availability and HBF
  hook-depth profiles on a common depth axis.

## Development

The package uses roxygen2 comments in `R/` to generate help pages and the
namespace. The standalone workshop driver remains available in the parent
`project_sources` directory as a reference for the CPUE index comparisons; the
first package version focuses on clean, documented operating-model functions.
