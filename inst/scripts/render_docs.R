# Render the vignettes to docs/ (GitHub Pages) from pre-run results ------------------
# Run the simulation scripts first; the vignettes load their results from
# rdata/ in the repository (not pushed). Run from the package folder:
#   source("inst/scripts/render_docs.R")

# Same folder as used by the simulation scripts (set the option before
# sourcing if the results are elsewhere, e.g. in your presentation folder)
options(cpueTargeting.rdata = normalizePath(getOption("cpueTargeting.rdata", "rdata"),
                                            mustWork = TRUE))
vignettes <- c(bathymetric = "vignettes/bathymetric-example.Rmd",
               flat_vs_contrast = "vignettes/tuna-flat-vs-contrast.Rmd",
               # these simulate their own data (no pre-run results needed)
               real_data = "vignettes/real-data-guide.Rmd",
               logbook_om = "vignettes/logbook-om-guide.Rmd",
               logbook = "vignettes/logbook-example.Rmd")
for (v in vignettes) {
  rmarkdown::render(v, output_dir = "docs", envir = new.env(),
                    output_options = list(self_contained = TRUE))
}
