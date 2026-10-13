# Render the vignettes to docs/ (GitHub Pages) from pre-run results ------------------
# Run the simulation scripts first; the vignettes load their results from
# ../rdata. Run from the package folder:
#   source("inst/scripts/render_docs.R")

options(cpueTargeting.rdata = normalizePath("../rdata", mustWork = TRUE))
vignettes <- c(bathymetric = "vignettes/bathymetric-example.Rmd",
               flat_vs_contrast = "vignettes/tuna-flat-vs-contrast.Rmd")
for (v in vignettes) {
  rmarkdown::render(v, output_dir = "docs", envir = new.env(),
                    output_options = list(self_contained = TRUE))
}
