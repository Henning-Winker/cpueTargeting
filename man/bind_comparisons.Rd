\name{bind_comparisons}
\alias{bind_comparisons}
\title{Combine several standardisation results into one}
\usage{
bind_comparisons(...)
}
\arguments{
\item{...}{Outputs of \code{\link{standardize_cpue}}.}
}
\value{An object like the output of \code{\link{standardize_cpue}}
with the \code{indices}, \code{diagnostics} and \code{formulas} of all inputs
(\code{data} and \code{pca} are taken from the first).}
\description{Useful to compare fits made with different settings (e.g.
spatial resolutions) in one plot. Give each call a different \code{label} so
that the model names are unique.}
\examples{
\dontrun{
a <- standardize_cpue(om, models = c("Base", "HBF+S"))
b <- standardize_cpue(om, models = "Base+S", spatial_res = 1,
                              label = " (1x1)")
plot_diagnostics(bind_comparisons(a, b))
}
}
