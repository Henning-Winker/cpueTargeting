\name{bind_comparisons}
\alias{bind_comparisons}
\title{Combine several standardisation results into one}
\usage{
bind_comparisons(...)
}
\arguments{
\item{...}{Outputs of \code{\link{standardize_tuna_indices}}.}
}
\value{An object like the output of \code{\link{standardize_tuna_indices}}
with the \code{indices}, \code{diagnostics} and \code{formulas} of all inputs
(\code{data} is taken from the first).}
\description{Useful to compare fits made with different settings (e.g.
spatial resolutions) in one plot. Give each call a different \code{label} so
that the model names are unique.}
\examples{
\dontrun{
a <- standardize_tuna_indices(om, models = c("Base", "HBF+S"))
b <- standardize_tuna_indices(om, models = "Base+S", spatial_res = 1,
                              label = " (1x1)")
plot_diagnostics(bind_comparisons(a, b))
}
}
