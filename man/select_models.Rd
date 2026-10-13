\name{select_models}
\alias{select_models}
\title{Keep only some models of a standardisation result}
\usage{
select_models(comparison, models)
}
\arguments{
\item{comparison}{Output of \code{\link{standardize_cpue}} or
\code{\link{bind_comparisons}}.}
\item{models}{Model names to keep.}
}
\value{The same structure, restricted to \code{models}.}
\description{Useful to reuse an expensive fit (e.g. a 1 x 1 spatial control)
and attach it to new fits with \code{\link{bind_comparisons}} instead of
refitting it.}
\examples{
\dontrun{
control <- select_models(old_fit, "Base+S (1x1)")
plot_indices(bind_comparisons(new_fit, control))
}
}
