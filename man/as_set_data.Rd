\name{as_set_data}
\alias{as_set_data}
\title{Reshape set-by-species records (long) to one row per set}
\usage{
as_set_data(data, value = "CPUE")
}
\arguments{
\item{data}{Long data with one row per set and species, e.g. \code{sim$Data}.}
\item{value}{Column to spread across species (default \code{"CPUE"}).}
}
\value{A data frame with set-level columns and one column per species.}
\description{Reshape set-by-species records (long) to one row per set.}
\examples{
\dontrun{
sim <- simulate_tuna(seed = 1)
sets <- as_set_data(sim$Data)
head(sets)
}
}
