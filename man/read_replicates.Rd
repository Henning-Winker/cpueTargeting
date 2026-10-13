\name{read_replicates}
\alias{read_replicates}
\title{Read replicate results saved one file per seed}
\usage{
read_replicates(dir, pattern = "\\\\.rds$", label = NULL)
}
\arguments{
\item{dir}{Folder with the \code{.rds} files.}
\item{pattern}{Regular expression for the file names.}
\item{label}{Optional label stored in a column \code{label} (e.g.
\code{"Flat trends"}), to combine runs of different scenarios.}
}
\value{A list with \code{diagnostics} and \code{indices}, as from
\code{\link{run_replicates}}.}
\description{Combines the \code{.rds} files written by
\code{\link{run_replicates}} (e.g. one file per seed from a parallel run) into
one result for \code{\link{plot_replicates}}.}
