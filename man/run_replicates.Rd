\name{run_replicates}
\alias{run_replicates}
\title{Run scenarios over several seeds (replicates)}
\usage{
run_replicates(scenarios, seeds, models, extra = NULL, file = NULL, ...)
}
\arguments{
\item{scenarios}{Named list of settings lists (see \code{\link{tuna_settings}}).}
\item{seeds}{Integer vector of seeds.}
\item{models}{Models passed to \code{\link{standardize_tuna_indices}}.}
\item{extra}{Optional list of extra standardisation calls, each a list of
arguments for \code{\link{standardize_tuna_indices}} (e.g. a 1 x 1 control
with a \code{label}); their results are added with
\code{\link{bind_comparisons}}.}
\item{file}{Optional \code{.rds} file to save and resume from.}
\item{...}{Further arguments for \code{\link{standardize_tuna_indices}}
(e.g. \code{spatial_res}, \code{n_pc}, \code{engine}).}
}
\value{A list with \code{diagnostics} and \code{indices}, each with columns
\code{scenario} and \code{seed} added.}
\description{Simulates each scenario for each seed, standardises the data and
keeps the indices and accuracy diagnostics. Results are saved to \code{file}
after every fit, and seeds/scenarios already in \code{file} are skipped, so a
long (e.g. overnight) run can be interrupted and resumed.}
\examples{
\dontrun{
base <- tuna_settings()
scenarios <- list(Fixed = base,
                  Moving = modifyList(base, list(thermocline_annual_sd = 15,
                                                 thermocline_annual_field_sd = 10)))
reps <- run_replicates(scenarios, seeds = 1:10,
                       models = c("Base+S", "Cluster+S", "DPC+S", "HBF+S"),
                       engine = "bam", file = "replicates.rds")
plot_replicates(reps)
}
}
