\name{select_penalty_grid}
\alias{select_penalty_grid}
\title{Select a positive penalty on an extendable grid}
\usage{select_penalty_grid(scales, evaluate, max_extensions = 16L,
  extension_factor = 10, extension_points = 4L)}
\arguments{
\item{scales}{Positive initial penalty scales.}
\item{evaluate}{Function returning losses for newly requested scales.}
\item{max_extensions}{Maximum extensions before an unresolved boundary fails.}
\item{extension_factor}{Multiplicative width of each extension.}
\item{extension_points}{New logarithmically spaced values per extended side.}
}
\description{Records a boundary minimum as a failed attempt, extends that side of the grid, and accepts only an interior minimum. An exact tie can select an interior minimizer. No zero penalty is used.}
\value{A list with \code{scale}, \code{penalty_cv}, \code{boundary_history}, and \code{status}. An unresolved boundary raises \code{cmbridge_penalty_boundary_error} carrying the evaluated grid and attempt history.}
