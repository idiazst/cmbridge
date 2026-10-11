\name{cell_moment_gram}
\alias{cell_moment_gram}
\title{Cell conditional-moment U- or V-statistic}
\usage{cell_moment_gram(residual, instrument, statistic = c("u", "v"))}
\arguments{
\item{residual}{Numeric vector or matrix, one candidate per column.}
\item{instrument}{Numeric conditioning variables; every column is retained.}
\item{statistic}{\code{"u"} (default) excludes self-products; \code{"v"} retains them.}
}
\description{U computes within-cell ordered cross-products over distinct observations, divided by the validation sample size times one less than the cell count. Singleton cells contribute zero, including a sample with no repeated cells. V includes self-products and divides by the validation sample size times the cell count. Singleton cells contribute their residual outer product divided by the sample size.}
\value{A symmetric matrix, PSD for V but possibly indefinite for U, with a \code{cell_scoring} attribute recording cell and singleton counts.}
