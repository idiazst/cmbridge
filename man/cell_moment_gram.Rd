\name{cell_moment_gram}
\alias{cell_moment_gram}
\title{Cell conditional-moment U-statistic}
\usage{cell_moment_gram(residual, instrument)}
\arguments{
\item{residual}{Numeric vector or matrix, one candidate per column.}
\item{instrument}{Numeric conditioning variables; every column is retained.}
}
\description{Computes within-cell ordered cross-products over distinct observations, divided by the validation sample size times one less than the cell count. Singleton cells contribute zero. If there is no repeated cell, scoring fails.}
\value{A symmetric matrix, possibly indefinite, with a \code{cell_scoring} attribute recording cell and singleton counts.}
