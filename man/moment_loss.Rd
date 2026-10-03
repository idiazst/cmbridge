\name{moment_loss}
\alias{moment_loss}
\title{Evaluate an empirical conditional-moment loss}
\usage{moment_loss(object, response = NULL, diagonal = NULL, target = NULL, instrument = NULL)}
\arguments{
\item{object}{A fitted \code{cmbridge_fit}.}
\item{response}{Optional new response.}
\item{diagonal}{Optional new diagonal.}
\item{target}{Optional new target variables.}
\item{instrument}{Optional new instrument variables.}
}
\description{Returns the training moment loss or an out-of-sample moment loss. Individual learners use their fitted instrument features and weighting. Ensembles use their common scoring kernel, with its full-training-sample scaling, bandwidth, and optional Nystr\"om centers held fixed.}
\value{A nonnegative numeric scalar.}
