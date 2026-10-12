\name{fit_bridge}
\alias{fit_bridge}
\title{Fit an inverse-probability bridge}
\usage{fit_bridge(B, V, M, method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list())}
\arguments{
\item{B}{Shadow-side variables.}
\item{V}{Bridge-side variables; may be missing where M is zero.}
\item{M}{Visit/response indicator.}
\item{method}{Estimation method.}
\item{control}{Method-specific control list.}
}
\description{Fits beta satisfying E[M beta(V)-1 given B]=0.}
\details{Optional \code{control$penalty_selection} chooses \code{"gcv"},
\code{"lcurve"}, or \code{"quasiopt"}; \code{penalty_selection_control}
passes path controls to \code{\link{select_cm_penalty}}. The choice is made
using only rows supplied to this fit.}
\value{An object of class \code{cmbridge_fit}.}
