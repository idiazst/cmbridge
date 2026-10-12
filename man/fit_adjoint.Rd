\name{fit_adjoint}
\alias{fit_adjoint}
\title{Fit an adjoint representer}
\usage{fit_adjoint(B, V, M, phi, method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list())}
\arguments{
\item{B}{Variables indexing lambda.}
\item{V}{Variables conditioning the adjoint equation.}
\item{M}{Visit/response indicator.}
\item{phi}{Target-specific loading.}
\item{method}{Estimation method.}
\item{control}{Method-specific control list.}
}
\description{Fits lambda satisfying E[M(lambda(B)-phi(V)) given V]=0 using its complete-case conditional-moment form.}
\details{Optional \code{control$penalty_selection} chooses \code{"gcv"},
\code{"lcurve"}, or \code{"quasiopt"}; \code{penalty_selection_control}
passes path controls to \code{\link{select_cm_penalty}}. The choice is made
using only rows supplied to this fit.}
\value{An object of class \code{cmbridge_fit}.}
