\name{fit_adjoint}
\alias{fit_adjoint}
\title{Fit an adjoint representer}
\usage{fit_adjoint(B, V, M, phi, method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list(),
  selection = c("fixed", "robins", "syrgkanis"), selection_control = list(), fold_id = NULL, seed = 1L)}
\arguments{
\item{B}{Variables indexing lambda.}
\item{V}{Variables conditioning the adjoint equation.}
\item{M}{Visit/response indicator.}
\item{phi}{Target-specific loading.}
\item{method}{Estimation method.}
\item{control}{Method-specific control list.}
\item{selection}{Fixed controls or experimental adaptive selection; see \code{fit_cm}.}
\item{selection_control}{Adaptive selection controls; see \code{fit_bridge_ensemble}.}
\item{fold_id}{Optional folds on the full input rows; restricted to complete cases.}
\item{seed}{Nonnegative random seed.}
}
\description{Fits lambda satisfying E[M(lambda(B)-phi(V)) given V]=0 using its complete-case conditional-moment form.}
\value{An object of class \code{cmbridge_fit}.}
