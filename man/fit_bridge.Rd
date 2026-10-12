\name{fit_bridge}
\alias{fit_bridge}
\title{Fit an inverse-probability bridge}
\usage{fit_bridge(B, V, M, method = c("sieve_md", "landweber", "pmmr", "saturated_l1"), control = list(),
  selection = c("fixed", "robins", "syrgkanis"), selection_control = list(), fold_id = NULL, seed = 1L)}
\arguments{
\item{B}{Shadow-side variables.}
\item{V}{Bridge-side variables; may be missing where M is zero.}
\item{M}{Visit/response indicator.}
\item{method}{Estimation method.}
\item{control}{Method-specific control list.}
\item{selection}{Fixed controls or experimental adaptive selection; see \code{fit_cm}.}
\item{selection_control}{Adaptive selection controls; see \code{fit_bridge_ensemble}.}
\item{fold_id}{Optional selection folds.}
\item{seed}{Nonnegative random seed.}
}
\description{Fits beta satisfying E[M beta(V)-1 given B]=0.}
\value{An object of class \code{cmbridge_fit}.}
