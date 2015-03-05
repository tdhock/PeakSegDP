\name{regularized.interval.regression}
\alias{regularized.interval.regression}
\title{regularized interval regression}
\description{Filter zero-variance features, scale features, filter flat limits,
then perform a path of increasingly regularized interval
regressions, returning a results list. We start at a small
regularization parameter specified as gamma.initial. Then we use
warm restarts, i.e. once we find the solution for one gamma, we
use that as a starting point for the optimization problem with a
larger gamma. We stop when gamma is so large that all the
coefficients are 0 at the optimum.}
\usage{regularized.interval.regression(features, limits, gamma.initial = 1e-04, 
    gamma.factor = 1.5, ...)}
\arguments{
  \item{features}{Matrix n x p of inputs: n signals, each with p features that will
be scaled. Zero-variance features will be ignored.}
  \item{limits}{Matrix n x 2 of output log(lambda). Each row corresponds to the
lower and upper bound of an interval on the lambda which is
optimal with respect to annotation error. Lower bound can be -Inf
and upper bound can be Inf, which correspond to zero asymptotic
cost.}
  \item{gamma.initial}{First regularization parameter in the path.}
  \item{gamma.factor}{Multiplicative factor to increase gamma between steps in the path.}
  \item{\dots}{Passed to optimization code in find.solution. You must specify
calc.grad and calc.loss.}
}

\value{List of solver results. For a feature matrix X with p columns, you
can use list$predict(X) to get model estimates of log(lambda).}

\author{Toby Dylan Hocking, Guillem Rigaill}





