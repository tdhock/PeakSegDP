\name{smooth.interval.regression}
\alias{smooth.interval.regression}
\title{smooth interval regression}
\description{Scale features and filter flat limits, then perform one interval
regression and return a results list. The precise optimization
problem is described in find.solution.}
\usage{smooth.interval.regression(features, limits, ...)}
\arguments{
  \item{features}{Matrix n x p of inputs: n signals, each with p features that will
be scaled.}
  \item{limits}{Matrix n x 2 of output log(lambda). Each row corresponds to the
lower and upper bound of an interval on the lambda which is
optimal with respect to annotation error. Lower bound can be -Inf
and upper bound can be Inf, which correspond to zero asymptotic
cost.}
  \item{\dots}{Passed to optimization code in find.solution. You must specify
calc.grad and calc.loss.}
}

\value{List of solver results. For a feature matrix X with p columns, you
can use list$predict(X) to get model estimates of lambda.}

\author{Toby Dylan Hocking, Guillem Rigaill}





