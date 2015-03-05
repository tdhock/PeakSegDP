\name{hinge.interval.regression}
\alias{hinge.interval.regression}
\title{hinge interval regression}
\description{Support vector interval regression using a quadratic programming
(QP) solver. The idea is that we first normalize the feature
matrix, giving normalized features x_i in R^p. Then we use a
linear function f(x_i) = w'x_i + b to predict a log(lambda) that
falls between the log limits L_i^left and L_i^right. So the
optimization problem is: min_f ||f|| + sum_i C*hinge(L_i^left,
L_i^right, f(x_i)). Since we assume f is linear the problem
becomes min_{w,b,z} w'w + sum_i C*z_i, with these constraints for
all relevant i: z_i >= 0, z_i >= 1 + L_i^left - b - w'x_i, z_i >=
1 - b - w'x + L_i^right. We call z_i slack, w weights, b intercept.}
\usage{hinge.interval.regression(features, limits, tune.C = 1, verbose = 0, 
    ...)}
\arguments{
  \item{features}{Matrix n x p of inputs: n signals, each with p features. We will
scale these internally.}
  \item{limits}{Matrix n x 2 of output lambda. Each row corresponds to the lower
and upper bound of an interval on the log(lambda) which is optimal
with respect to annotation error. Lower bound can be -Inf and
upper bound can be Inf, which correspond to zero asymptotic
cost. }
  \item{tune.C}{
}
  \item{verbose}{
}
  \item{\dots}{ignored.}
}

\value{List of solver results. For a feature matrix X with p columns, you
can use list$predict(X) to get model estimates of log(lambda).}

\author{Toby Dylan Hocking, Guillem Rigaill}





