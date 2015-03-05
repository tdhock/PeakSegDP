\name{find.solution}
\alias{find.solution}
\title{find solution}
\description{Interval regression using a smooth loss to relax the annotation
loss, and an accelerated proximal gradient descent solver. The
idea is that we first normalize the feature matrix, giving
normalized features x_i in R^p. Then we use a linear function
f(x_i) = w'x_i + b to predict a log(lambda) that falls between the
log limits L_i^left and L_i^right. So the optimization problem is:
min_{w,b} gamma*||w||_1 + 1/n * sum_i phi[f(x_i)-L_i^left] +
phi[L_i^right-f(x_i)], where phi(L) is a convex relaxation of the
annotation loss, and should be specified in the calc.grad and
calc.loss arguments. The optimization stops when we find an
optimization variable for which each dimension is within a
threshold of the subgradient optimality condition.}
\usage{find.solution(features, limits, gamma = 0, starting.iterate = rep(0, 
    l = ncol(features) + 1), threshold = 0.001, verbose = 1, 
    max.iterations = 10000, calc.grad = stop("must specify calc.grad(x,features,limits)"), 
    calc.loss = stop("must specify calc.loss(x,features,limits)"), 
    L0 = ncol(features) + sqrt(ncol(features)), step = "constant")}
\arguments{
  \item{features}{feature matrix n x p which we assume has already been scaled.}
  \item{limits}{limit matrix n x 2.}
  \item{gamma}{regularization >= 0, by default 0.}
  \item{starting.iterate}{Where to start the optimization, by default at the origin.}
  \item{threshold}{When the stopping criterion gets below this threshold, the
solution is optimal.}
  \item{verbose}{print optimization status?}
  \item{max.iterations}{exit with an error if we haven't converged after this many
iterations.}
  \item{calc.grad}{Function x,features,limits => gradient vector. This will be used
for the optimization.}
  \item{calc.loss}{Function to display the calculate cost function, necessary for the
backtracking line search.}
  \item{L0}{Lipshitz constant for step size.}
  \item{step}{constant or linesearch.}
}



\author{Toby Dylan Hocking, Guillem Rigaill}





