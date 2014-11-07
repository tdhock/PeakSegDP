\name{largestContinuousMinimum}
\alias{largestContinuousMinimum}
\title{largestContinuousMinimum}
\description{Find the run of minimum cost with the largest size.}
\usage{largestContinuousMinimum(cost, size)}
\arguments{
  \item{cost}{
}
  \item{size}{
}
}



\author{Toby Dylan Hocking, Guillem Rigaill}




\examples{
data(chr11ChIPseq)
one <- subset(chr11ChIPseq$coverage, sample.id=="McGill0322")
fit <- PeakSegDP(one, 5L)
## First compute the optimal number of peaks function.
exact.df <- with(fit$error, exactModelSelection(error, peaks))
## Then compute the PeakError of these models with respect to the
## annotated regions.
regions <- subset(chr11ChIPseq$regions, sample.id=="McGill0322")
peak.list <- split(fit$segments, fit$segments$peaks)
require(PeakError)
all.error <- NULL
for(peaks.chr in names(peak.list)){
  peaks <- subset(peak.list[[peaks.chr]], status=="peak")
  error <- PeakErrorChrom(peaks, regions)
  all.error <- rbind(all.error, {
    data.frame(peaks=as.integer(peaks.chr),
               errors=with(error, sum(fp+fn)))
  })
}
rownames(all.error) <- all.error$peaks
exact.df$errors <-
  all.error[as.character(exact.df$model.complexity), "errors"]
indices <- with(exact.df, {
  largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
})
## The target interval (min.log.lambda, max.log.lambda) is the
## largest continuous interval such that the error is minimal.
target.interval <- with(indices, {
  data.frame(min.log.lambda=exact.df$min.log.lambda[start],
             max.log.lambda=exact.df$max.log.lambda[end])
})
print(target.interval)
}
