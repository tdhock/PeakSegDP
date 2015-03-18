\name{H3K4me3.TDH.immune.chunk12.cluster4}
\alias{H3K4me3.TDH.immune.chunk12.cluster4}
\docType{data}
\title{
Histone ChIP-seq data, 26 samples, chr1 subset
}
\description{
26 samples, each with the same overlapping peak(s).
}
\usage{data("H3K4me3.TDH.immune.chunk12.cluster4")}
\format{
  A data frame with 20199 observations on the following 6 variables.
  \describe{
    \item{\code{cell.type}}{a factor with levels \code{monocyte} \code{bcell} \code{tcell}}
    \item{\code{sample.id}}{a factor with 26 levels}
    \item{\code{chromStart}}{integer vector}
    \item{\code{chromEnd}}{integer vector}
    \item{\code{count}}{integer vector}
  }
}
\source{
http://cbio.ensmp.fr/~thocking/chip-seq-chunk-db H3K4me3_TDH_immune data
set, chunk.id=12.
}
\keyword{datasets}
