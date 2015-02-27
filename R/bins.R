binSum <- function
### Compute sum of compressed coverage profile in bins, using fast C
### code.
(compressed,
### data.frame with integer columns chromStart, chromEnd, count.
 bin.chromStart=0L,
### Base before first bin.
 bin.size=1L,
### Bin size.
 n.bins=2000L
### Number of bins.
 ){
  stopifnot(is.integer(compressed$chromStart))
  stopifnot(is.integer(compressed$chromEnd))
  stopifnot(is.integer(compressed$count))
  stopifnot(is.integer(bin.chromStart))
  stopifnot(length(bin.chromStart) == 1)
  stopifnot(is.integer(bin.size))
  stopifnot(length(bin.size) == 1)
  stopifnot(is.integer(n.bins))
  stopifnot(length(n.bins) == 1)
  result <- 
  .C("binSum_interface",
     profile.chromStart=as.integer(compressed$chromStart),
     profile.chromEnd=as.integer(compressed$chromEnd),
     profile.coverage=as.integer(compressed$count),
     n.profiles=as.integer(nrow(compressed)),
     bin.total=integer(n.bins),
     bin.size=as.integer(bin.size),
     n.bins=as.integer(n.bins),
     bin.chromStart=as.integer(bin.chromStart),
     package="PeakSegDP")
  total <- result$bin.total
  chromStart <- seq(bin.chromStart, by=bin.size, l=n.bins)
  chromEnd <- chromStart + bin.size
  data.frame(chromStart, chromEnd, total)
### data.frame with n.bins rows and columns chromStart, chromEnd,
### total.
}
