PoissonSeg <- structure(function
### A dynamic programming algorithm can be used to compute the best
### segmentation with respect to the Poisson likelihood, subject to a
### constraint on the number of segments, and the changes which must
### alternate: up, down, up, down, ...
(count,
### Integer vector of count data to segment.
 weight=rep(1, length(count)),
### Data weights (normally this is the number of base pairs).
 maxSegments
### Maximum number of segments to consider.
 ){
  stopifnot(is.numeric(count))
  stopifnot(is.numeric(weight))
  stopifnot(length(count) == length(weight))
  stopifnot(length(maxSegments) == 1)
  nData <- length(count)
  A <- .C("cPoisLinProgDyn",
          count=as.integer(count),
          weight=as.integer(weight),
          nData=as.integer(nData),
          maxSegments=as.integer(maxSegments),
          res1=double(maxSegments*nData),
          res2=integer(maxSegments*nData),
          res3=double(maxSegments*nData),
          PACKAGE="PeakSegDP")
  A$res1 <- matrix(A$res1, nrow=maxSegments, byrow=TRUE)
  A$res2 <- matrix(A$res2, nrow=maxSegments, byrow=TRUE)
  A$res3 <- matrix(A$res3, nrow=maxSegments, byrow=TRUE)
  return(A)
}, ex=function(){
  fit <- PoissonSeg(c(0, 10, 11, 1), maxSegments=3)
  stopifnot(fit$res2[3,4] == 3)
  stopifnot(fit$res2[2,3] == 1)
})

### Compute the weighted Poisson loss function, which is seg.mean -
### count * log(seg.mean). The edge case is when the mean is zero, in
### which case the probability mass function takes a value of 1 when
### the data is 0 (and 0 otherwise). Thus the log-likelihood of a
### maximum likelihood segment with mean zero must be zero.
PoissonLoss <- structure(function(count, seg.mean, weight=1){
  stopifnot(is.numeric(count))
  stopifnot(is.numeric(seg.mean))
  stopifnot(is.numeric(weight))
  n.data <- length(count)
  if(length(seg.mean) == 1){
    seg.mean <- rep(seg.mean, n.data)
  }
  if(length(weight) == 1){
    weight <- rep(weight, n.data)
  }
  stopifnot(n.data == length(seg.mean))
  stopifnot(n.data == length(weight))
  if(any(weight < 0)){
    stop("PoissonLoss undefined for negative weight")
  }
  if(any(seg.mean < 0)){
    stop("PoissonLoss undefined for negative segment mean")
  }
  not.integer <- round(count) != count
  not.positive <- count < 0
  loss <-
    ifelse(not.integer|not.positive, Inf,
           ifelse(seg.mean == 0,
                  ifelse(count == 0, 0, Inf),
                  seg.mean - count * log(seg.mean)))
  sum(loss*weight)
}, ex=function(){
  PoissonLoss(1, 1)
  PoissonLoss(0, 0)
  PoissonLoss(1, 0)
  PoissonLoss(0, 1)
})

### Extract endpoint matrix from PoissonSeg result.
getPath <- function(A){
  n <- ncol(A$res2)
  res3 <- matrix(NA, nrow=nrow(A$res2), ncol=nrow(A$res2))
  res3[1, 1] <- 0
  for(i in 2: nrow(A$res2)){
    res3[i, i-1] <- A$res2[i, n]
    if(res3[i, i-1] > 0){
      for(k in 1:(i-1)){
        res3[i, i-1-k] <- A$res2[i-k, res3[i, i-k]]
      }
    }
  }
  diag(res3) <- ncol(A$res2)
  return(res3)
}

PeakSegDP <- structure(function
### Compute the PeakSeg model on a data.frame of compressed sequence
### reads.
(compressed,
### data.frame with columns chromStart, chromEnd, count.
 maxPeaks
### maximum number of peaks to consider.
 ){
  stopifnot(diff(compressed$chromStart) > 0)
  count <- compressed$count
  weight <- compressed$bases <- with(compressed, chromEnd-chromStart)
  stopifnot(is.integer(count))
  stopifnot(is.integer(weight))
  stopifnot(count >= 0)
  stopifnot(weight > 0)
  stopifnot(is.integer(maxPeaks))
  stopifnot(length(maxPeaks) == 1)
  stopifnot(length(count) == length(weight))
  maxSegments <- maxPeaks * 2 + 1
  stopifnot(maxSegments > 0)
  stopifnot(maxSegments <= nrow(compressed))
  fit <- PoissonSeg(count, weight, maxSegments)
  segment.ends <- getPath(fit)
  results <- list()
  for(peaks in 0:maxPeaks){
    peak.list <- list()
    segments <- as.integer(peaks*2 + 1)
    model.i <- peaks * 2 + 1
    last.i <- as.integer(segment.ends[model.i, 1:model.i])
    break.after <- last.i[-model.i]
    first.i <- as.integer(c(1, break.after+1))
    model.error <- 0
    if(length(break.after)){
      results$breaks[[paste(peaks)]] <-
        data.frame(peaks, segments, break.after,
                   chromEnd=compressed$chromEnd[break.after])
    }
    for(segment.i in seq_along(last.i)){
      status <- ifelse(segment.i %% 2, "background", "peak")
      first <- first.i[[segment.i]]
      last <- last.i[[segment.i]]
      seg.data <- compressed[first:last,]
      seg.mean <- with(seg.data, sum(count*bases)/sum(bases))
      model.error <- model.error + with(seg.data, {
        PoissonLoss(count, seg.mean, bases)
      })
      chromStart <- seg.data$chromStart[1]
      chromEnd <- seg.data$chromEnd[nrow(seg.data)]
      results$segments[[paste(peaks, segment.i)]] <- 
        data.frame(mean=seg.mean,
                   first,
                   last,
                   chromStart,
                   chromEnd,
                   status,
                   peaks,
                   segments)
      if(status == "peak"){
        peak.list[[paste(segment.i)]] <-
          data.frame(first, last,
                     chromStart, chromEnd,
                     peaks, segments)
      }
    }#segment.i
    results$peaks[[as.character(peaks)]] <- do.call(rbind, peak.list)
    results$error[[as.character(peaks)]] <- 
      data.frame(segments=model.i, peaks, error=model.error)
  }
  results$peaks <- c(list("0"=dp.fit$peaks[[1]][0,]), results$peaks)
  results$error <- do.call(rbind, results$error)
  results$segments <- do.call(rbind, results$segments)
  results$breaks <- do.call(rbind, results$breaks)
  results
}, ex=function(){
  data(chr11ChIPseq)
  one <- subset(chr11ChIPseq$coverage, sample.id=="McGill0002")
  fit <- PeakSegDP(one, 3L)
  library(ggplot2)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=one)+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 data=fit$segments, color="green")+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 data=subset(fit$segments, status=="peak"),
                 size=3, color="deepskyblue")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(peaks ~ ., scales="free", labeller=function(var, val){
      s <- ifelse(val==1, "", "s")
      paste0(val, " peak", s)
    })
})
