multiSampleSegHeuristic <- structure(function
### Find one peak common to several samples.
(profiles,
### List of data.frames with columns chromStart, chromEnd, count, or
### single data.frame with additional column sample.id.
 bin.factor=2L
### Size of bin pyramid.
 ){
  stopifnot(is.numeric(bin.factor))
  stopifnot(length(bin.factor)==1)
  if(is.data.frame(profiles)){
    profiles <- split(profiles, profiles$sample.id, drop=TRUE)
  }
  stopifnot(is.list(profiles))
  for(profile.i in seq_along(profiles)){
    df <- profiles[[profile.i]]
    stopifnot(is.data.frame(df))
    stopifnot(is.integer(df$chromStart))
    stopifnot(is.integer(df$chromEnd))
    stopifnot(is.integer(df$count))
    profiles[[profile.i]] <- df[, c("chromStart", "chromEnd", "count")]
  }
  chromStartEnd <-
    .Call("multiSampleSegHeuristic_interface",
          profiles,
          as.integer(bin.factor),
          PACKAGE="PeakSegDP")
  data.frame(chromStart=chromStartEnd[1],
             chromEnd=chromStartEnd[2])
}, ex=function(){
  library(PeakSegDP)
  data(chr11ChIPseq)
  two <- subset(chr11ChIPseq$coverage,
                118090000 < chromStart &
                chromEnd < 118100000 &
                sample.id %in% c("McGill0002", "McGill0004"))
  ## Find the best peak location across 2 samples.

  ## optimal.seconds <- system.time({
  ##   optimal <- multiSampleSegOptimal(two)
  ## })[["elapsed"]]
  optimal <- data.frame(chromStart=NA, chromEnd=NA)
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(two, 2)
  })[["elapsed"]]
  rbind(heuristic.seconds, optimal.seconds)
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  library(ggplot2)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=two)+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_segment(aes(chromStart/1e3, 0,
                     color=model, size=model,
                     xend=chromEnd/1e3, yend=0),
                 data=peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")

  ## Heuristic in R code for debugging.
  two.list <- split(two, two$sample.id, drop=TRUE)
  max.chromStart <- max(sapply(two.list, with, chromStart[1]))
  min.chromEnd <- min(sapply(two.list, with, chromEnd[length(chromEnd)]))
  bases <- min.chromEnd-max.chromStart
  bin.factor <- 2L
  bases.per.bin <- 1L
  while(bases/bases.per.bin/bin.factor >= 4){
    bases.per.bin <- bases.per.bin * bin.factor
  }
  n.bins <- as.integer(floor(bases/bases.per.bin))
  n.samples <- length(two.list)
  cumsum.mat <- bin.mat <- matrix(NA, n.bins, n.samples)
  bin.list <- list()
  norm.list <- list()
  for(sample.i in seq_along(two.list)){
    sample.id <- names(two.list)[sample.i]
    one <- two.list[[sample.i]]
    max.count <- max(one$count)
    bins <- binSum(one, max.chromStart, bases.per.bin, n.bins)
    bins$mean.norm <- bins$mean/max.count
    bin.list[[sample.id]] <- data.frame(sample.id, rbind(bins, NA))
    bin.mat[, sample.i] <- bins$count
    cumsum.mat[, sample.i] <- cumsum(bins$count)
    one$count.norm <- one$count/max.count
    norm.list[[sample.i]] <- one
  }
  bin.df <- do.call(rbind, bin.list)
  norm.df <- do.call(rbind, norm.list)
  OptimalPoissonLoss <- function(cumsum.value, mean.value){
    ifelse(cumsum.value == 0, 0, cumsum.value * (1-log(mean.value)))
  }
  loss.list <- list()
  for(seg1.last in 1:(n.bins-2)){
    seg1.cumsums <- cumsum.mat[seg1.last, ]
    seg1.means <- seg1.cumsums/seg1.last/bases.per.bin
    seg1.loss <- OptimalPoissonLoss(seg1.cumsums, seg1.means)
    for(seg2.last in (seg1.last+1):(n.bins-1)){
      seg2.cumsums <- cumsum.mat[cbind(seg2.last, 1:n.samples)]-seg1.cumsums
      seg2.bases <- seg2.last - seg1.last
      seg2.means <- seg2.cumsums/seg2.bases/bases.per.bin
      seg2.loss <- OptimalPoissonLoss(seg2.cumsums, seg2.means)
      seg3.cumsums <- cumsum.mat[cbind(n.bins, 1:n.samples)]-seg2.cumsums
      seg3.bases <- n.bins-seg2.last
      seg3.means <- seg3.cumsums/seg3.bases/bases.per.bin
      seg3.loss <- OptimalPoissonLoss(seg3.cumsums, seg3.means)
      total.loss <- seg1.loss + seg2.loss + seg3.loss
      loss.list[[paste(seg1.last, seg2.last)]] <-
        data.frame(seg1.last, seg2.last, total.loss,
                   peakStart=seg1.last*bases.per.bin+max.chromStart,
                   peakEnd=seg2.last*bases.per.bin+max.chromStart)
    }
  }
  loss.df <- do.call(rbind, loss.list)
  loss.best <- loss.df[which.min(loss.df$total.loss), ]
  
  ggplot()+
    scale_size_manual(values=c(data=2, bins=1))+
    scale_color_manual(values=c(data="grey50", bins="black", peak="green"))+
    geom_step(aes(chromStart/1e3, count.norm, color=what),
              data=data.frame(norm.df, what="data"))+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm,
                     color=what),
                 data=data.frame(bin.df, what="bins"))+
    geom_segment(aes(peakStart/1e3, 0,
                     color=what,
                     xend=peakEnd/1e3, yend=0),
                 data=data.frame(loss.best, what="peak"))+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ .)
  peakStart <- loss.best$peakStart
  peakEnd <- loss.best$peakEnd
  left.cumsum.vec <- if(loss.best$seg1.last == 1){
    rep(0, n.samples)
  }else{
    cumsum.mat[loss.best$seg1.last-1, ]
  }
  right.cumsum.vec <- cumsum.mat[loss.best$seg2.last-1, ]
  n.bins.zoom <- bin.factor * 2L
  n.cumsum <- n.bins.zoom + 1L
  while(bases.per.bin > 1){
    left.chromStart <- peakStart - bases.per.bin
    right.chromStart <- peakEnd-bases.per.bin
    bases.per.bin <- as.integer(bases.per.bin / bin.factor)

    right.cumsum.mat <- left.cumsum.mat <- matrix(NA, n.cumsum, n.samples)
    right.limits <- bases.per.bin*(0:(n.cumsum-1))+right.chromStart
    right.intervals <-
      paste0(right.limits[-length(right.limits)], "-", right.limits[-1])
    rownames(right.cumsum.mat) <- c("before", right.intervals)
    left.limits <- bases.per.bin*(0:(n.cumsum-1))+left.chromStart
    left.intervals <-
      paste0(left.limits[-length(left.limits)], "-", left.limits[-1])
    rownames(left.cumsum.mat) <- c("before", left.intervals)
    for(sample.i in seq_along(two.list)){
      one <- two.list[[sample.i]]
      left.bins <- binSum(one, left.chromStart, bases.per.bin, n.bins.zoom)
      right.bins <- binSum(one, right.chromStart, bases.per.bin, n.bins.zoom)
      left.count <- c(left.cumsum.vec[sample.i], left.bins$count)
      left.cumsum.mat[, sample.i] <- cumsum(left.count)
      right.count <- c(right.cumsum.vec[sample.i], right.bins$count)
      right.cumsum.mat[, sample.i] <- cumsum(right.count)
    }
  }
  four <- subset(chr11ChIPseq$coverage,
                 118120000 < chromStart &
                 chromEnd < 118126000) 
  ## Find the best peak location across 4 samples.
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(four)
  })[["elapsed"]]
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(four, 2)
  })[["elapsed"]]
  rbind(heuristic.seconds, optimal.seconds)
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=four)+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_segment(aes(chromStart/1e3, 0,
                     color=model, size=model,
                     xend=chromEnd/1e3, yend=0),
                 data=peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")

  ## A fake data set with two profiles with very different scales for
  ## the count variable, showing that the profile with the small
  ## counts will be basically ignored when computing the optimal
  ## peak.
  count <-
    c(rep(c(0, 1), each=100),
      rep(c(10, 11), each=200),
      rep(c(0, 1), each=100))
  chromEnd <- seq_along(count)
  chromStart <- chromEnd-1L
  offset <- 50L
  multi.scale <- 
  rbind(data.frame(sample.id="low", chromStart, chromEnd,
                   count=as.integer(count)),
        data.frame(sample.id="hi",
                   chromStart=chromStart+offset,
                   chromEnd=chromEnd+offset,
                   count=as.integer(count*1000)))
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(multi.scale, 2)
  })[["elapsed"]]
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(multi.scale)
  })[["elapsed"]]
  rbind(heuristic.seconds, optimal.seconds)
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  ggplot()+
    geom_step(aes(chromStart, count), data=multi.scale)+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_segment(aes(chromStart, 0,
                     color=model, size=model,
                     xend=chromEnd, yend=0),
                 data=peaks)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")
  

})

multiSampleSegOptimal <- function
### Find one peak common to several samples.
(profiles
### List of data.frames with columns chromStart, chromEnd, count, or
### single data.frame with additional column sample.id.
 ){
  if(is.data.frame(profiles)){
    profiles <- split(profiles, profiles$sample.id, drop=TRUE)
  }
  stopifnot(is.list(profiles))
  for(profile.i in seq_along(profiles)){
    df <- profiles[[profile.i]]
    stopifnot(is.data.frame(df))
    stopifnot(is.integer(df$chromStart))
    stopifnot(is.integer(df$chromEnd))
    stopifnot(is.integer(df$count))
    profiles[[profile.i]] <- df[, c("chromStart", "chromEnd", "count")]
  }
  chromStartEnd <-
    .Call("multiSampleSegOptimal_interface",
          profiles,
          PACKAGE="PeakSegDP")
  data.frame(chromStart=chromStartEnd[1],
             chromEnd=chromStartEnd[2])
}
