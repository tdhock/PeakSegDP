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
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(two)
  })[["elapsed"]]
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
      seg2.bases <- (seg2.last - seg1.last)*bases.per.bin
      seg2.means <- seg2.cumsums/seg2.bases
      seg2.loss <- OptimalPoissonLoss(seg2.cumsums, seg2.means)
      seg3.cumsums <- cumsum.mat[cbind(n.bins, 1:n.samples)]-seg2.cumsums
      seg3.bases <- (n.bins-seg2.last)*bases.per.bin
      seg3.means <- seg3.cumsums/seg3.bases
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

  ## These need to be computed once and then they never change.
  last.cumsum.vec <- cumsum.mat[nrow(cumsum.mat), ]
  last.chromEnd <- n.bins * bases.per.bin + max.chromStart
  n.bins.zoom <- bin.factor * 2L
  n.cumsum <- n.bins.zoom + 1L

  ## These will change at the end of each iteration.
  peakStart <- loss.best$peakStart
  peakEnd <- loss.best$peakEnd
  left.cumsum.vec <- if(loss.best$seg1.last == 1){
    rep(0, n.samples)
  }else{
    cumsum.mat[loss.best$seg1.last-1, ]
  }
  right.cumsum.vec <- cumsum.mat[loss.best$seg2.last-1, ]
  while(bases.per.bin > 1){
    left.chromStart <- peakStart - bases.per.bin
    right.chromStart <- peakEnd-bases.per.bin
    bases.per.bin <- as.integer(bases.per.bin / bin.factor)

    right.cumsum.mat <- left.cumsum.mat <- matrix(NA, n.cumsum, n.samples)
    right.limits <- bases.per.bin*(0:(n.cumsum-1))+right.chromStart
    right.chromStart.vec <- right.limits[-length(right.limits)]
    right.chromEnd.vec <- right.limits[-1]
    right.intervals <-
      paste0(right.chromStart.vec, "-", right.chromEnd.vec)
    rownames(right.cumsum.mat) <- c("before", right.intervals)
    left.limits <- bases.per.bin*(0:(n.cumsum-1))+left.chromStart
    left.chromStart.vec <- left.limits[-length(left.limits)]
    left.chromEnd.vec <- left.limits[-1]
    left.intervals <-
      paste0(left.chromStart.vec, "-", left.chromEnd.vec)
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
    possible.grid <- 
    expand.grid(left.cumsum.row=3:n.cumsum, right.cumsum.row=2:n.cumsum)
    possible.grid$left.chromStart <-
      left.chromStart.vec[possible.grid$left.cumsum.row-1]
    possible.grid$left.chromEnd <-
      left.chromEnd.vec[possible.grid$left.cumsum.row-1]
    possible.grid$right.chromStart <-
      right.chromStart.vec[possible.grid$right.cumsum.row-1]
    possible.grid$right.chromEnd <-
      right.chromEnd.vec[possible.grid$right.cumsum.row-1]
    feasible.grid <-
      subset(possible.grid,
             left.chromEnd <= right.chromStart &
             right.chromEnd < last.chromEnd)
    feasible.grid$model.i <- 1:nrow(feasible.grid)
    model.list <- list()
    seg.list <- list()
    for(model.i in feasible.grid$model.i){
      model.row <- feasible.grid[model.i, ]
      
      seg1.cumsums <- left.cumsum.mat[model.row$left.cumsum.row-1, ]
      seg1.chromEnd <- left.chromStart.vec[model.row$left.cumsum.row-1]
      seg1.bases <- seg1.chromEnd-max.chromStart
      seg1.means <- seg1.cumsums/seg1.bases
      seg1.loss <- OptimalPoissonLoss(seg1.cumsums, seg1.means)
      seg.list[[paste(model.i, 1)]] <-
        data.frame(chromStart=max.chromStart, chromEnd=seg1.chromEnd,
                   mean=seg1.means, sample.id=names(two.list),
                   model.i)
      
      right.cumsum.row <- right.cumsum.mat[model.row$right.cumsum.row, ]
      seg2.cumsums <- right.cumsum.row-seg1.cumsums
      seg2.chromEnd <- right.chromEnd.vec[model.row$right.cumsum.row-1]
      seg2.bases <- seg2.chromEnd-seg1.chromEnd
      seg2.means <- seg2.cumsums/seg2.bases
      seg2.loss <- OptimalPoissonLoss(seg2.cumsums, seg2.means)
      seg.list[[paste(model.i, 2)]] <-
        data.frame(chromStart=seg1.chromEnd, chromEnd=seg2.chromEnd,
                   mean=seg2.means, sample.id=names(two.list),
                   model.i)
      
      seg3.cumsums <- last.cumsum.vec-right.cumsum.row
      seg3.bases <- last.chromEnd-seg2.chromEnd
      seg3.means <- seg3.cumsums/seg3.bases
      seg3.loss <- OptimalPoissonLoss(seg3.cumsums, seg3.means)
      seg.list[[paste(model.i, 3)]] <-
        data.frame(chromStart=seg2.chromEnd, chromEnd=last.chromEnd,
                   mean=seg3.means, sample.id=names(two.list),
                   model.i)

      total.bases <- sum(seg1.bases + seg2.bases + seg3.bases)
      stopifnot(all.equal(total.bases, last.chromEnd - max.chromStart))

      total.loss <- sum(seg1.loss + seg2.loss + seg3.loss)
      model.list[[model.i]] <- data.frame(model.row, total.loss)
    }
    ## Plot the segment means as a reality check.
    seg.df <- do.call(rbind, seg.list)
    ggplot()+
    geom_step(aes(chromStart/1e3, count),
              data=data.frame(norm.df, what="data"),
              color="grey")+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 data=data.frame(seg.df, what="models"),
                 size=4, color="green")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ model.i, scales="free")

    ## Then plot the peaks only, colored by total cost of the model.
    model.df <- do.call(rbind, model.list)
    model.df$y <- -model.df$model.i * 0.1
    best.model <- model.df[which.min(model.df$total.loss), ]

    ggplot()+
    geom_step(aes(chromStart/1e3, count.norm),
              data=data.frame(norm.df, what="data"),
              color="grey")+
    geom_segment(aes(chromStart/1e3, mean.norm,
                     xend=chromEnd/1e3, yend=mean.norm),
                 data=data.frame(bin.df, what="bins"),
                 color="black")+
    geom_segment(aes(peakStart/1e3, 0,
                     xend=peakEnd/1e3, yend=0),
                 data=data.frame(loss.best, what="peak"),
                 color="green")+
    geom_segment(aes(left.chromStart/1e3, y,
                     color=total.loss,
                     xend=right.chromEnd/1e3, yend=y),
                 data=data.frame(model.df, what="models"),
                 size=4)+
    geom_text(aes(left.chromStart/1e3, y,
                  label="optimal "),
              data=best.model,
              hjust=1)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ .)

    peakStart <- best.model$left.chromStart
    peakEnd <- best.model$right.chromEnd
    left.cumsum.vec <- left.cumsum.mat[best.model$left.cumsum.row-2, ]
    right.cumsum.vec <- right.cumsum.mat[best.model$right.cumsum.row-1, ]
  }

  ## Find the best peak location across 4 samples.
  four <- subset(chr11ChIPseq$coverage,
                 118120000 < chromStart &
                 chromEnd < 118126000) 
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
