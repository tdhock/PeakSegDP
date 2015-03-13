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
