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
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(two, 2)
  })[["elapsed"]]
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(two)
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
  heuristic.seconds <- system.time({
    heuristic <- multiSampleSegHeuristic(four, 2)
  })[["elapsed"]]
  optimal.seconds <- system.time({
    optimal <- multiSampleSegOptimal(four)
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
