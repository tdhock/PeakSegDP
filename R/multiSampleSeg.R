multiSampleSegHeuristic <- structure(function
### Find one peak common to several samples.
(profiles,
 n.bins=100L
 ){
  stopifnot(is.numeric(n.bins))
  stopifnot(length(n.bins)==1)
  stopifnot(n.bins >= 3)
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
          as.integer(n.bins),
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
  heuristic <- multiSampleSegHeuristic(two, 150)
  optimal <- multiSampleSegOptimal(two)
  peaks <-
    rbind(data.frame(optimal, model="optimal"),
          data.frame(heuristic, model="heuristic"))
  library(ggplot2)
  ggplot()+
    scale_size_manual(values=c(optimal=2, heuristic=1))+
    geom_step(aes(chromStart/1e3, count), data=two)+
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
  peak <- multiSampleSeg(four)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=four)+
    geom_segment(aes(chromStart/1e3, 0,
                     xend=chromEnd/1e3, yend=0),
                 color="green",
                 data=peak)+
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
