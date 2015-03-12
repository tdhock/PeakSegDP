multiSampleSeg <- structure(function
### Find one peak common to several samples.
(profiles,
 n.bins=100L
 ){
  if(is.data.frame(profiles)){
    profiles <- split(profiles, profiles$sample.id)
  }
  stopifnot(is.list(profiles))
  for(profile.i in seq_along(profiles)){
    profiles[[profile.i]] <-
      data.frame(profiles[[profile.i]])[, c("chromStart", "chromEnd", "count")]
  }
  .Call("multiSampleSeg_interface",
        profiles,
        n.bins,
        PACKAGE="PeakSegDP")
}, ex=function(){
  library(PeakSegDP)
  data(chr11ChIPseq)
  two <- subset(chr11ChIPseq$coverage,
                118090000 < chromStart &
                chromEnd < 118100000 &
                sample.id %in% c("McGill0002", "McGill0004"))
  library(ggplot2)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=two)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")
  four <- subset(chr11ChIPseq$coverage,
                 118120000 < chromStart &
                 chromEnd < 118126000)
  ggplot()+
    geom_step(aes(chromStart/1e3, count), data=four)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    facet_grid(sample.id ~ ., scales="free")
})
