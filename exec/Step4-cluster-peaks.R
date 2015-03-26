library(data.table)
library(xtable)
library(PeakSegDP)
library(PeakError)

exData <- "~/exampleData"
exData <- system.file("exampleData", package="PeakSegDP")
argv <- file.path(exData, "learned.model.RData")

argv <- commandArgs(trailingOnly=TRUE)

if(length(argv) != 1){
  stop("usage: Rscript Step4.R path/to/dir
where there are dir/*/*_peaks.bed and dir/*/*_residuals.RData files")
}

data.dir <- normalizePath(argv[1], mustWork=TRUE)

bed.files <- Sys.glob(file.path(data.dir, "*", "*_peaks.bed"))

chunk.list <- list()
for(bed.file in bed.files){
  peaks <- fread(bed.file)
  setnames(peaks, c("chrom", "chromStart", "chromEnd"))
  sample.id <- sub("_peaks.bed$", "", basename(bed.file))
  residuals.file <- sub("peaks.bed$", "residuals.RData", bed.file)
  robjs <- load(residuals.file)
  setkey(peaks, chrom, chromStart, chromEnd)
  setkey(chunk.dt, chrom, chunkStart, chunkEnd)
  joined <- foverlaps(peaks, chunk.dt, nomatch = 0L)
  setkey(joined, chunk.id, chromStart, chromEnd)
  for(chunk.id in names(chunk.info)){
    one.chunk <- chunk.info[[chunk.id]]
    one.chunk$peaks <- joined[chunk.id]
    one.chunk$error <- with(one.chunk, PeakError(peaks, regions))
    for(data.type in c("bins", "peaks", "error")){
      one.dt <- one.chunk[[data.type]]
      chunk.list[[chunk.id]][[data.type]][[sample.id]] <-
        data.table(sample.id, one.dt)
    }
  }
  peak.list[[sample.id]] <- data.table(sample.id, peaks)
  region.list[[sample.id]] <- data.table(sample.id, regions)
}

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

figure.dir <- "predicted-figures"
figure.path <- file.path(data.dir, figure.dir)
dir.create(figure.path, showWarnings = FALSE)

table.list <- list()
for(chunk.id in names(chunk.list)){
  one <- 
    lapply(chunk.list[[chunk.id]], function(type.list){
      do.call(rbind, type.list)
    })
  png.base <- sprintf("chunk%s.png", chunk.id)
  png.path <- file.path(figure.path, png.base)
  thumb.base <- sprintf("chunk%s-thumb.png", chunk.id)
  thumb.path <- file.path(figure.path, thumb.base)
  chunk.totals <- 
  one$error[, .(
    errors=sum(fp+fn),
    regions=length(fp),
    fp=sum(fp),
    possible.fp=sum(possible.fp),
    fn=sum(fn),
    possible.fn=sum(possible.tp)
    )]
  chunk.totals$href <-
    sprintf('<a href="%s"><img src="%s" /></a>',
            png.base, thumb.base)
  
  peaks <- clusterPeaks(one$peaks)
  tit <- with(chunk.totals, {
    sprintf("%d/%d train errors (%d/%d false positives, %d/%d false negatives)",
            errors, regions, fp, possible.fp, fn, possible.fn)
  })
  chunkPlot <- 
  ggplot()+
    ggtitle(tit)+
    xlab(paste("position on", one$error$chrom[1], "(kilo bases = kb)"))+
    scale_fill_manual(values=ann.colors)+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  color="grey",
                  alpha=0.5,
                  data=one$error)+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    ylab("aligned read coverage")+
    facet_grid(sample.id ~ ., scales="free")+
    coord_cartesian(xlim=with(one$bins, c(chromStart[1], max(chromEnd))/1e3))+
    geom_step(aes(chromStart/1e3, mean), data=one$bins, color="grey50")+
    scale_linetype_manual("error type",
                          values=c(correct=0, "false positive"=1,
                            "false negative"=3))+
    guides(linetype=guide_legend(order=2,
             override.aes=list(fill="white")))+
    scale_color_discrete("cluster")+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      linetype=status),
                  fill=NA,
                  color="black",
                  size=1,
                  data=one$error)+
    geom_segment(aes(chromStart/1e3, 0,
                     color=factor(cluster),
                     xend=chromEnd/1e3, yend=0),
                 ##color="deepskyblue",
                 size=2,
                 data=peaks)

  png(png.path, width=15, h=10, units="in", res=200)
  print(chunkPlot)
  dev.off()

  cmd <- sprintf("convert %s -resize 230 %s", png.path, thumb.path)
  system(cmd)

  table.list[[chunk.id]] <- chunk.totals
}

xdt <- do.call(rbind, table.list)
xt <- xtable(xdt)
html.file <- file.path(figure.path, "index.html")
print(xt, type="html", file=html.file,
      include.rownames=FALSE, sanitize.text.function=identity)
