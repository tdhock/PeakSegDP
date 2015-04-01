library(data.table)
library(PeakSegDP)

exData <- system.file("exampleData", package="PeakSegDP")
argv <-
  file.path(exData, c("learned.model.RData", "bcell/McGill0091.bedGraph"))

argv <- commandArgs(trailingOnly=TRUE)

if(length(argv) != 2){
  stop("usage: Rscript Step3.R learned.model.RData profile.bedGraph")
}

model.RData.path <- normalizePath(argv[1], mustWork=TRUE)
bedGraph.path <- normalizeDir(argv[2])
base.dir <- sub("[.][a-zA-Z]*$", "", bedGraph.path)

cat("Reading learned model ", model.RData.path, "\n", sep="")
model.object.names <- load(model.RData.path)

cat("Reading profile coverage ", bedGraph.path, "\n", sep="")
sample.coverage <- fread(bedGraph.path)
setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
setkey(sample.coverage, chrom, chromStart, chromEnd)

## Set segmentation problem size (in bins/problem). For example 2000
## bins/problem means that 2000 adjacent bins will be used as the data
## to segment in 1 segmentation problem. NOTE: the cDPA takes
## quadratic O(bases.per.problem ^2) time, and is empirically faster
## than the log-linear time pDPA algorithm when bins.per.problem <
## 4000.
bins.per.problem <- 2000L
bases.per.problem <- bases.per.bin * bins.per.problem

chroms <- unique(sample.coverage$chrom)
problem.list <- list()
for(chrom in chroms){
  chrom.dir <- file.path(base.dir, bases.per.bin, chrom)
  chrom.coverage <- sample.coverage[chrom]
  chrom.bases <- max(chrom.coverage$chromEnd)
  first.chromStart <- chrom.coverage$chromStart[1]
  problems.df <-
    chromProblems(chrom, first.chromStart, chrom.bases, bases.per.problem)
  problem.list[[chrom]] <- data.table(problems.df)
}
problems.per.chrom <- sapply(problem.list, nrow)
total.problems <- sum(problems.per.chrom)
seconds.per.problem <- 1 # about 1 second per problem.
estimated.seconds <- total.problems * seconds.per.problem
seconds2hhmmss <- function(seconds){
  minutes <- seconds %/% 60
  seconds.extra <- seconds %% 60
  hours <- minutes %/% 60
  minutes.extra <- minutes %% 60
  sprintf("%02d:%02d:%02d", hours, minutes.extra, seconds.extra)
}
cat("Predicting peaks for ", total.problems,
    " segmentation problems, estimated time\n",
    seconds2hhmmss(estimated.seconds),
    " (", seconds.per.problem, " sec/problem)\n",
    sep="")

sample.peak.list <- list()
problem.time.list <- list()
for(chrom in chroms){
  chrom.dir <- file.path(base.dir, bases.per.bin, chrom)
  some.problems <- problem.list[[chrom]]
  chrom.coverage <- sample.coverage[chrom]
  setkey(some.problems, chromStart, chromEnd)
  setkey(chrom.coverage, chromStart, chromEnd)
  chrom.peak.list <- list()
  for(problem.i in 1:nrow(some.problems)){
    problem <- some.problems[problem.i, ]
    problem.coverage <- foverlaps(chrom.coverage, problem, nomatch = 0L)
    problem.name <- problem$problem.name
    problem.time.list[[problem.name]] <- system.time({
      RData.base <- paste0(problem.name, ".RData")
      RData.path <- file.path(chrom.dir, RData.base)
      RData.dir <- dirname(RData.path)
      dir.create(RData.dir, showWarnings = FALSE, recursive = TRUE)
      if(file.exists(RData.path)){
        load(RData.path)
      }else{
        cat(RData.path, "\n")

        seg.info <-
          segmentBins(chrom.coverage,
                      problem$chromStart,
                      bases.per.bin,
                      bins.per.problem)

        save(seg.info, file=RData.path)
      }#if(file.exists(RData.path))/else

      fmat <- rbind(seg.info$features)
      ## Set any non-finite features to 0 so we can predict a finite
      ## penalty value for all possible test data.
      fmat[!is.finite(fmat)] <- 0
      lmat <- learned.model$predict(fmat)
      pred.log.lambda <- as.numeric(lmat)

      selected <-
        subset(seg.info$modelSelection,
               min.log.lambda < pred.log.lambda &
                 pred.log.lambda < max.log.lambda)
      stopifnot(nrow(selected) == 1)
      param.name <- paste(selected$peaks)

      peaks <- data.table(seg.info$fit$peaks[[param.name]])
      if(nrow(peaks)){
        is.before <- peaks$chromEnd < problem$peakStart
        is.after <- problem$peakEnd < peaks$chromStart
        not.in.region <- is.before | is.after
        chrom.peak.list[[problem.name]] <- peaks[!not.in.region, ]
      }
    })[["elapsed"]]
  }#problem.i
  all.chrom.peaks <- do.call(rbind, chrom.peak.list)
  if(!is.null(all.chrom.peaks)){
    clustered.peaks <- clusterPeaks(all.chrom.peaks)
    table(clustered.peaks$cluster)
    clustered.peak.list <- split(clustered.peaks, clustered.peaks$cluster)
    reduced.peak.list <- list()
    for(cluster.str in names(clustered.peak.list)){
      peaks <- clustered.peak.list[[cluster.str]]
      reduced.peak.list[[cluster.str]] <- 
        data.frame(chromStart=peaks$chromStart[1],
                   chromEnd=peaks$chromEnd[nrow(peaks)])
    }
    reduced.peaks <- do.call(rbind, reduced.peak.list)
    sample.peak.list[[chrom]] <- data.table(chrom, reduced.peaks)
  }
}#chrom
elapsed.seconds <- sum(unlist(problem.time.list))
print(elapsed.seconds)
cat(seconds2hhmmss(as.integer(elapsed.seconds)),
    " for ",
    total.problems,
    " segmentation problems \n",
    sep="")

sample.peaks <- do.call(rbind, sample.peak.list)
cat("Peak counts per chromosome:\n")
table(sample.peaks$chrom)
out.path <- sub("[.]bedGraph$", "_peaks.bed", bedGraph.path)
cat("Writing ", nrow(sample.peaks),
    " peaks to ", out.path, "\n",
    sep="")

write.table(sample.peaks, file=out.path, quote=FALSE,
            sep="\t", row.names=FALSE, col.names=FALSE)
