library(data.table)
library(PeakSegDP)

exData <- system.file("exampleData", package="PeakSegDP")
argv <-
  file.path(exData, c("learned.model.RData", "bcell/McGill0322.bedGraph"))

argv <- commandArgs(trailingOnly=TRUE)

if(length(argv) != 2){
  stop("usage: Rscript Step3.R learned.model.RData profile.bedGraph")
}

RData.path <- normalizePath(argv[1], mustWork=TRUE)
bedGraph.path <- normalizePath(argv[2], mustWork=TRUE)

cat("Reading learned model ", RData.path, "\n", sep="")
model.object.names <- load(RData.path)

cat("Reading profile coverage ", bedGraph.path, "\n", sep="")
sample.coverage <- fread(bedGraph.path)
setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
setkey(sample.coverage, chrom, chromStart, chromEnd)

out.path <- sub("[.]bedGraph$", "_peaks.bed", bedGraph.path)

## Set segmentation problem size (in bins/problem). For example 2000
## bins/problem means that 2000 adjacent bins will be used as the data
## to segment in 1 segmentation problem. NOTE: the cDPA takes
## quadratic O(bases.per.problem ^2) time, and is empirically faster
## than the log-linear time pDPA algorithm when bins.per.problem <
## 4000.
bins.per.problem <- 2000L
bases.per.problem <- bases.per.bin * bins.per.problem

base.dir <- sub("[.][a-zA-Z]*$", "", bedGraph.path)
chroms <- unique(sample.coverage$chrom)
for(chrom in chroms){
  chrom.dir <- file.path(base.dir, bases.per.bin, chrom)
  chrom.coverage <- sample.coverage[chrom]
  chrom.bases <- max(chrom.coverage$chromEnd)
  problemEnd <-
    as.integer(seq(0, chrom.bases, by=bases.per.problem/2)[-(1:3)])
  problemStart <- as.integer(problemEnd-bases.per.problem)
  peakStart <- as.integer(problemStart + bases.per.problem/4)
  peakEnd <- as.integer(problemEnd - bases.per.problem/4)
  problems <- 
    data.table(chromStart=as.integer(c(0,
                 problemStart, problemEnd[length(problemEnd)-1])),
               peakStart=as.integer(c(0,
                 peakStart, peakEnd[length(peakEnd)])),
               peakEnd=as.integer(c(peakStart[1], peakEnd, chrom.bases)),
               chromEnd=as.integer(c(bases.per.problem,
                 problemEnd, chrom.bases)))
  problems[, problem.name :=
             sprintf("%s:%09d-%09d", chrom, chromStart, chromEnd)]
  first.chromStart <- chrom.coverage$chromStart[1]
  problem.before <- problems$chromEnd < first.chromStart
  some.problems <- problems[!problem.before, ]
  setkey(some.problems, chromStart, chromEnd)
  setkey(chrom.coverage, chromStart, chromEnd)
  for(problem.i in 1:nrow(some.problems)){
    problem <- some.problems[problem.i, ]
    problem.coverage <- foverlaps(chrom.coverage, problem, nomatch = 0L)
    problem.name <- problem$problem.name
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
    lmat <- learned.model$predict(fmat)
    pred.log.lambda <- as.numeric(lmat)

    selected <- subset(exact,
           min.log.lambda < pred.log.lambda &
             pred.log.lambda < max.log.lambda)
    stopifnot(nrow(selected) == 1)
    param.name <- paste(selected$peaks)

    pred.peaks <- fit$peaks[[param.name]]

    
  }#problem.i
  stop("TODO: combine peaks on this chrom")
}#chrom

stop("TODO: save all peaks to _peaks.bed file")
