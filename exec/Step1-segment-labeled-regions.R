## Set grid of problem resolutions (in bases/bin). For example 20
## bases/bin means that the sum of aligned reads over 20 bases will be
## used to create 1 bin, which will be used as one data point in each
## segmentation problem.
resolution.grid <- 
  expand.grid(x=c(1, 2, 4, 8),
              y=c(1, 10, 100, 1000))
bases.per.bin.grid <- as.integer(with(resolution.grid, x * y))

## Set segmentation problem size (in bins/problem). For example 2000
## bins/problem means that 2000 adjacent bins will be used as the data
## to segment in 1 segmentation problem. NOTE: the cDPA takes
## quadratic O(bases.per.problem ^2) time, and is empirically faster
## than the log-linear time pDPA algorithm when bins.per.problem <
## 4000.
bins.per.problem <- 2000L

library(PeakError)
library(PeakSegDP)
library(ggplot2)
library(data.table) #"1.9.4",

argv <-
  c("/home/thocking/genomecov/H3K36me3/McGill0322.bedGraph",
                    # -> output files /McGill0322/2/chr1/chr1:0000-4000.RData
    "/home/thocking/genomelabels/H3K36me3_TDH_immune/McGill0322.bed")
                                  # -> output file  /McGill0322.RData

argv <-
  c("/home/thocking/genomecov/H3K4me3/McGill0001.bedGraph",
    "/home/thocking/genomelabels/H3K4me3_PGP_immune/McGill0001.bed")

base <-
  file.path(system.file("exampleData", package="PeakSegDP"),
            "bcell", "McGill0091")
base <-
  file.path(system.file("exampleData", package="PeakSegDP"),
            "bcell", "McGill0322")
argv <-
  c(paste0(base, ".bedGraph"),
    paste0(base, "_labels.bed"))

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

if(length(argv) != 2){
  stop("usage: Rscript Step1.R profile.bedGraph labels.bed")
}

bedGraph.path <- normalizeDir(argv[1])
bed.path <- normalizeDir(argv[2])

regions.wide <- fread(bed.path)
region.names <- c("chrom", "chromStart", "chromEnd", "annotation", "chunk.id")
regions <- regions.wide[, seq_along(region.names), with=FALSE]
setnames(regions, region.names)
regions[, bases := chromEnd - chromStart]

peak.ranges <- 
  regions[annotation %in% c("peakStart", "peakEnd"),
          .(min.bases=min(bases), max.bases=max(bases)) ]
candidate.resolutions <- 
  data.table(peak.ranges,
             bases.per.problem=bases.per.bin.grid * bins.per.problem,
             bases.per.bin=bases.per.bin.grid)
feasible.resolutions <-
  candidate.resolutions[max.bases < bases.per.problem &
                        bases.per.bin < min.bases, ]
cat("Using min/max size of peakStart/peakEnd regions,\n",
    "feasible resolutions are:\n", sep="")
print(feasible.resolutions)

regions.by.chrom <- split(regions, regions$chrom)

sample.coverage <- fread(bedGraph.path)
setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
setkey(sample.coverage, chrom, chromStart, chromEnd)

base.dir <- sub("[.][a-zA-Z]*$", "", bedGraph.path)

features.list <- list()
bin.list <- list()
modelSelection.list <- list() # for computing test error.
peaks.list <- list() # for computing test error.
weighted.error.mat.list <- list()
limits.list <- list()
error.list <- list()
error.region.list <- list()
all.problem.list <- list()# for deterimining where to keep peaks.
for(chrom in names(regions.by.chrom)){
  chrom.regions <- regions.by.chrom[[chrom]]
  chrom.coverage <- sample.coverage[chrom]
  first.chromStart <- chrom.coverage$chromStart[1]
  chrom.bases <- max(chrom.coverage$chromEnd)
  for(problem.size.i in 1:nrow(feasible.resolutions)){
    res.row <- feasible.resolutions[problem.size.i, ]
    bases.per.problem <- as.integer(res.row$bases.per.problem)
    bases.per.bin <- as.integer(res.row$bases.per.bin)
    chrom.dir <- file.path(base.dir, bases.per.bin, chrom)

    problems <-
      chromProblems(chrom, first.chromStart, chrom.bases, bases.per.problem)

    all.problem.list[[paste(bases.per.bin)]][[chrom]] <- problems

    ## Figure out which problems overlap which regions.
    setkey(problems, chromStart, chromEnd)
    setkey(chrom.regions, chromStart, chromEnd)
    chrom.regions$region.i <- 1:nrow(chrom.regions)
    ov.regions <- foverlaps(chrom.regions, problems)
    stopifnot(chrom.regions$region.i %in% ov.regions$region.i)
    region.weight <- 1/table(ov.regions$region.i)
    ov.regions$weight <- region.weight[paste(ov.regions$region.i)]
    stopifnot(sum(ov.regions$weight) == nrow(chrom.regions))
    problem.list <- split(ov.regions, ov.regions$problem.name)
    data.list <- list()
    check.list <- list()
    for(problem.name in names(problem.list)){
      RData.base <- paste0(problem.name, ".RData")
      RData.path <- file.path(chrom.dir, RData.base)
      problem.regions <- problem.list[[problem.name]]
      RData.dir <- dirname(RData.path)
      dir.create(RData.dir, showWarnings = FALSE, recursive = TRUE)
      if(file.exists(RData.path)){
        load(RData.path)
      }else{
        cat(RData.path, "\n")

        seg.info <-
          segmentBins(chrom.coverage,
                      problem.regions$chromStart[1],                      
                      bases.per.bin,
                      bins.per.problem)
        
        save(seg.info, file=RData.path)
      }#if(file.exists(RData.path))/else

      features.list[[paste(bases.per.bin)]][[problem.name]] <- seg.info$features
      bin.list[[paste(bases.per.bin)]][[problem.name]] <- seg.info$bins

      ## To map a predicted penalty value back to peaks, we need to
      ## store the exact model selection function and the list of
      ## peaks.
      exact <- seg.info$modelSelection
      modelSelection.list[[paste(bases.per.bin)]][[problem.name]] <- exact
      peaks.list[[paste(bases.per.bin)]][[problem.name]] <-
        seg.info$fit$peaks
      
      exact$weighted.error <- NA
      rownames(exact) <- exact$peaks
      problem.regions[, chromStart := i.chromStart ]
      problem.regions[, chromEnd := i.chromEnd ]
      chunk.list <- split(problem.regions, problem.regions$chunk.id)
      weighted.error.mat <-
        matrix(NA, nrow(exact), length(chunk.list),
               dimnames=list(peaks=exact$peaks, chunk.id=names(chunk.list)))
      problem.error.list <- list()
      for(chunk.id in names(chunk.list)){
        chunk.regions <- chunk.list[[chunk.id]]
        for(param.name in rownames(exact)){
          param.peaks <- seg.info$fit$peaks[[param.name]]
          error <- PeakErrorChrom(param.peaks, chunk.regions)
          error$weight <- chunk.regions$weight
          error$weighted.error <- with(error, weight * (fp+fn))
          problem.error.list[[chunk.id]][[param.name]] <-
            data.frame(chunk.id, peaks=param.name, error)
          exact[param.name, "weighted.error"] <- sum(error$weighted.error)
        }
        weighted.error.mat[, chunk.id] <- exact$weighted.error
        ## For visualization, store the minimum error model with the
        ## least peaks.
        min.exact <- subset(exact, weighted.error == min(weighted.error))

        ## For learning, store the optimal interval of penalty values.
        limits <- with(exact, {
          largestContinuousMinimum(weighted.error,
                                   max.log.lambda-min.log.lambda)
        })

        weighted.error.mat.list[[paste(bases.per.bin)]][[problem.name]] <-
          weighted.error.mat

        limits.list[[paste(bases.per.bin)]][[chunk.id]][[problem.name]] <- 
          c(exact$min.log.lambda[limits$start],
            exact$max.log.lambda[limits$end])

        result.dt <-
          data.table(chrom,
                     problem.name,
                     chunk.id,
                     bases.per.bin,
                     bases.per.problem,
                     binSum.seconds=seg.info$binSum.seconds,
                     cDPA.seconds=seg.info$cDPA.seconds,
                     min.peaks=exact$peaks[limits$end],
                     max.peaks=exact$peaks[limits$start],
                     weighted.error=min(exact$weighted.error),
                     total.weight=sum(chunk.regions$weight),
                     regions=nrow(chunk.regions))

        error.list[[paste(bases.per.bin, problem.name, chunk.id)]] <-
          check.list[[paste(bases.per.bin, problem.name, chunk.id)]] <-
            result.dt
      }#chunk.id
      error.region.list[[paste(bases.per.bin)]][[problem.name]] <-
        problem.error.list
    }#problem.name
    check.dt <- do.call(rbind, check.list)
    stopifnot(all.equal(sum(check.dt$total.weight),
                        nrow(chrom.regions)))
  }#problem.size.i
}#chrom

errors <- do.call(rbind, error.list)

res.errors <-
  errors[,
         .(weighted.error=sum(weighted.error),
           total.weight=sum(total.weight)),
         by=bases.per.bin]
exp.weight <- rep(res.errors$total.weight[1], nrow(res.errors))
stopifnot(all.equal(res.errors$total.weight, exp.weight))

sampleError <- 
ggplot(res.errors, aes(bases.per.bin, weighted.error))+
  geom_line()+
  geom_point()+
  scale_x_log10()
## The error function depends on the labels, so the error plot should
## be saved with the labels/bed file (which may or may not be in the
## same directory as the bedGraph file).
error.png <- sub("labels[.]bed$", "weightedError.png", bed.path)
png(error.png, units="in", res=200, width=12, height=7)
print(sampleError)
dev.off()

features.limits <- list()
for(bases.per.bin.str in names(features.list)){
  chunk.list <- limits.list[[bases.per.bin.str]]
  chunk.mats <- list()
  for(chunk.id in names(chunk.list)){
    chunk.mats[[chunk.id]] <- do.call(rbind, chunk.list[[chunk.id]])
  }
  features.limits[[bases.per.bin.str]] <-
    list(features=do.call(rbind, features.list[[bases.per.bin.str]]),
         ## matrix[problem, feature]
         limits=chunk.mats, #matrix[problem, ]
         bins=bin.list[[bases.per.bin.str]],
         weighted.error.mats=weighted.error.mat.list[[bases.per.bin.str]],
         problems=all.problem.list[[bases.per.bin.str]],
         error.regions=error.region.list[[bases.per.bin.str]],
         modelSelection=modelSelection.list[[bases.per.bin.str]],#list[problem]
         peaks=peaks.list[[bases.per.bin.str]])#list[problem]
}

## Finally, we will want to display model predictions for each chunk
## in Step4, so save some data in those regions for display later.
region.list <- split(regions, regions$chunk.id)
chunk.info <- list()
chunk.dt.list <- list()
for(chunk.id in names(region.list)){
  chunk.regions <- region.list[[chunk.id]]
  chrom <- paste(chunk.regions$chrom[1])
  chrom.coverage <- sample.coverage[chrom]
  chunkStart <- min(chunk.regions$chromStart)
  chunkEnd <- max(chunk.regions$chromEnd)
  chunk.bases <- chunkEnd-chunkStart
  n.bins <- 2000L
  bases.per.bin <- as.integer(chunk.bases/n.bins)
  bins <- binSum(chrom.coverage, chunkStart, bases.per.bin, n.bins)
  chunk <- data.table(chunk.id, chrom, chunkStart, chunkEnd)
  chunk.dt.list[[chunk.id]] <- chunk
  chunk.info[[chunk.id]] <-
    list(regions=chunk.regions,
         bins=data.table(bins),
         chunk=chunk)
}
chunk.dt <- do.call(rbind, chunk.dt.list)

out.RData <- sub("labels[.]bed$", "residuals.RData", bed.path)
save(features.limits, # used for the learning/training.
     ## (limits separated by chunk for cross-validation).
     errors, # used for selecting the best resolution before training.
     regions,
     chunk.dt,
     chunk.info, # used for model visualization.
     ## For this model we will need to compute per-chromosome errors, for
     ## several different test sets of chunk ids. So save regions so we can
     ## compute test error later.
     file=out.RData)
