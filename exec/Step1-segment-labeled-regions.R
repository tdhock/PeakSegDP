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

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

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

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

bedGraph.path <- argv[1]
bed.path <- argv[2]

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
limits.list <- list()
error.list <- list()
for(chrom in names(regions.by.chrom)){
  chrom.regions <- regions.by.chrom[[chrom]]
  chrom.coverage <- sample.coverage[chrom]
  chrom.bases <- max(chrom.coverage$chromEnd)
  for(problem.size.i in 1:nrow(feasible.resolutions)){
    res.row <- feasible.resolutions[problem.size.i, ]
    bases.per.problem <- as.integer(res.row$bases.per.problem)
    bases.per.bin <- as.integer(res.row$bases.per.bin)
    chrom.dir <- file.path(base.dir, bases.per.bin, chrom)

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
      if(file.exists(RData.path)){
        load(RData.path)
      }else{
        cat(RData.path, "\n")
        binSum.seconds <- system.time({
          bins <- 
            binSum(chrom.coverage,
                   bin.chromStart=problem.regions$chromStart[1],
                   bin.size=bases.per.bin)
        })[["elapsed"]]
        
        ggplot()+
          geom_tallrect(aes(xmin=i.chromStart/1e3, xmax=i.chromEnd/1e3,
                            fill=annotation),
                        color="grey",
                        data=problem.regions)+
          scale_fill_manual(values=ann.colors)+
          geom_step(aes(chromStart/1e3, count),
                    data=bins)+
          theme_bw()+
          theme(panel.margin=grid::unit(0, "cm"))+
          facet_grid(chunk.id ~ .)
        
        cDPA.seconds <- system.time({
          fit <- PeakSegDP(bins, maxPeaks = 9L)
        })[["elapsed"]]

        ## Compute feature vector for learning using this segmentation
        ## problem.
        bases <- with(bins, chromEnd-chromStart)
        long <- rep(bins$count, bases)
        n.bases <- sum(bases)
        n.data <- nrow(bins)

        feature.vec <-
          c(unweighted.quartile=quantile(bins$count),
            weighted.quartile=quantile(long),
            unweighted.mean=mean(bins$count),
            weighted.mean=mean(long),
            bases=n.bases,
            data=n.data)        
        features <-
          c(feature.vec,
            `log+1`=log(feature.vec+1),
            log=log(feature.vec),
            log.log=log(log(feature.vec)))

        RData.dir <- dirname(RData.path)
        dir.create(RData.dir, showWarnings = FALSE, recursive = TRUE)
        save(fit, features,
             binSum.seconds, cDPA.seconds,
             file=RData.path)
      }#if(file.exists(RData.path))/else

      features.list[[paste(bases.per.bin)]][[problem.name]] <- features
      
      all.loss <- data.frame(fit$error)
      all.loss$cummin <- cummin(all.loss$error)
      loss <- subset(all.loss, error == cummin)
      rownames(loss) <- loss$segments
      exact <- with(loss, exactModelSelection(error, segments))
      exact$segments <- exact$model.complexity
      exact$weighted.error <- NA
      rownames(exact) <- exact$segments
      exact$peaks <- loss[paste(exact$segments), "peaks"]
      rownames(exact) <- exact$peaks
      problem.regions[, chromStart := i.chromStart ]
      problem.regions[, chromEnd := i.chromEnd ]
      chunk.list <- split(problem.regions, problem.regions$chunk.id)
      for(chunk.id in names(chunk.list)){
        chunk.regions <- chunk.list[[chunk.id]]
        for(param.name in rownames(exact)){
          param.peaks <- fit$peaks[[param.name]]
          error <- PeakErrorChrom(param.peaks, chunk.regions)
          error$weight <- chunk.regions$weight
          error$weighted.error <- with(error, weight * (fp+fn))
          exact[param.name, "weighted.error"] <- sum(error$weighted.error)
        }
        limits <- with(exact, {
          largestContinuousMinimum(weighted.error,
                                   max.log.lambda-min.log.lambda)
        })

        limits.list[[paste(bases.per.bin)]][[chunk.id]][[problem.name]] <- 
          c(exact$min.log.lambda[limits$start],
            exact$max.log.lambda[limits$end])

        result.dt <-
          data.table(chrom,
                     problem.name,
                     chunk.id,
                     bases.per.bin,
                     bases.per.problem,
                     binSum.seconds,
                     cDPA.seconds,
                     min.peaks=exact$peaks[limits$end],
                     max.peaks=exact$peaks[limits$start],
                     weighted.error=min(exact$weighted.error),
                     total.weight=sum(chunk.regions$weight),
                     regions=nrow(chunk.regions))

        error.list[[paste(bases.per.bin, problem.name, chunk.id)]] <-
          check.list[[paste(bases.per.bin, problem.name, chunk.id)]] <-
            result.dt
      }#chunk.id
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

ggplot(res.errors, aes(bases.per.bin, weighted.error))+
  geom_line()+
  geom_point()+
  scale_x_log10()

features.limits <- list()
for(bases.per.bin.str in names(features.list)){
  chunk.list <- limits.list[[bases.per.bin.str]]
  chunk.mats <- list()
  for(chunk.id in names(chunk.list)){
    chunk.mats[[chunk.id]] <- do.call(rbind, chunk.list[[chunk.id]])
  }
  features.limits[[bases.per.bin.str]] <-
    list(features=do.call(rbind, features.list[[bases.per.bin.str]]),
         limits=chunk.mats)
}

out.RData <- sub("[^.]*$", "RData", bed.path)
save(features.limits, # used for the learning/training.
     ## (limits separated by chunk for cross-validation).
     errors, # used for selecting the best resolution before training.
     regions,
     ## For this model we will need to compute per-chromosome errors, for
     ## several different test sets of chunk ids. So save regions so we can
     ## compute test error later.
     file=out.RData)
