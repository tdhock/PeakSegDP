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
library(animint)
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
argv <-
  c(paste0(base, ".bedGraph"),
    paste0(base, "_labels.bed"))

argv <- commandArgs(trailingOnly=TRUE)

print(argv)

bedGraph.path <- normalizePath(argv[1], mustWork=TRUE)
bed.path <- normalizePath(argv[2], mustWork=TRUE)

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

sample.quartile <- quantile(sample.coverage$count)

features.list <- list()
modelSelection.list <- list() # for computing test error.
peaks.list <- list() # for computing test error.
limits.list <- list()
error.list <- list()
all.problem.list <- list()# for deterimining where to keep peaks.
for(chrom in names(regions.by.chrom)){
  chrom.regions <- regions.by.chrom[[chrom]]
  chrom.coverage <- sample.coverage[chrom]
  chrom.quartile <- quantile(chrom.coverage$count)
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
        binSum.seconds <- system.time({
          bins <- 
            binSum(chrom.coverage,
                   bin.chromStart=problem.regions$chromStart[1],
                   bin.size=bases.per.bin)
        })[["elapsed"]]

        cDPA.seconds <- system.time({
          fit <- PeakSegDP(bins, maxPeaks = 9L)
        })[["elapsed"]]

        ## Compute feature vector for learning using this segmentation
        ## problem.
        bases <- with(bins, chromEnd-chromStart)
        long <- rep(bins$count, bases)
        n.bases <- sum(bases)
        n.data <- nrow(bins)

        uq <- quantile(bins$count)
        feature.vec <-
          c(unweighted.quartile=uq,
            weighted.quartile=quantile(long),
            sample.quantile.ratio=uq/sample.quartile,
            chrom.quantile.ratio=uq/chrom.quartile,
            sample.quantile=sample.quartile,
            chrom.quantile=chrom.quartile,
            unweighted.mean=mean(bins$count),
            weighted.mean=mean(long),
            bases=n.bases,
            data=n.data)
        suppressWarnings({
          features <-
            c(feature.vec,
              `log+1`=log(feature.vec+1),
              log=log(feature.vec),
              log.log=log(log(feature.vec)))
        })

        save(fit, features, bins,
             binSum.seconds, cDPA.seconds,
             file=RData.path)
      }#if(file.exists(RData.path))/else

      features.list[[paste(bases.per.bin)]][[problem.name]] <- features

      all.loss <- data.frame(fit$error)
      all.loss$cummin <- cummin(all.loss$error)
      loss <- subset(all.loss, error == cummin)
      rownames(loss) <- loss$segments
      bases <- with(fit$segments[1,], chromEnd-chromStart)
      in.sqrt <- 1.1 + log(bases / loss$segments)
      in.square <- 1 + 4 * sqrt(in.sqrt)
      complexity <- in.square * in.square * loss$segments
      
      exact <- with(loss, exactModelSelection(error, complexity, peaks))

      ## To map a predicted penalty value back to peaks, we need to
      ## store the exact model selection function and the list of
      ## peaks.
      modelSelection.list[[paste(bases.per.bin)]][[problem.name]] <- exact
      peaks.list[[paste(bases.per.bin)]][[problem.name]] <- fit$peaks
      
      exact$weighted.error <- NA
      rownames(exact) <- exact$peaks
      problem.regions[, chromStart := i.chromStart ]
      problem.regions[, chromEnd := i.chromEnd ]
      chunk.list <- split(problem.regions, problem.regions$chunk.id)

      show.error.list <- list()
      show.param.list <- list()
      show.exact.list <- list()
      for(chunk.id in names(chunk.list)){
        chunk.regions <- chunk.list[[chunk.id]]
        for(param.name in rownames(exact)){
          param.peaks <- fit$peaks[[param.name]]
          error <- PeakErrorChrom(param.peaks, chunk.regions)
          error$weight <- chunk.regions$weight
          error$weighted.error <- with(error, weight * (fp+fn))
          show.error.list[[paste(chunk.id, param.name)]] <-
            data.frame(chunk.id, peaks=param.name, error)
          exact[param.name, "weighted.error"] <- sum(error$weighted.error)
        }
        show.exact.list[[chunk.id]] <- data.frame(chunk.id, exact)
        ## For visualization, store the minimum error model with the
        ## least peaks.
        min.exact <- subset(exact, weighted.error == min(weighted.error))
        show.param.list[[chunk.id]] <- paste(min(min.exact$peaks))

        ## For learning, store the optimal interval of penalty values.
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
      show.exact <- do.call(rbind, show.exact.list)
      show.peaks <- do.call(rbind, fit$peaks)
      show.errors <- do.call(rbind, show.error.list)
      tit <- sub(dirname(dirname(bedGraph.path)), "", RData.path)
      viz <- 
        list(profile=ggplot()+
               ggtitle(tit)+
               xlab(paste("position on",
                          chrom,
                          "(kilo bases = kb)"))+
               guides(linetype=guide_legend(order=2,
                        override.aes=list(fill="white")))+
               scale_linetype_manual("error type", 
                                     values=c(correct=0,
                                       "false negative"=3,
                                       "false positive"=1))+
               geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                                 showSelected=peaks,
                                 linetype=status,
                                 fill=annotation),
                             alpha=0.5,
                             data=show.errors)+
               scale_fill_manual(values=ann.colors)+
               geom_line(aes(chromStart/1e3, count),
                         data=bins)+
               facet_grid(chunk.id ~ ., labeller=function(var, val){
                 paste("chunk", val)
               })+
               theme_animint(width=1500),

             title=tit,

             first=list(peaks=show.param.list[[1]]),

             selection=ggplot()+
             geom_segment(aes(min.log.lambda, weighted.error,
                              xend=max.log.lambda, yend=weighted.error),
                          data=show.exact)+
             scale_y_continuous(limits=c(0, NA))+
             facet_grid(chunk.id ~ ., labeller=function(var, val){
               paste("chunk", val)
             })+
             xlab("penalty log(lambda)")+
             geom_tallrect(aes(xmin=min.log.lambda, xmax=max.log.lambda,
                               clickSelects=peaks),
                           alpha=0.5,
                           data=show.exact))
      png.path <- sub("RData", "png", RData.path)
      png(png.path, units="in", res=200, width=15, height=4)
      print(viz$profile)
      dev.off()
      thumb.path <- sub("[.]png", "-thumb.png", png.path)
      cmd <- sprintf("convert %s -resize 230 %s", png.path, thumb.path)
      system(cmd)
      if(is.data.frame(show.peaks) && nrow(show.peaks)){
        viz$profile <- viz$profile+
          geom_segment(aes(chromStart/1e3, 0,
                           showSelected=peaks,
                           xend=chromEnd/1e3, yend=0),
                       data=show.peaks,
                       size=4,
                       color="deepskyblue")+
          geom_point(aes(chromStart/1e3, 0,
                         showSelected=peaks),
                     data=show.peaks,
                     size=5,
                     color="deepskyblue")
      }
      viz.path <- sub("[.]RData$", "", RData.path)
      animint2dir(viz, out.dir=viz.path)
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
error.png <- sub("[.]bedGraph$", "_weightedError.png", bedGraph.path)
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
         problems=all.problem.list[[bases.per.bin.str]],
         modelSelection=modelSelection.list[[bases.per.bin.str]],#list[problem]
         peaks=peaks.list[[bases.per.bin.str]])#list[problem]
}

out.RData <- sub("[.]bedGraph$", "_residuals.RData", bedGraph.path)
save(features.limits, # used for the learning/training.
     ## (limits separated by chunk for cross-validation).
     errors, # used for selecting the best resolution before training.
     regions,
     ## For this model we will need to compute per-chromosome errors, for
     ## several different test sets of chunk ids. So save regions so we can
     ## compute test error later.
     file=out.RData)
