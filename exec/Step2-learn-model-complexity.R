library(animint)
library(data.table)
library(xtable)
library(ggplot2)
library(PeakSegDP)

argv <- c("/home/thocking/genomelabels/H3K36me3_TDH_immune/")
argv <- system.file("exampleData", package="PeakSegDP")

argv <- commandArgs(trailingOnly=TRUE)

no.trailing <- normalizePath(argv, mustWork=TRUE)
RData.files <- Sys.glob(file.path(no.trailing, "*", "*.RData"))
cat("Input files:\n")
print(RData.files)
stopifnot(length(RData.files) > 0)
names(RData.files) <- sub("_residuals[.]RData$", "", basename(RData.files))
out.RData <- file.path(no.trailing, "learned.model.RData")
cat("Output file:\n", out.RData, "\n", sep="")

error.list <- list()
features.limits.list <- list()
for(sample.id in names(RData.files)){
  RData.file <- RData.files[[sample.id]]
  objs <- load(RData.file)
  error.list[[sample.id]] <- data.table(sample.id, errors)
  features.limits.list[[sample.id]] <- features.limits
}

errors <- do.call(rbind, error.list)
res.errors <- 
  errors[,
         .(weighted.error=sum(weighted.error),
           total.weight=sum(total.weight)),
         by=bases.per.bin]

## consider only resolutions with max weights.
max.weight.value <- max(res.errors$total.weight)
max.weight <- res.errors[total.weight == max.weight.value, ]

## There could be several minimum error resolutions.
min.err.value <- min(max.weight$weighted.error)
min.err <- max.weight[weighted.error == min.err.value, ]

## The largest resolution will take the least amount of CPU time.
bases.per.bin.dt <- min.err[which.max(bases.per.bin), ]
bases.per.bin <- bases.per.bin.dt$bases.per.bin
bases.per.bin.str <- paste(bases.per.bin)

y.range <- max(max.weight$weighted.error)-min.err.value
errPlot <- 
ggplot()+
  ggtitle("selecting the best resolution")+
  geom_line(aes(bases.per.bin, weighted.error), data=max.weight)+
  geom_point(aes(bases.per.bin, weighted.error), data=bases.per.bin.dt)+
  geom_text(aes(bases.per.bin, weighted.error-y.range/20,
                label=paste(bases.per.bin, "bases/bin")),
            data=bases.per.bin.dt)+
  scale_x_log10()+
  theme_grey()

## idea is to make train-figures directory with index.html, several
## png files and several animint data viz, which is easily copyable
## for sharing results.
figures.dir <- "train-figures"
error.base <- "figure-weightedError.png"
figures.path <- file.path(no.trailing, figures.dir)
dir.create(figures.path, showWarnings = FALSE)
error.png <- file.path(figures.path, error.base)
png(error.png, units="in", res=100, width=5, height=3)
print(errPlot)
dev.off()

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")

## Now that we know what is the best resolution, make some plots that
## show the training model residuals at that resolution.
for(sample.id in names(features.limits.list)){
  sample.data <- features.limits.list[[sample.id]][[bases.per.bin.str]]
  exact.list <- sample.data$modelSelection
  for(problem.name in names(exact.list)){
    exact <- exact.list[[problem.name]]
    weighted.error.mat <- sample.data$weighted.error.mats[[problem.name]]
    show.error.list <- list()
    show.param.list <- list()
    show.exact.list <- list()
    for(chunk.id in colnames(weighted.error.mat)){
      exact$weighted.error <- weighted.error.mat[, chunk.id]
      show.exact.list[[chunk.id]] <- data.frame(chunk.id, exact)
      best <- subset(exact, weighted.error == min(weighted.error))
      show.param.list[[chunk.id]] <- paste(min(best$peaks))
      error.region.list <- sample.data$error.regions[[problem.name]][[chunk.id]]
      show.error.list[[chunk.id]] <- do.call(rbind, error.region.list)
    }#chunk.id
    show.exact <- do.call(rbind, show.exact.list)
    show.peaks <- do.call(rbind, sample.data$peaks[[problem.name]])
    show.errors <- do.call(rbind, show.error.list)
    bins <- sample.data$bins[[problem.name]]
    problem.path <- file.path(figures.path, sample.id, problem.name)
    dir.create(problem.path, showWarnings = FALSE, recursive = TRUE)
    tit <- paste(sample.id, problem.name)
    chrom <- sub(":.*", "", problem.name)
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

    png.path <- paste0(problem.path, ".png")
    png(png.path, units="in", res=100, width=8, height=4)
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
    animint2dir(viz, out.dir=problem.path, open.browser=FALSE)
  }#problem.name
}#sample.id

errors.by.res <- split(errors, errors$bases.per.bin)
errors.best.res <- errors.by.res[[bases.per.bin.str]]
errors.ordered <- errors.best.res[order(weighted.error, decreasing = TRUE), ]
errors.ordered[, sample.dir := sample.id]
errors.ordered[, problem.dir := file.path(sample.dir, problem.name)]
errors.ordered[, thumb.png := paste0(problem.dir, "-thumb.png")]
errors.ordered[, index := file.path(problem.dir, "index.html")]
errors.ordered[, thumb.href := sprintf('<a href="%s"><img src="%s" /></a>',
                              index, thumb.png)]
errors.ordered[, chrom.fac := factor(chrom, paste0("chr", c(1:22, "X", "Y")))]
errors.ordered[, chromRange := sub(".*:", "", problem.name)]
errors.ordered[, chromEnd := as.integer(sub(".*-", "", chromRange))]
errors.ordered[, chromStart := as.integer(sub("-.*", "", chromRange))]
errors.ordered[, errors := ifelse(weighted.error == 0, "none", "some")]

## TODO: add chrom diagram viz with 
## ggplot()+
##   geom_point(aes(chromEnd, chrom.fac,
##   href=index, color=errors), data=errors.ordered)+
##   scale_color_manual(values=c(none="grey", some="black"))+
##   facet_grid(sample.id ~ .)
rcopy <- data.table(regions)
rcopy[, chrom.fac := factor(paste(chrom), paste0("chr", c(1:22, "X", "Y")))]
tit <-
  paste0("training data:\n",
         nrow(regions), " labels per sample (color),\n",
         length(RData.files), " samples (", sample.id, " shown),\n",
         length(unique(regions$chunk.id)), " chunks,\n",
         length(unique(errors.ordered$problem.name)),
         " segmentation problems (black)")
chromPlot <- 
ggplot()+
  ggtitle(tit)+
  theme_grey()+
  ylab("chromosome")+
  xlab("position on chromosome (kilo bases = kb)")+
  geom_segment(aes(chromStart/1e3, chrom.fac,
                   xend=chromEnd/1e3, yend=chrom.fac),
               size=2,
               data=errors.ordered)+
  geom_point(aes(chromEnd/1e3, chrom.fac), data=errors.ordered,
             size=3)+
  scale_color_manual(values=ann.colors)+
  geom_segment(aes(chromStart/1e3, chrom.fac,
                   color=annotation,
                   xend=chromEnd/1e3, yend=chrom.fac),
               size=1,
               data=rcopy)+
  geom_point(aes(chromEnd/1e3, chrom.fac, color=annotation),
             size=2,
             data=rcopy)
chrom.base <- "figure-training-data-chroms.png"
chrom.path <- file.path(figures.path, chrom.base)
png(chrom.path, units="in", res=100, width=10, height=8)
print(chromPlot)
dev.off()


xt.cols <-
  c("weighted.error", "total.weight",
    "thumb.href", "sample.id", "problem.name",
    "min.peaks", "max.peaks")
errors.df <- data.frame(errors.ordered)[, xt.cols]
xt <- xtable(errors.df)

report.file <- file.path(figures.path, "index.html")
cat(sprintf('<img src="%s" /> <br />', chrom.base),
    sprintf('<img src="%s" /> <br />', error.base),
    file=report.file)
print(xt, type="html", file=report.file, append=TRUE,
      sanitize.text.function=identity)
cat("Wrote result report to ", report.file, "\n", sep="")
## interesting cases for training:

## only 1 sample, but several chunks => train on a few chunks,
## validate on the others.

## several samples, several chunks => train on a few chunks across all
## samples, test on some other chunks.

## theoretical: several samples, only 1 chunk => train on some
## samples, test on the other samples.

chunk.sample.list <- list()
for(sample.id in names(features.limits.list)){
  one.res <- features.limits.list[[sample.id]][[bases.per.bin.str]]
  for(chunk.id in names(one.res$limits)){
    limits <- one.res$limits[[chunk.id]]
    features <- one.res$features[rownames(limits), ,drop=FALSE]
    rownames(limits) <- rownames(features) <-
      paste(sample.id, chunk.id, rownames(limits))
    chunk.sample.list[[chunk.id]][[sample.id]] <-
      list(features=features, limits=limits)
  }
}

chunk.mats <- list()
for(chunk.id in names(chunk.sample.list)){
  one.chunk <- chunk.sample.list[[chunk.id]]
  data.types <- list()
  for(sample.id in names(one.chunk)){
    one.sample <- one.chunk[[sample.id]]
    for(data.type in c("features", "limits")){
      data.types[[data.type]][[sample.id]] <- one.sample[[data.type]]
    }
  }
  for(data.type in c("features", "limits")){
    chunk.mats[[chunk.id]][[data.type]] <-
      do.call(rbind, data.types[[data.type]])
  }
}

## chunk.mats is now a named list with an element for each chunk, that
## we can split on for cross-validation.

set.seed(1)
n.folds <- min(length(chunk.mats), 5)
folds <- 1:n.folds
fold.id <- sample(rep(folds, l=length(chunk.mats)))
eval.error.list <- list()
for(validation.fold in folds){
  sets <-
    list(validation=fold.id == validation.fold,
         train=fold.id != validation.fold)
  set.data <- list()
  for(set.name in names(sets)){
    is.set <- sets[[set.name]]
    for(data.type in c("features", "limits")){
      mats <- lapply(chunk.mats[is.set], "[[", data.type)
      set.data[[set.name]][[data.type]] <- do.call(rbind, mats)
    }
  }
  all.train.features <- set.data$train$features
  all.features <- rbind(all.train.features, set.data$validation$features)
  is.na.col <- apply(is.na(all.features), 2, any)
  is.inf.col <- apply(!is.finite(all.features), 2, any)
  ignore <- is.na.col | is.inf.col
  not.na <- all.train.features[, !ignore]
  fit <-
  regularized.interval.regression(features=not.na,
                                  limits=set.data$train$limits,
                                  L0=ncol(not.na)*1.5,
                                  calc.grad=calc.grad.list$square,
                                  calc.loss=calc.loss.list$square)
  for(set.name in names(sets)){
    eval.data <- set.data[[set.name]]
    pred.penalty.mat <- fit$predict(eval.data$features)
    for(complexity.i in 1:ncol(pred.penalty.mat)){
      pred.penalty <- pred.penalty.mat[, complexity.i]
      too.lo <- pred.penalty < eval.data$limits[,1]
      too.hi <- eval.data$limits[,2] < pred.penalty
      is.error <- too.lo | too.hi
      errors <- sum(is.error)
      ## TODO: compute real annotation error?
      gamma <- fit$gamma.seq[[complexity.i]]
      possible.errors <- length(is.error)
      percent.error <- errors/possible.errors
      eval.error.list[[paste(validation.fold, set.name, complexity.i)]] <-
        data.table(validation.fold, set.name, gamma,
                   errors, possible.errors, percent.error)
    }
  }
}

eval.error <- do.call(rbind, eval.error.list)
setkey(eval.error, set.name)
vali.error <- eval.error["validation"]
vali.error[, min.err := min(errors), by=validation.fold]
best.err <- vali.error[min.err == errors, ]
gamma.dt <- 
  best.err[, .(min.gamma=min(gamma),
               max.gamma=max(gamma)),
           by=validation.fold]

no.facets <- 
ggplot()+
  geom_vline(aes(xintercept=-log10(max.gamma)), data=gamma.dt)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_line(aes(-log10(gamma), percent.error,
                group=interaction(set.name, validation.fold),
                linetype=set.name),
            data=eval.error)
with.facets <-
  no.facets+
  facet_grid(validation.fold ~ ., labeller=function(var, val){
    paste("fold", val)
  })
modelSelection.base <- "figure-modelSelection.png"
modelSelection.png <- file.path(no.trailing, modelSelection.base)
png(modelSelection.png, units="in", res=200, width=15, height=4)
print(with.facets)
dev.off()

gamma.vals <- gamma.dt$max.gamma

best.gamma <- mean(gamma.vals)

## Fit the model on the entire data set.
all.features <- do.call(rbind, lapply(chunk.mats, "[[", "features"))
is.na.col <- apply(is.na(all.features), 2, any)
is.inf.col <- apply(!is.finite(all.features), 2, any)
ignore <- is.na.col | is.inf.col
not.na <- all.features[, !ignore]
all.limits <- do.call(rbind, lapply(chunk.mats, "[[", "limits"))
stopifnot(rownames(not.na) == rownames(all.limits))
learned.model <-
  smooth.interval.regression(features=not.na,
                             limits=all.limits,
                             gamma=best.gamma,
                             L0=ncol(not.na)*1.5,
                             calc.grad=calc.grad.list$square,
                             calc.loss=calc.loss.list$square)

## Save the learned resolution (bases.per.bin) and the inteval
## regression penalty function prediction model.
save(learned.model, bases.per.bin, file=out.RData)
