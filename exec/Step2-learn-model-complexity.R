library(data.table)
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
names(RData.files) <- sub("[.]RData$", "", basename(RData.files))
out.RData <- file.path(no.trailing, "learned.model.RData")
cat("Output file:\n", out.RData, "\n", sep="")

error.list <- list()
features.limits.list <- list()
for(sample.id in names(RData.files)){
  RData.file <- RData.files[[sample.id]]
  load(RData.file)
  error.list[[sample.id]] <- data.table(sample.id, errors)
  features.limits.list[[sample.id]] <- features.limits
}

errors <- do.call(rbind, error.list)
res.errors <- 
  errors[,
         .(weighted.error=sum(weighted.error),
           total.weight=sum(total.weight)),
         by=bases.per.bin]

errPlot <- 
ggplot(res.errors, aes(bases.per.bin, weighted.error))+
  geom_line()+
  scale_x_log10()+
  geom_point()

error.pdf <- file.path(no.trailing, "figure-weightedError.pdf")
pdf(error.pdf)
print(errPlot)
dev.off()

## consider only resolutions with max weights.
max.weight.value <- max(res.errors$total.weight)
max.weight <- res.errors[total.weight == max.weight.value, ]

## There could be several minimum error resolutions.
min.err.value <- min(max.weight$weighted.error)
min.err <- max.weight[weighted.error == min.err.value, ]

## The largest resolution will take the least amount of CPU time.
bases.per.bin <- max(min.err$bases.per.bin)
bases.per.bin.str <- paste(bases.per.bin)

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
      ## TODO: compute real annotation error.
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
modelSelection.pdf <- file.path(no.trailing, "figure-modelSelection.pdf")
pdf(modelSelection.pdf)
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
