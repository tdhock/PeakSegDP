if(!require(PeakSegDP))
  devtools::install_github("tdhock/PeakSegDP")
if(!require(PeakError))
  devtools::install_github("tdhock/PeakError")

## Make and run qsub scripts for all steps of the PeakSegDP pipeline.
R.bin <- R.home("bin")
Rscript <- file.path(R.bin, "Rscript")
labels.txt.file <- # interactive default for debugging.
  system.file(file.path("exampleData", "manually_annotated_region_labels.txt"),
              package="PeakSegDP")
argv <- commandArgs(trailingOnly=TRUE)
if(length(argv) != 1){
  stop("usage: AllSteps_qsub.R path/to/labels.txt
where there are path/to/*/*.bedGraph files")
}
labels.txt.file <- normalizePath(argv[1], mustWork=TRUE)
data.dir <- dirname(labels.txt.file) 

## Step0 does not take very long so we can just do it live in this
## script.
Step0 <-
  system.file(file.path("exec", "Step0-convert-labels.R"),
              package="PeakSegDP")
cmd <- paste(Rscript, Step0, labels.txt.file)
system(cmd)

## Step1 takes a while so we run it with qsub.
data.dir <- normalizePath(data.dir, mustWork=TRUE)#NEEDS TO BE ABSOLUTE!
Step1 <-
  system.file(file.path("exec", "Step1-segment-labeled-regions.R"),
              package="PeakSegDP")
labels.files <- Sys.glob(file.path(data.dir, "*", "*_labels.bed"))
residual.qsub.id.list <- list()
for(labels.file in labels.files){
  base.path <- sub("_labels.bed$", "", labels.file)
  sample.id <- basename(base.path)
  bedGraph.path <- paste0(base.path, ".bedGraph")
  residuals.base <- paste0(base.path, "_residuals")
  script.txt <-
    paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00                      
#PBS -A bws-221-ae                             
#PBS -o ", residuals.base, ".out
#PBS -e ", residuals.base, ".err
#PBS -V                                        
#PBS -N ", sample.id, "
", Rscript, " ", Step1, " ", bedGraph.path, " ", labels.file)
  script.file <- paste0(residuals.base, ".sh")
  cat(script.txt, file=script.file)
  cmd <- paste("qsub", script.file)
  qsub.out <- system(cmd, intern=TRUE)
  residual.qsub.id <- sub("[.].*", "", qsub.out)
  cat("submitted job ", residual.qsub.id, "\n", sep="")
  residual.qsub.id.list[[script.file]] <- residual.qsub.id
}

## Step2 learning depends on Step1.
Step2 <-
  system.file(file.path("exec", "Step2-learn-model-complexity.R"),
              package="PeakSegDP")
residual.qsub.id.txt <- paste(residual.qsub.id.list, collapse=":")
learned.base <- file.path(data.dir, "learned.model")
script.txt <-
  paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00                      
#PBS -A bws-221-ae
#PBS -W depend=afterok:", residual.qsub.id.txt, "
#PBS -o ", learned.base, ".out
#PBS -e ", learned.base, ".err
#PBS -V                                        
#PBS -N learned.model
", Rscript, " ", Step2, " ", data.dir)
script.file <- paste0(learned.base, ".sh")
cat(script.txt, file=script.file)
cmd <- paste("qsub", script.file)
qsub.out <- system(cmd, intern=TRUE)
learned.qsub.id <- sub("[.].*", "", qsub.out)
cat("submitted job ", learned.qsub.id, "\n", sep="")

## Step3: genome-wide peak prediction.
Step3 <-
  system.file(file.path("exec", "Step3-predict-peaks.R"),
              package="PeakSegDP")
bedGraph.files <- Sys.glob(file.path(data.dir, "*", "*.bedGraph"))
RData.file <- file.path(data.dir, "learned.model.RData")
peak.qsub.id.list <- list()
for(bedGraph.file in bedGraph.files){
  peaks.base <- sub("[.]bedGraph$", "_peaks", bedGraph.file)
  script.txt <-
    paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00                      
#PBS -A bws-221-ae
#PBS -W depend=afterok:", learned.qsub.id, "
#PBS -o ", peaks.base, ".out
#PBS -e ", peaks.base, ".err
#PBS -V                                        
#PBS -N ", basename(peaks.base), "
", Rscript, " ", Step3, " ", RData.file, " ", bedGraph.file)
  script.file <- paste0(peaks.base, ".sh")
  cat(script.txt, file=script.file)
  cmd <- paste("qsub", script.file)
  qsub.out <- system(cmd, intern=TRUE)
  peak.qsub.id <- sub("[.].*", "", qsub.out)
  cat("submitted job ", peak.qsub.id, "\n", sep="")
  peak.qsub.id.list[[script.file]] <- peak.qsub.id
}

## Step4: model visualization and peak clustering.
Step4 <-
  system.file(file.path("exec", "Step4-cluster-peaks.R"),
              package="PeakSegDP")
peak.qsub.id.txt <- paste(peak.qsub.id.list, collapse=":")
predicted.base <- file.path(data.dir, "predicted-figures")
script.txt <-
  paste0("#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=12:00:00                      
#PBS -A bws-221-ae
#PBS -W depend=afterok:", peak.qsub.id.txt, "
#PBS -o ", predicted.base, ".out
#PBS -e ", predicted.base, ".err
#PBS -V                                        
#PBS -N predicted-figures
", Rscript, " ", Step4, " ", data.dir)
script.file <- paste0(predicted.base, ".sh")
cat(script.txt, file=script.file)
cmd <- paste("qsub", script.file)
qsub.out <- system(cmd, intern=TRUE)
cluster.qsub.id <- sub("[.].*", "", qsub.out)
cat("submitted job ", cluster.qsub.id, "\n", sep="")

