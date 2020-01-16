#!/bin/Rscript

if (!require(plyr, quietly = TRUE)) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}

absDiff <- function(arg1,arg2) {
  return(sqrt((arg1-arg2)^2))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Input data file name not provided\n", call. = FALSE)
}

in_original <- args[1]
in_acr <- args[2]

simTableOriginal <- read.table(in_original, col.names = c("xpos", "ypos", "density", "velx", "vely"))
simTableAcr <- read.table(in_acr, col.names = c("xpos", "ypos", "density", "velx", "vely"))

if (!all.equal(simTableOriginal[,c("xpos", "ypos")],simTableAcr[,c("xpos", "ypos")])) {
  stop("The data grid of the two files does not match\n", call. = FALSE)
}

maxdens <- max(simTableOriginal$density)
mindens <- min(simTableOriginal$density)
range  <- absDiff(maxdens,mindens)

summaryData <- data.frame(unclass(summary((absDiff(simTableAcr$density, simTableOriginal$density)) / range * 100)))
names(summaryData) <- ""
print(summaryData, digits = 3, zero.print = "0")
cat(sprintf("Range %.3f\n", range))
