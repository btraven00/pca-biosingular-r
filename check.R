#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(BiocSingular)
  library(HDF5Array)
  library(DelayedArray)
  library(Matrix)
  library(argparser)
  library(jsonlite)
  library(RhpcBLASctl)
  library(data.table)
})
cat("OK\n")
