#!/usr/bin/env Rscript

packageList <- c(
  "circlize",
  "clusterProfiler",
  "ComplexHeatmap",
  "DESeq2",
  "DOSE",
  "dplyr",
  "fdrtool",
  "fgsea",
  "ggplot2",
  "ggrepel",
  "GO.db",
  "grid",
  "GSEABase",
  "limma",
  "pvclust"
)

packagesToInstall <- setdiff(packageList, installed.packages()[,"Package"])
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(packagesToInstall, updates = TRUE)
