#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))
render(normalizePath(snakemake@input[["script"]]), output_file=normalizePath(snakemake@output[[1]],mustWork = FALSE), knit_root_dir=normalizePath(getwd()), quiet=TRUE, params = list(rmd=normalizePath(snakemake@input[["script"]])),output_dir = normalizePath(dirname(snakemake@output[[1]])))
