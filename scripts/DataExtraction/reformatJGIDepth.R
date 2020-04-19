#!/usr/bin/env Rscript

##
## reformatJGIDepth.R - Add full contig description line to the JGI depth script data
##
## Knutt.org/KnuttReads2Bins

options(warn=2)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))

rawdetailsfile <- snakemake@input[["dat"]]
seqfile <- snakemake@input[["seq"]]
detailsfile <- snakemake@output[[1]]




dat <- fread(rawdetailsfile, sep="\t")
colnames(dat) <- c("contigid", "len", "avgcov", "avgcov2", "var")
dat[, avgcov2:=NULL]
seqs <- readFasta(seqfile)
seqs <- data.table(contig=as.character(id(seqs)))
seqs[, contigid:=tstrsplit(contig, " ", fixed=T)[[1]]]
setkey(dat, "contigid")
setkey(seqs, "contigid")
dat <- seqs[dat]
dat[, contigid:=NULL]
setcolorder(dat, "contig")
write.table(dat, detailsfile, sep="\t", row.names = F, quote=F)


