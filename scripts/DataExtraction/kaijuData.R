#!/usr/bin/env Rscript

##
## kaijuData.R - Reformat kaijus result files a little bit better
##
## Knutt.org/KnuttReads2Bins

options(warn = 2)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))

datafile <- snakemake@input$kaiju
seqfile <- snakemake@input$seq
outputfile <- snakemake@output[[1]]

data <- fread(datafile, fill=T, header=F, sep="\t")
colnames(data) <- c("classified","qid","taxid","lengthorscore","match_taxids","match_ascs","match_fragments","tax")
setkey(data, qid)
seqs <- data.table(qname=as.character(id(readFastq(seqfile))))
seqs[, qid:=tstrsplit(qname, " ", fixed=T)[[1]]]
setkey(seqs, qid)
data <- seqs[data]
data[, qid:=NULL]

fwrite(data, outputfile, sep="\t", quote=F)