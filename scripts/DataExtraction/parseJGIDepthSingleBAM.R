#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

logfiles <- c("")
outputfile <- "testing.tsv"

if(exists("snakemake")){
  logfiles <- snakemake@input[["files"]]
  samplenames <- snakemake@params[["samplenames"]]
  outputfile <- snakemake@output[["details"]]
}

names(logfiles) <- samplenames

sumdf <- rbindlist(lapply(samplenames,function(s)cbind(sample=s,fread(logfiles[[s]],sep="\t",fill=T,col.names = c("contigName","contigLen","totalAvgDepth","bamvar"),select = c(1:3,5)))),fill=T)
fwrite(sumdf,outputfile,sep="\t",row.names = F)


