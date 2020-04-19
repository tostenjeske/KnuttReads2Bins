#!/usr/bin/env Rscript

##
## sourmashSummaryToKrona.R - Fix sourmash abdundances for Krona
##
## Knutt.org/KnuttReads2Bins
#
# Read the given file and construct a krona file.

options(warn=2)
suppressPackageStartupMessages(library(data.table))


inputfile <- snakemake@input[[1]]
outputfile <- snakemake@output[[1]]
taxcols <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain")

dat <- fread(inputfile)
dat[is.na(dat)] <- ""
dat[, depth:=apply(.SD, 1, function(row)sum(row!="")), .SD=taxcols]

for(curdepth in 2:length(taxcols)){
  subtract <- dat[depth==curdepth, ]
  if(nrow(subtract)==0)
    next
  toreplace <- taxcols[curdepth:length(taxcols)]
  subtract[, c(toreplace):=""]
  subtract <- subtract[, .(minus=sum(count)), by=taxcols]
  dat <- merge(dat, subtract, on=taxcols, all.x=TRUE)
  dat[is.na(minus), minus:=0]
  dat[, count:=count-minus]
  dat[, minus:=NULL]
  dat <- dat[count!=0,]
}

dat[, depth:=NULL]
setcolorder(dat, c("count"))
dat[,(taxcols):=lapply(.SD, function(col)sub(".__(.+)","\\1",col)),.SD=taxcols]
fwrite(dat, outputfile, col.names=F, sep="\t", quote=FALSE)
system(paste0("sed 's/[[:blank:]]*$//' -i ", outputfile))