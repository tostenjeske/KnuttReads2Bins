#!/usr/bin/env Rscript

##
## binMap.R - Create a bin to contig mapping file
##
## Knutt.org/KnuttReads2Bins
#

options(warn=2)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))

bindir <- snakemake@input[["bindir"]]
additionalfiles <- snakemake@input[["additional"]]
binpattern <- "(.+)\\.fa"
assembly <- snakemake@input[["assembly"]]
outputfile <- snakemake@output[[1]]

binfiles <- c(list.files(bindir, binpattern, full.names=T), additionalfiles)
binfiles <- gsub("//", "/", binfiles, fixed=T)

readBin <- function(file){
    contigid <- tstrsplit(as.character(id(readFasta(file))), " ", fixed=T)
    if(length(contigid)==0)
        return(data.table(file=character(), contigid=character()))
    data.table(file, contigid=contigid[[1]])
}

ids <- rbindlist(lapply(binfiles, readBin))
ids[, bin:=sub(binpattern, "\\1", basename(file))]
seqs <- data.table(contig=as.character(id(readFasta(assembly))))
seqs[, contigid:=tstrsplit(contig, " ", fixed=T)[[1]]]
setkey(seqs, contigid)
setkey(ids, contigid)
ids <- seqs[ids]
ids[, file:=NULL]
ids[, contigid:=NULL]
setcolorder(ids, c("bin","contig"))
fwrite(ids, outputfile, sep="\t", quote=F)