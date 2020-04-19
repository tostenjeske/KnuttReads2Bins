#!/usr/bin/env Rscript

##
## mergeDBintoBLASTXres.R - Add the database info
##
## Knutt.org/KnuttReads2Bins
#
# Merge the database info file into the BLASTX results and filter
# the best hit for every entry

options(warn=2)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))


if(exists("snakemake")){
  blastxcolnames <- snakemake@params[["blastxcolnames"]]
  inputfile <- snakemake@input[["blastxres"]]
  datafile <- snakemake@input[["datafile"]]
  seqfile <- snakemake@input[["seq"]]
  outputfile <- snakemake@output[[1]]
  threads <- snakemake@threads
}

setDTthreads(threads)



blastres <- fread(inputfile ,header=F)
colnames(blastres) <- blastxcolnames
data <- fread(datafile)
seqs <- data.table(qname=as.character(id(readFastq(seqfile))))
seqs[, qseqid:=tstrsplit(qname, " ", fixed=T)[[1]]]
setkey(seqs, qseqid)

# Filter the best result for each query
setorder(blastres, qseqid, evalue)
blastres <- blastres[, head(.SD, 1), by=qseqid]

# Add db data
colnames(data)[1] <- "sseqid"
blastres <- data[blastres, on="sseqid"]
setcolorder(blastres, c(blastxcolnames))

# Add full read title
blastres <- seqs[blastres, on="qseqid"]
setcolorder(blastres, c("qname","stitle"))
setnames(blastres,c("stitle"),c("sname"))
blastres[, qseqid:=NULL]
blastres[, sseqid:=NULL]

fwrite(blastres, file=outputfile, sep="\t", quote=F)
