#!/usr/bin/env Rscript

##
## CombineBBMapAndSinaData.R - Add the results from BBmap to the SINA file
##
## Knutt.org/KnuttReads2Bins

# This script parses the BAM file and merges this data with the SINA csv
# file
options(warn = 2)

suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(stringi))

bam.file <- "/media/data/ChinaMM2MG_bbmap_16S.bam"
sinahits.file <- "/media/data/ChinaMM2MG_sina_hits.csv"
outputfile <- "test.tsv"

if(exists("snakemake")){
  bam.file <- snakemake@input[["bam"]]
  sinahits.file <- snakemake@input[["sina"]]
  outputfile <- snakemake@output[[1]]
}



bam2Df <- function(bam.file){
  bam <- scanBam(bam.file, param=ScanBamParam(what=c("qname","rname","strand","pos","qwidth","mapq","cigar"), tag=c( "MD", "XM","NM"), flag=scanBamFlag(isUnmappedQuery=FALSE)))
  as.data.frame(bam)
}

extractFastaIDfromName <- function(names){
  split <- stri_split_fixed(names, pattern=" ", n=2)
  ids <- sapply(split, "[[", 1)
  descr <- sapply(split,function(elements)ifelse(length(elements)>1, elements[[2]], NA))
  data.frame(id=ids, descr=descr)
}


result.bbmap <- bam2Df(bam.file)
result.bbmap[c("qid", "qdescr")] <- extractFastaIDfromName(result.bbmap$qname)
result.bbmap[c("rid", "rdescr")] <- extractFastaIDfromName(result.bbmap$rname)

result.sina <- read.csv(sinahits.file)
colnames(result.sina)[which(colnames(result.sina)=="name")] <- "qid"

result <- merge(result.bbmap,result.sina, by=c("qid"), all=T)
result[c("qid", "qdescr", "rid", "rdescr", "full_name")] <- NULL
colnames(result)[which(colnames(result)=="rname")] <- "sname"

write.table(result,outputfile, row.names = F, sep = "\t", quote=F)
