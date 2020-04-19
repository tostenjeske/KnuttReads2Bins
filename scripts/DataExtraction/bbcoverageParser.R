#!/usr/bin/env Rscript

##
## bbcoverageParser.R - Get coverage statistics from pileup.sh script
##
## Knutt.org/KnuttReads2Bins
#
# The log contains information on the total read coung and proper pairing
# of the pairs, this information can't be extracted from the detail file.


options(warn=2)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))

logfile <- snakemake@input[["log"]]
rawdetailsfile <- snakemake@input[["details"]]
seqfile <- snakemake@input[["seq"]]
summaryfile <- snakemake@output[["summary"]]
detailsfile <- snakemake@output[["details"]]

patternsuffix <- ":\\s+\\t(\\d+\\.?\\d*)"
pat <- function(prefix){
 paste0(prefix,patternsuffix)
}

parseSummaryLogFile <- function(file){
  content <- readChar(file, file.info(file)$size)

  extract <- function(prefix){
    matches <- str_match_all(content,pat(prefix))
    as.numeric(matches[[1]][,2])
  }

  reads <- extract("Reads")
  mappedreads <- extract("Mapped reads")
  mappedbp <- extract("Mapped bases")
  contigs <- extract("Ref scaffolds")
  contigbp <- extract("Ref bases")
  properpairsperc <- extract("Percent proper pairs")
  avgcov <- extract("Average coverage")
  stddev <- extract("Standard deviation")
  contigswithanycovperc <- extract("Percent scaffolds with any coverage")
  bpswithanycovperc <- extract("Percent of reference bases covered")

  df<-data.frame(reads,mappedreads,mappedbp,contigs,contigbp,properpairsperc,avgcov,stddev,contigswithanycovperc,bpswithanycovperc)
  df
}

sumdf <- parseSummaryLogFile(logfile)
write.table(sumdf, summaryfile, sep="\t", row.names = F, quote=F)


newcols <- c("contigid", "avgcov", "len", "refgc", "covperc", "plus_reads", "minus_reads", "readgc", "mediancov", "stddev", "covbases")
dat <- fread(rawdetailsfile, sep="\t")
setnames(dat, c("#ID", "Avg_fold", "Length", "Ref_GC", "Covered_percent", "Plus_reads", "Minus_reads", "Read_GC", "Median_fold", "Std_Dev", "Covered_bases"), newcols)
newcols[[1]] <- c("contig")
dat[, refgc:=NULL]
seqs <- readFasta(seqfile)
seqs <- data.table(contig=as.character(id(seqs)), refgc=round(letterFrequency(sread(seqs), letters="CG", as.prob=T)[, 1], 3))
seqs[, contigid:=tstrsplit(contig, " ", fixed=T)[[1]]]
setkey(seqs, contigid)
setkey(dat, contigid)
dat <- seqs[dat]
dat[, contigid:=NULL]
setcolorder(dat, newcols)
write.table(dat, detailsfile, sep="\t", row.names = F, quote=F)


