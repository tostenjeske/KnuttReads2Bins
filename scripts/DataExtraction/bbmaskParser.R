#!/usr/bin/env Rscript

##
## bbmergeParser.R - Read data from the BBmask report
##
## Knutt.org/KnuttReads2Bins
#


options(warn=2)
suppressPackageStartupMessages(library(stringr))

logfile <- c("output/ReadPrep/Merging_tr/ChinaNM2MG/ChinaNM2MG_merge_tr_merged_mask.log")

if(exists("snakemake")){
  logfile <- snakemake@input[[1]]
  outputfile <- snakemake@output[[1]]
}

parseSummaryLogFile <- function(file){
  content <- readChar(file, file.info(file)$size)
  matches <- str_match_all(content,"Total Bases Masked:\\s+(\\d+)")
  total_masked <- as.numeric(matches[[1]][,2])
  data.frame(total_masked)
}

logdata <- parseSummaryLogFile(logfile)

write.table(logdata,outputfile, row.names = F, sep = "\t", quote = FALSE) 
