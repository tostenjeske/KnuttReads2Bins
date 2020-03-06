#!/usr/bin/env Rscript

##
## bbmergeParser.R - Read data from the BBmerge report
##
## Knutt.org/Knutt2Reads2Bins
#
# The log parsing approach is currently necessary, as the insert details
# doesn't report the ambigous merges, it only reports the flags P,I and F,
# F also contains the ambigous pairs. This script also adds the detected
# adapter from the FASTA file to the data file.
# Should be replaced with a Python script

options(warn=2)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ShortRead))

logfile <- c("output/ReadPrep/merging_tr/ChinaNM2MG/ChinaNM2MG_tr_merge.txt")
adaptersfile <- c("output/ReadPrep/merging_tr/ChinaNM2MG/ChinaNM2MG_tr_adapters.fa")

if(exists("snakemake")){
  logfile <- snakemake@input[["log"]]
  adaptersfile <- snakemake@input[["adapter"]]
  outputfile <- snakemake@output[[1]]
}

parseSummaryLogFile <- function(file){
  content <- readChar(file, file.info(file)$size)
  matches <- str_match_all(content,"Adapters counted:\\s+(\\d+)")
  adaptercount <- as.numeric(matches[[1]][,2])
  matches <- str_match_all(content,"Pairs:\\s+(\\d+)")
  pairstomerge <- as.numeric(matches[[1]][,2])
  matches <- str_match_all(content,"Joined:\\s+(\\d+)")
  joined <- as.numeric(matches[[1]][,2])
  matches <- str_match_all(content,"Ambiguous:\\s+(\\d+)")
  ambiguous <- as.numeric(matches[[1]][,2])
  matches <- str_match_all(content,"No Solution:\\s+(\\d+)")
  nosolution <- as.numeric(matches[[1]][,2])
  matches <- str_match_all(content,"Too Short:\\s+(\\d+)")
  tooshort <- as.numeric(matches[[1]][,2])
  matches <- str_match_all(content,"Insert range:\\s+(\\d+)\\s+-\\s+(\\d+)")
  shortestinsert <- as.numeric(matches[[1]][,2])
  longestinsert <- as.numeric(matches[[1]][,3])
  matches <- str_match_all(content,"50th percentile:\\s+(\\d+)")
  medianinsert <- as.numeric(matches[[1]][,2])
  data.frame(pairstomerge,joined,ambiguous,nosolution,tooshort,shortestinsert,medianinsert,longestinsert,adaptercount)
}

readAdapters <- function(file){
  fasta <- readFasta(file)
  reads <- as.character(sread(fasta))
  names(reads) <- id(fasta)
  data.frame(Read1_adapter=reads[["Read1_adapter"]],Read2_adapter=reads[["Read2_adapter"]])
}

logdata <- parseSummaryLogFile(logfile)
adapters <- readAdapters(adaptersfile)
result <- cbind(logdata,adapters)

write.table(result,outputfile,row.names = F,sep = "\t")

