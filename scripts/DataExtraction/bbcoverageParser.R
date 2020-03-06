#!/usr/bin/env Rscript

##
## bbcoverageParser.R - Get coverage statistics from pileup.sh script
##
## Knutt.org/Knutt2Reads2Bins
#
# The log contains information on the total read coung and proper pairing
# of the pairs, this information can't be extracted from the detail file.


options(warn=2)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))

logfile <- "ChinaNM2MG_pileup_cov.log"
outputfile <- "testing.tsv"

if(exists("snakemake")){
  logfile <- snakemake@input[[1]]
  outputfile <- snakemake@output[[1]]
}

patternsuffix = ":\\s+\\t(\\d+\\.?\\d*)"
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
  mappedbases <- extract("Mapped bases")
  refscaffs <- extract("Ref scaffolds")
  refbases <- extract("Ref bases")
  percproperpairs <- extract("Percent proper pairs")
  avcov <- extract("Average coverage")
  stddev <- extract("Standard deviation")
  scaffwithanycov <- extract("Percent scaffolds with any coverage")
  refbasescov <- extract("Percent of reference bases covered")

  df<-data.frame(reads,mappedreads,mappedbases,refscaffs,refbases,percproperpairs,avcov,stddev,scaffwithanycov,refbasescov)
  df
}

sumdf=parseSummaryLogFile(logfile)
write.table(sumdf,outputfile,sep="\t",row.names = F)
