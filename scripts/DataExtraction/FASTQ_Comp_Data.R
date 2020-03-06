#!/usr/bin/env Rscript

##
## FASTQ_Comp_Data.R - Compare the input and output of a FASTQ step
##
## Knutt.org/Knutt2Reads2Bins

# Takes two FASTQ files and tries to find the matching entries (same id)
# and reports changes in average quality,length and entry count.

options(warn = 2)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))



if (exists("snakemake")) {
  reads_before <- snakemake@input[["before"]]
  reads_after <- snakemake@input[["after"]]
  threads <- snakemake@threads
  summary_file <- snakemake@output[[1]]
} else{
  reads_before <- "input/ChinaNM2MG_R1.fastq.gz"
  reads_after <-
    "output/ReadPrep/trimming/ChinaNM2MG/ChinaNM2MG_R1_adptr_tr.fastq.gz"
  summary_file <- "test.tsv"
  threads <- 4
}

options(mc.cores = threads)
setDTthreads(threads = threads)

compareReads <-
  function(before.fastq.path,
           after.fastq.path,
           verbose = F,
           entries = 1000000,
           details_points = 512) {
    before.fastq.stream <-
      FastqStreamer(before.fastq.path, verbose = verbose, n = entries)
    after.fastq.stream <-
      FastqStreamer(after.fastq.path, verbose = verbose, n = entries)
    on.exit(close(before.fastq.stream), add = T)
    on.exit(close(after.fastq.stream), add = T)
    before.fastq.going <- T
    after.fastq.going <- T
    unmatchedentries.before <-
      data.table(id = character(),
                 length = integer(),
                 avgquality = double())
    unmatchedentries.after <-
      data.table(id = character(),
                 length = integer(),
                 avgquality = double())
    lengthchanges <- integer()
    qualchanges <- double()
    while (before.fastq.going || after.fastq.going) {
      if (before.fastq.going) {
        if (verbose)
          print("Reading before ...")
        before.entry <- yield(before.fastq.stream)
        if (length(before.entry) == 0)
          before.fastq.going = F
        else{
          if (verbose)
            print("Crunching before ...")
          unmatchedentries.before = rbind(
            unmatchedentries.before,
            data.table(
              id = as.character(id(before.entry)),
              length = width(before.entry),
              avgquality = apply(as(quality(
                before.entry
              ), "matrix"), 1, mean, na.rm = T)
            )
          )
        }
      }
      if (after.fastq.going) {
        if (verbose)
          print("Reading after ...")
        after.entry <- yield(after.fastq.stream)
        if (length(after.entry) == 0)
          after.fastq.going = F
        else{
          if (verbose)
            print("Crunching after ...")
          unmatchedentries.after <-
            rbind(
              unmatchedentries.after,
              data.table(
                id = as.character(id(after.entry)),
                length = width(after.entry),
                avgquality = apply(as(quality(
                  after.entry
                ), "matrix"), 1, mean, na.rm = T)
              )
            )
        }
      }
      if (verbose)
        print("Merging ...")
      newmatchedentries <-
        merge(
          unmatchedentries.before,
          unmatchedentries.after,
          by = "id",
          suffixes = c(".before", ".after")
        )
      if (verbose)
        print("Finding unmerged entries ...")
      unmatchedentries.before <-
        unmatchedentries.before[!id %in% newmatchedentries$id, ]
      unmatchedentries.after <-
        unmatchedentries.after[!id %in% newmatchedentries$id, ]
      if (verbose)
        print("Appending new values ...")
      if (nrow(newmatchedentries) > 0) {
        lengthchanges <-
          c(lengthchanges, newmatchedentries[, length.after - length.before])
        qualchanges <-
          c(qualchanges, newmatchedentries[, avgquality.after - avgquality.before])
      }
      #    before.fastq.going = F
      #    after.fastq.going = F
    }
    if (verbose)
      print("Calculating ...")
    newentries <- nrow(unmatchedentries.after)
    deletedentries <- nrow(unmatchedentries.before)
    data.table(
      minReadLenChange = min(lengthchanges),
      quantile25ReadLenChange = quantile(lengthchanges, 0.25),
      medianReadLenChange = median(lengthchanges),
      quantile75ReadLenChange = quantile(lengthchanges, 0.75),
      maxReadLenChange = max(lengthchanges),
      meanReadLenChange = mean(lengthchanges),
      minAvgQualChange = min(qualchanges),
      quantile25AvgQualChange = quantile(qualchanges, 0.25),
      medianAvgQualChange = median(qualchanges),
      quantile75AvgQualChange = quantile(qualchanges, 0.75),
      maxAvgQualChange = max(qualchanges),
      meanAvgQualChange = mean(qualchanges),
      newentries,
      deletedentries,
      matchingentries = length(lengthchanges)
    )
  }


summary <- compareReads(reads_before, reads_after)

fwrite(summary, summary_file, sep = "\t")
