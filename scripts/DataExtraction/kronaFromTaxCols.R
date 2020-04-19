#!/usr/bin/env Rscript

##
## kronaFromTaxCols.R - Convert a file with taxcols into a krona compatible one
##
## Knutt.org/KnuttReads2Bins
#
# Read the given file and construct a krona file.

options(warn=2)
suppressPackageStartupMessages(library(data.table))

taxcols <- c("superkingdom","phylum","class","order","family","genus","species")
inputfile <- "/run/user/2496/gvfs/sftp:host=highmem1/data1/2019_LJ_data/bachelor_thesis/metadiet/output/ReadAnnotation/Larrelt/Larrelt_prot_pilus_blast.tsv"
outputfile <- "testing.tsv"


if(exists("snakemake")){
  inputfile <- snakemake@input[[1]]
  outputfile <- snakemake@output[[1]]
}


dataset <- fread(inputfile)
kronadata <- cbind(freq=1, All="All", dataset[, taxcols, with=F])
kronadata[is.na(kronadata)] <- "Unavailable"
fwrite(kronadata, outputfile, col.names=F, sep="\t")

