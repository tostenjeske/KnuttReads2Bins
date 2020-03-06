#!/usr/bin/env Rscript

##
## addTaxToCustomQueryDB.R - Add tax columns to the databases from Uniprot
##
## Knutt.org/Knutt2Reads2Bins

# Uses the "Organism ID" field in the input file and the pregenerated
# ncbitaxhelper file. Also produces a file for Krona.

options(warn=2)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))


taxfield <- "Organism ID"
taxcols <-
  c("superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species")

if (exists("snakemake")) {
  database_file <- snakemake@input[["db"]]
  translator_file <- snakemake@input[["translator"]]
  output_file <- snakemake@output[["db"]]
  krona_file <- snakemake@output[["krona"]]
  threads <- snakemake@threads
}else{
  database_file <- "pilus.tsv"
  translator_file <- "reference_data/ncbitax/ncbi_tax.RData"
  output_file <- "test.tsv"
  krona_file <- "test2.tsv"
  threads <- detectCores()
}

setDTthreads(threads)
options(mc.cores = threads)

database <- fread(database_file)
taxids <- database[[taxfield]]


load(translator_file)
tax <-
  lookup(taxids, taxlevels = taxcols, threads = threads)
database <- cbind(database, tax)
database[, query := NULL]
fwrite(database, file = output_file, sep = "\t")
kronadata <- cbind(freq = 1, All = "All", database[, taxcols, with = F])
fwrite(kronadata, krona_file, col.names = F, sep = "\t")

