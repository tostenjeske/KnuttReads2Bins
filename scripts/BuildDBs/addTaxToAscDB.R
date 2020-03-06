#!/usr/bin/env Rscript
##
## addTaxToCustomQueryDB.R - Add tax columns to the databases with NCBI Ascs.
##
## Knutt.org/Knutt2Reads2Bins

# Very similar to addTaxToCustomQueryDB.R, it uses external grep to prefilter
# the database asc to taxid file to reduce the memory footprint during
# merging.

options(warn=2)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

ascfield = "sseqid"
taxcols = c("superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species")

if (exists("snakemake")) {
  database_file <- snakemake@input[["db"]]
  asctax_file <- snakemake@input[["asctax"]]
  translator_file <- snakemake@input[["translator"]]
  output_file <- snakemake@output[["db"]]
  krona_file <- snakemake@output[["krona"]]
  threads <- snakemake@threads
}else{
  database_file <- "databases/readanno/prot/raw/hyddb.tsv"
  translator_file <- "databases/ncbitax/ncbi_tax.RData"
  asctax_file <- "databases/ncbitax/ncbi_taxid_acc_prot.dmp"
  output_file <- "test.tsv"
  krona_file <- "test2.tsv"
  threads <- detectCores()
}

setDTthreads(threads)
options(mc.cores = threads)

database <- fread(database_file)
fieldid <- which(colnames(database) == ascfield)
# Use grep to filter only the asc entries that are actually needed
# Reduces the memory footprint of the merge operation.
cmd <- paste(
  "export LC_ALL=C && cut -f",
  fieldid,
  database_file,
  " | tr -d '\"' |  grep -F -f -",
  asctax_file
)
taxidtable <- fread(
  cmd = cmd,
  header = F,
  col.names = c(ascfield, "taxid")
)
database <- taxidtable[database, on = ascfield]
taxids <- database[["taxid"]]

# Same as addTaxToCustomQueryDB:
load(translator_file)
tax <- lookup(taxids, taxlevels = taxcols, threads = threads)
database <- cbind(database, tax)
database[, query := NULL]
fwrite(database, file = output_file, sep = "\t")
kronadata <- cbind(freq = 1, All = "All", database[, taxcols, with = F])
fwrite(kronadata, krona_file, col.names = F, sep = "\t")
