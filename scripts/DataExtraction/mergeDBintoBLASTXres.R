#!/usr/bin/env Rscript

##
## mergeDBintoBLASTXres.R - Add the database info
##
## Knutt.org/Knutt2Reads2Bins
#
# Merge the database info file into the BLASTX results and filter
# the best hit for every entry

options(warn=2)
suppressPackageStartupMessages(library(data.table))

taxcols=c("superkingdom","phylum","class","order","family","genus","species")
blastxcolnames = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","staxids","sscinames")
inputfile = "/run/user/2496/gvfs/sftp:host=highmem1/data1/2019_LJ_data/bachelor_thesis/metadiet/output/ReadAnnotation/Larrelt/Larrelt_prot_pilus_blast.tsv"
datafile = "/run/user/2496/gvfs/sftp:host=highmem1/data1/2019_LJ_data/bachelor_thesis/metadiet/output/databases/readanno/prot/pilus.tsv"
outputfile = "testing.tsv"
threads=1


if(exists("snakemake")){
  blastxcolnames = snakemake@params[["blastxcolnames"]]
  inputfile = snakemake@input[["blastxres"]]
  datafile = snakemake@input[["datafile"]]
  outputfile = snakemake@output[[1]]
  threads = snakemake@threads
}

setDTthreads(threads)



blastres = fread(inputfile,header=F)
colnames(blastres) = blastxcolnames
data = fread(datafile)

# Filter the best result for each query
setorder(blastres,qseqid,evalue)
blastres = blastres[, head(.SD, 1), by=qseqid]
colnames(data)[[1]] = "sseqid"
blastres = data[blastres,on="sseqid"]
fwrite(blastres,file=outputfile,sep="\t")
