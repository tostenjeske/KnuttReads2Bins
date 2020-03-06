#!/usr/bin/env Rscript

##
## FASTQ_Data.R - Generate plotable data from a FASTQ file
##
## Knutt.org/Knutt2Reads2Bins

# Outputs on overview and a detail file
#
# The detail file combines multiple x,y,(z) datasets:
# gccontent_density: GC content density distribution (over all reads)
# avg_seq_qual_dens: Density of the quality averages
# cycle_qual_counts: Occurence(z) of the quality values (y) in each
#                    cycle (x).
#
# readlength_density: Density distribution of the read lengths
# 
# There is something strange happening with the datatypes in the details.
# Writing and reading the toplot data.table fixes the incorrect numericals
# during plotting.


options(warn=2)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))


if(exists("snakemake")){
  input_file <- snakemake@input[[1]]
  threads <- snakemake@threads  
  overview_file <- snakemake@output[["overview"]]
  plot_data_file <- snakemake@output[["toplot"]]
}else{
  # List of named (sample names) lists with the R1 and R2 files
  input_file="input/ChinaNM2MG_R1.fastq.gz"
  threads <- 4
  overview_file <- "overview.tsv"
  plot_data_file <- "toplot.tsv"
}

options(mc.cores=threads)
setDTthreads(threads)

# Get the overview and plot dataset for a single file
# sample: Sample name
# readdirection: R1/R2
# verbose: Print progress
# setsize: Number of sequences to draw in on iteration
# details_points: Number of points in (cumulative) density plots
getOverviewAndDetails <- function(fastq.path,verbose=F,setsize=1000000,details_points=512){
  # sample=names(raw_reads)[1];readdirection="R1"
  # fastq.path=raw_reads[[sample]][["R1"]];verbose=T
  # setsize=1000000;details_points=512
  fastq.stream <- FastqStreamer(fastq.path,verbose = verbose,n = setsize)
  on.exit(close(fastq.stream),add = T)
  readlengths <- c()
  gccontents <- c()
  avg_seq_quals <- c()
  # Frequency table of the cycle, quality combinations over all reads
  cycle_quals <- NULL
  repeat {
    entry = yield(fastq.stream)
    if(length(entry) == 0)
      break
    readlengths <- c(readlengths, width(entry))
    gccontents <- c(gccontents,unlist(letterFrequency(sread(entry),"GC",as.prob = T)))
    # Rows: Reads, Cols: Cycles
    qualmat<-as(quality(entry),"matrix") 
    if(verbose)
      print("Averaging ...")
    avg_seq_quals <- unlist(lapply(1:nrow(qualmat),function(i)mean(qualmat[i,],na.rm=T)))
    if(verbose)
      print("Counting ...")
    qual_occurence_chunk <- lapply(1:ncol(qualmat),function(i)data.table(cycle=i,as.data.frame(table(qual=qualmat[,i]))))
    if(verbose)
      print("Combining ...")
    if(!is.null(cycle_quals)){
      cycle_quals <- rbindlist(c(qual_occurence_chunk,list(cycle_quals)))[,.(Freq=sum(Freq)),by=c("cycle","qual")]
    }
    else
      cycle_quals<-rbindlist(qual_occurence_chunk)[,.(Freq=sum(Freq)),by=c("cycle","qual")]
  }
  if(verbose)
    print("Results ...")
  overview <- data.table(minReadLen=min(readlengths),quantile25ReadLen=quantile(readlengths,0.25),medianReadLen=median(readlengths),
                        quantile75ReadLen=quantile(readlengths,0.75),maxReadLen=max(readlengths),meanReadLen=mean(readlengths),reads=length(readlengths),basepairs=sum(readlengths))
  gccontent <- data.table(type="gccontent_density",as.data.table((density(gccontents,from=min(gccontents),to=max(gccontents),n=details_points)[c("x","y")])))
  lendens <- data.table(type="readlength_density",as.data.table((density(readlengths,from=min(readlengths),to=max(readlengths),n=details_points)[c("x","y")])))
  qualdens <- data.table(type="avg_seq_qual_dens",as.data.table((density(avg_seq_quals,from=min(avg_seq_quals),to=max(avg_seq_quals),n=details_points)[c("x","y")])))
  cycle_qual_sum <- data.table(type="cycle_qual_counts",x=cycle_quals$cycle,y=cycle_quals$qual,z=cycle_quals$Freq)
  if(verbose)
    print("Done!")
  list(overview=overview,plotdata=rbindlist(list(gccontent,cycle_qual_sum,lendens,qualdens),use.names = T,fill = T))
}


result = getOverviewAndDetails(input_file)

fwrite(result$overview,overview_file,sep="\t")
fwrite(result$plotdata,plot_data_file,sep="\t")
