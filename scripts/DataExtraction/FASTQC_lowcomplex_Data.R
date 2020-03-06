

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(data.table))


fastq_unmasked_file <- "output/ReadPrep/ClassificationReads/sample/sample_tr_classreads_unsmpld.fastq.gz"
fastq_masked_file <- "output/ReadPrep/ClassificationReads/sample/sample_tr_classreads_unsmpld_masked.fastq.gz"


getNcounts <- function(fastq.path,verbose=F,setsize=1000000,details_points=512){
  fastq.stream <- FastqStreamer(fastq.path,verbose = verbose,n = setsize)
  on.exit(close(fastq.stream),add = T)
  Ncount <- c()
  ids <- c()
  widths <- c()
  repeat {
    entry <- yield(fastq.stream)
    if(length(entry) == 0)
      break
    Ncount <- c(Ncount, alphabetFrequency(sread(entry))[,"N"])
    ids <- c(ids, as.character(id(entry)))
    widths <- c(widths,width(sread(entry)))
  }
  res <- data.table(id=ids,Ncount,widths)
  setkey(res,id)
}

result <- merge(getNcounts(fastq_unmasked_file),getNcounts(fastq_masked_file),suffixes=c(".bef",".af"))
setnafill(result,fill=0, cols=c(2:ncol(result)))
result[, Ncount:=Ncount.af-Ncount.bef]
result[, Nfrac:=(Ncount.bef/widths.af)]

