#!/usr/bin/env Rscript

##
## CAZyDBReadAnnoforKrona.R - Create functional CAZyDB BLASTX summary for Krona
##
## Knutt.org/KnuttReads2Bins
#
# Read the given files and construct a krona file.

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


infile <- snakemake@input[[1]]
outputfile <- snakemake@output[[1]]
taxcols <- c("superkingdom","phylum","class","order","family","genus","species")

extractCAZy <- function(datacol){
  require(data.table)
  res=data.table(data=unique(datacol))
  res[,dataspl:=strsplit(data," ",fixed=T)]
  ECpattern="[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+\\.[\\d-n]+"
  res[,ECs:=rep(list(),.N)]
  res[grepl(".",dataspl,fixed=T),ECs:=lapply(dataspl,grep,pattern=ECpattern,perl=T,value=T)]
  res[,CazySubClasses:=lapply(dataspl,grep,pattern=ECpattern,invert=T,perl=T,value=T)]
  res[,CazyClasses:=lapply(CazySubClasses,sub,pattern="(\\D{2,3}).+",replacement="\\1")]
  res=res[datacol,on="data"]
  res[,dataspl:=NULL]
  res[,data:=NULL]
  res
}

unnestSingleColUnique<- function(dat,col){
  lookup = dat[[col]] %>% unique %>% tibble(datacol=sapply(.,paste0,collapse=";"),tounnest=.) %>% unnest("tounnest")
  result = dat %>% mutate(datacol=sapply(get(col),paste0,collapse=";")) %>% select(-all_of(col)) %>% full_join(lookup,by="datacol") %>% select(-c(datacol))%>% rename(!!col:="tounnest")
  result
}

hits <- fread(infile)
hits[,(c("ECs","CazySubClasses","CazyClasses")):=extractCAZy(CAZyECs)]
cazydat <- hits[ ,c(taxcols,"CazySubClasses"),with=F]
cazydat <- as.data.table(unnestSingleColUnique(cazydat,"CazySubClasses"))
cazydat[,CazyClasses:=lapply(CazySubClasses,sub,pattern="(\\D{2,3}).+",replacement="\\1")]
cazydat[,CazyClasses:=sapply(CazyClasses,"[[",1)]
cazydat <- cazydat[,.N,by=c("CazyClasses","CazySubClasses","superkingdom", "phylum")]
setcolorder(cazydat, "N")
fwrite(cazydat, outputfile, col.names=F, sep="\t", quote=FALSE)
