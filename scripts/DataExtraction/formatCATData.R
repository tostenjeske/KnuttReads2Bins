
##
## formatCATData.R - Convert CAT Data into a more user friendly format
##
## Knutt.org/KnuttReads2Bins
#
options(warn=2)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ShortRead))

infile <- snakemake@input[["cat"]]
inseqfile <- snakemake@input[["seq"]]
userfile <- snakemake@output[["dat"]]
kronafile <- snakemake@output[["krona"]]

cat <- fread(infile, sep="\t", fill=T)
setnames(cat, c("# contig"), c("contigid"))
seqs <- data.table(contig=as.character(id(readFasta(inseqfile))))
seqs[, contigid:=tstrsplit(contig, " ", fixed=T)[[1]]]
setkey(seqs, contigid)
setkey(cat, contigid)
cat <- seqs[cat]
cat[, contigid:=NULL]

taxcols <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
old <- paste0(taxcols, "_")
scores <- paste0(taxcols, "_score")
monocladic <- paste0(taxcols, "_monocladic")
setnames(cat, taxcols, old)
cat[, (old):=lapply(.SD, function(x)sub("not classified", ": ", x, fixed=T)), .SD=old]
cat[, (taxcols):=lapply(.SD, function(x)sub("^(.*): .*$", "\\1", x)), .SD=old]
cat[, (scores):=lapply(.SD, function(x)as.numeric(sub("^.*: (.*)$", "\\1", x))), .SD=old]
cat[, (monocladic):=lapply(.SD, function(x)grepl("^.+\\*$", x)), .SD=taxcols]
cat[, (taxcols):=lapply(.SD, function(x)sub("^(.+)\\*$","\\1", x)), .SD=taxcols]
for (taxi in rev(seq_along(taxcols)[-1])) {
    currentlevel <- taxcols[taxi]
    upperlevel <- taxcols[taxi - 1]
    cat[!is.na(get(currentlevel)) & is.na(get(upperlevel)), (upperlevel) := "(Taxonomy Gap)"]
}
cat[, (old):=NULL]
cat[, classified:=ifelse(classification=="unclassified","Unclassified","Classified")]
cat[, classification:=NULL]
setcolorder(cat, c("contig", "classified"))
krona <- cat[, .N, by=c("classified", taxcols)]
setcolorder(krona, c("N", "classified"))
fwrite(krona, kronafile, col.names = F, sep = "\t", quote = F)
system(paste0("sed 's/[[:blank:]]*$//' -i ", kronafile))
fwrite(cat, userfile, sep="\t", quote=F)

