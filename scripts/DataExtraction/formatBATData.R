
##
## formatCATData.R - Convert CAT Data into a more user friendly format
##
## Knutt.org/KnuttReads2Bins
#
options(warn=2)

suppressPackageStartupMessages(library(data.table))

batfile <- snakemake@input[["bat"]]
catfile <- snakemake@input[["cat"]]
binmapfile <- snakemake@input[["binmap"]]
userfile <- snakemake@output[["dat"]]
kronafile <- snakemake@output[["krona"]]

bat <- fread(batfile, sep="\t", fill=T)
setnames(bat, c("# bin"), c("bin"))
taxcols <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
old <- paste0(taxcols, "_")
scores <- paste0(taxcols, "_score")
monocladic <- paste0(taxcols, "_monocladic")
setnames(bat, taxcols, old)
bat[, (old):=lapply(.SD, function(x)sub("not classified", ": ", x, fixed=T)), .SD=old]
bat[, (taxcols):=lapply(.SD, function(x)sub("^(.*): .*$", "\\1", x)), .SD=old]
bat[, (scores):=lapply(.SD, function(x)as.numeric(sub("^.*: (.*)$", "\\1", x))), .SD=old]
bat[, (monocladic):=lapply(.SD, function(x)grepl("^.+\\*$", x)), .SD=taxcols]
bat[, (taxcols):=lapply(.SD, function(x)sub("^(.+)\\*$","\\1", x)), .SD=taxcols]
for (taxi in rev(seq_along(taxcols)[-1])) {
    currentlevel <- taxcols[taxi]
    upperlevel <- taxcols[taxi - 1]
    bat[!is.na(get(currentlevel)) & is.na(get(upperlevel)), (upperlevel) := "(Taxonomy Gap)"]
}
bat[, (old):=NULL]
bat[, classified:=ifelse(classification=="unclassified","Unclassified","Classified")]
bat[, classification:=NULL]
setcolorder(bat, c("bin", "classified"))
fwrite(bat, userfile, sep="\t", quote=F)


cat <- fread(catfile, key="contig")
binmap <- fread(binmapfile, key="contig")
cat <- cat[binmap]
krona <- cat[, .N, by=c("bin" ,"classified", taxcols)]
setcolorder(krona, c("N", "bin", "classified"))
fwrite(krona, kronafile, col.names = F, sep = "\t", quote = F)
system(paste0("sed 's/[[:blank:]]*$//' -i ", kronafile))


