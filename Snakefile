##
## Snakefile - The main executable Snakefile for KnuttReads2Bins
##
## Knutt.org/KnuttReads2Bins


# It is just a proxy for the Snakefile_KnuttReads2Bins. This is done to
# make the workflow callable without speciying the Snakefile location
# and also allow the "include"ing of this workflow without conflicts
#  on the "all" rule




include: "Snakefile_0KnuttReads2Bins"
