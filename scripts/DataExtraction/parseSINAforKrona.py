#!/usr/bin/env python3

##
## parseSINAforKrona.py - Use the SINA output to generate a Krona text file
##
## Knutt.org/KnuttReads2Bins

# This scripts converts the lineage information in the SINA results
# into a format for ktImportText 

import csv

try:
	snakemake.input
except NameError:
	silvafile = "scripts/parsers/testdata/silvatax.csv"
	kronaoutfile = "test.tsv"
	taxonomyfield = "lca_tax_slv"
else:
	silvafile = snakemake.input["sinacsv"]
	kronaoutfile = snakemake.output["kronatext"]
	taxonomyfield = snakemake.params["taxfield"]



def formatforkrona(classi):
	els = classi.split(";")
	return "\t".join(els)

with open(silvafile,'r') as silva:
	reader = csv.DictReader(silva)
	with open(kronaoutfile, 'w') as outf:
		for row in reader:
				print("1\t"+"\t".join(row[taxonomyfield].split(";")), file=outf)
