#!/usr/bin/env python3

##
## SILVAforKrona.py - Use the SILVA db fasta to generate a Krona text file
##
## Knutt.org/KnuttReads2Bins


import csv

try:
	snakemake.input
except NameError:
	silvafile = "reference_data/SSU/SILVA_132_SSURef_Nr99_tax_silva.fasta"
	kronaoutfile = "test.tsv"
else:
	silvafile = snakemake.input[0]
	kronaoutfile = snakemake.output[0]



def formatforkrona(classi):
	els = classi.split(";")
	return "\t".join(els)

with open(silvafile,'r') as silva:
	with open(kronaoutfile, 'w') as outf:
		for row in silva:
			row = row.strip()
			if row.startswith(">"):
				row = "".join(row.split(" ")[1:])
				print("1\t"+"\t".join(row.split(";")),file=outf)
