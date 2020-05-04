#!/usr/bin/env python3

import sourmash
import os
import subprocess
import tempfile
import sys

sigfile = snakemake.input["sig"]
dbfile = snakemake.input["db"]
params = snakemake.params["params"]
outfile = snakemake.output[0]
header = snakemake.params["header"]
logfile = snakemake.log[0]


with tempfile.NamedTemporaryFile('r+t') as tmpsig, tempfile.NamedTemporaryFile('r+t') as tmpcsv:
    with open(sigfile, 'rt') as sigfile, open(outfile, 'wt') as outfile, open(logfile, 'wt') as logfile:
     sigiter = sourmash.signature.load_signatures(sigfile)
     outfile.write(header)
     outfile.write('\n')
     for sig in sigiter:
            binfile = sig.name()
            binname = os.path.basename(binfile)
            sourmash.signature.save_signatures([sig],tmpsig)
            cmd = params + ["-o", tmpcsv.name, tmpsig.name, dbfile]
            p = subprocess.Popen(cmd,  stdout=logfile, stderr=logfile)
            p.wait()
            firstline = True
            for line in tmpcsv:
                if firstline:
                    firstline = False
                    continue
                outfile.write(binname)
                outfile.write("\t") 
                outfile.write(line.replace(",","\t"))
            tmpcsv.truncate()
            tmpcsv.seek(0)
            tmpsig.truncate()
            tmpsig.seek(0)

