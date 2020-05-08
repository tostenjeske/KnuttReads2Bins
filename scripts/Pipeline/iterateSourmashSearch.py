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



with open(sigfile, 'rt') as sigfile, open(outfile, 'wt') as outfile, open(logfile, 'wt') as logfile:
    sigiter = sourmash.signature.load_signatures(sigfile)
    outfile.write(header)
    outfile.write('\n')
    for sig in sigiter:
        binfile = sig.name()
        binname = os.path.basename(binfile)
        with tempfile.NamedTemporaryFile('wt') as tmpsig, tempfile.NamedTemporaryFile('rt') as tmpcsv:
            sourmash.signature.save_signatures([sig],tmpsig)
            tmpsig.flush()
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

