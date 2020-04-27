#!/usr/bin/env python3

import re
import sys

labelfile = sys.argv[1]
regex = re.compile(sys.argv[2])
joiner = '\n'

with open(labelfile, 'rt') as labels:
    labels = [label for label in labels]

with open(labelfile, 'wt') as out:
    for label in labels:
        out.write(regex.search(label.strip()).group(1))
        out.write('\n')