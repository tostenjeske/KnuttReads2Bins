#!/usr/bin/env python3

##
## combineBandageInfoFiles.py - Fix weird format of Bandage info tsv files
##
## Knutt.org/KnuttReads2Bins

import fileinput

con = [list(map(str.strip, row.replace("%", "").split("\t"))) for row in fileinput.input()]
newcon =  con[0][1:] + con[1][0:1]
print("node_count\tedge_count\tedge_overlap_min\tedge_overlap_max\ttotal_length\ttotal_length_no_overlaps\tdead_ends\tperc_dead_ends\tconnected_components\tlargest_component\ttotal_length_orphans\tN50\tnode_len_min\tnode_len_25th\tnode_len_median\tnode_len_75th\tnode_len_max\tmedian_base_depth\testimated_total_length")
print("\t".join(newcon))