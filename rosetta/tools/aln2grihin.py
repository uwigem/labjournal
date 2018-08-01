#!/usr/bin/env python
"""
Convert clustal alignment files to grishin for use in rosetta

Assumes first protein is target protein

Args:
    clustal alignment file

@author: Ed
"""

import sys

aln = open(sys.argv[1])
proteins = []

for i, line in enumerate(aln):
    if i == 0 or line == '\n' or line[0] == ' ':
        continue
    words = line.split()
    skip = 0
    for protein in proteins:
        if protein[0] == words[0]:
            protein[1] += words[1]
            skip = 1
            continue
    if not skip:
        proteins.append([words[0], words[1]])

target = proteins[0]

for protein in proteins[1:]:
    grishin = open(target[0] + "_" + protein[0] + ".grishin", "w")
    grishin.write("## %s %s_thread\n#\nscores from program: 0\n0 %s\n0 %s\n" %
                  (target[0], protein[0], target[1], protein[1]))
