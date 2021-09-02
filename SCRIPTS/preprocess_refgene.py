from __future__ import print_function
import os
import sys

if len(sys.argv) != 2:
    print("Usage: python " + sys.argv[0] + " <refgene file>",
            file=sys.stderr)
    sys.exit()
    
refgene_file = sys.argv[1]

gene_to_refseq = dict()

with open(refgene_file, "rU") as f:
    for line in f:
        items = line.split("\t")
        if items[12] not in gene_to_refseq:
            gene_to_refseq[items[12]] = set()
        gene_to_refseq[items[12]].add(items[1])

for gene in sorted(gene_to_refseq.keys()):
    print(gene + "\t" + ",".join(sorted(gene_to_refseq[gene])))
