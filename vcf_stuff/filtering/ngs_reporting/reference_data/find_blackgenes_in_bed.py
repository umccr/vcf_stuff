#!/usr/bin/env python

import sys

gene_names = [l.strip().split('\t')[4] for l in open(sys.argv[1]) if l.strip()]


def read_genes(fpath):
    return [f for f in open(fpath).read().split() if f]


key_genes = read_genes('/Users/vlad/vagrant/NGS_Reporting/az/reference_data/az_key_genes.823.txt')

for fpath in sys.argv[1:]:
    blacklist_genes = read_genes(fpath)
    blacklist_key_genes = [g for g in key_genes if g in blacklist_genes]
    for bl_g in blacklist_key_genes:
        print(fpath + ': ' + bl_g)
