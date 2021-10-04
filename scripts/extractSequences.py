#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import pandas as pd
from Bio import SeqIO

def ReverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','t':'A','g':'C','c':'G','a':'T','N':'N'}
    s=s[::-1]
    seq=''
    for i in s:
        seq+=complement[i]
    return seq
       
PARSER = argparse.ArgumentParser(description='Extract Sequences of ROI')
PARSER.add_argument('-o', dest='outDir', type=str, required=True,
                    help="Output directory.")
PARSER.add_argument('-f', dest='fasta', type=str, required=True,
                    help="Fasta file containing reference genome")
PARSER.add_argument('-w', dest='window', type=int, default=100,
                    help="Length of flanking sequences")
PARSER.add_argument('--hybrid', dest='hybrid', type=int, default=10,
                    help="Length of the DNA-RNA hybrid")
args = PARSER.parse_args()

half_footprint=25

example_bed=pd.read_csv(args.outDir+'/exampleSet.csv')

genome=SeqIO.index(args.fasta,"fasta")
chrm=''
hybrid_seqs, downstream_seqs, upstream_seqs, bacterial, stacking, all_seqs, footprint = [], [], [], [], [], [], []
for loc in range(example_bed.shape[0]):
    row=example_bed.iloc[loc]
    if row["chrom"]!=chrm:
        chrm=str(row["chrom"])
        chrom_sequence=genome[chrm].seq[:]
    if row["strand"]=='+':
        upstream_seqs.append(chrom_sequence[row['chromStart']+1-args.hybrid-args.window:row['chromStart']+1-args.hybrid].upper())
        bacterial.append(chrom_sequence[row['chromStart']-args.hybrid:row['chromStart']+2].upper())
        downstream_seqs.append(chrom_sequence[row['chromStart']+1:row['chromStart']+1+args.window].upper()) 
        stacking.append(chrom_sequence[row['chromStart']:row['chromStart']+2].upper())
        hybrid_seqs.append(chrom_sequence[row['chromStart']+1-args.hybrid:row['chromStart']+1].upper())
        all_seqs.append(chrom_sequence[row['chromStart']-args.window-1:row['chromStart']+args.window+1].upper())
        footprint.append(chrom_sequence[row['chromStart']-half_footprint-1:row['chromStart']+half_footprint+1].upper())
    else:
        hybrid_seqs.append(chrom_sequence[row['chromStart']:args.hybrid+row['chromStart']].reverse_complement().upper())
        all_seqs.append(chrom_sequence[row['chromStart']-args.window:row['chromStart']+args.window+2].reverse_complement().upper())
        footprint.append(chrom_sequence[row['chromStart']-half_footprint:row['chromStart']+half_footprint+2].reverse_complement().upper())
        upstream_seqs.append(chrom_sequence[row['chromStart']+args.hybrid:args.window+row['chromStart']+args.hybrid].reverse_complement().upper())
        downstream_seqs.append(chrom_sequence[row['chromStart']-args.window:row['chromStart']].reverse_complement().upper())
        bacterial.append(chrom_sequence[row['chromStart']-1:args.hybrid+row['chromStart']+1].reverse_complement().upper())
        stacking.append(chrom_sequence[row['chromStart']-1:row['chromStart']+1].reverse_complement().upper())

example_bed["neighbourhood"]=all_seqs
example_bed["footprint"]=footprint
example_bed["up"]=upstream_seqs
example_bed["down"]=downstream_seqs
example_bed["hyb"]=hybrid_seqs
example_bed["stacking"]=stacking
example_bed["+1_template"]=[i[-1] for i in bacterial]
example_bed["-1_template"]=[i[-2] for i in bacterial]
example_bed["-2_template"]=[i[-3] for i in bacterial]
example_bed["-3_template"]=[i[-4] for i in bacterial]
example_bed["-10_template"]=[i[-1-args.hybrid] for i in bacterial]
example_bed["-11_template"]=[i[-2-args.hybrid] for i in bacterial]
example_bed.to_csv(args.outDir+'/exampleSequences.csv',index=None)