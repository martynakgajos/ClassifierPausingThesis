#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import pandas as pd
import sys

def Skew(x,y,seq,w):
    l=[]
    for i in range(len(seq)-w):
        xc=seq[i:i+w].count(x)
        yc=seq[i:i+w].count(y)
        try:
            v=1.*(xc-yc)/(xc+yc)
        except:
            v=0
        l.append(v)
    return np.array(l)

PARSER = argparse.ArgumentParser(description='Script for modeling using RF.')
PARSER.add_argument('-o', dest='outDir', type=str, required=True,
                    help="Output directory.")
PARSER.add_argument('-w', dest='window', type=int, default=100,
                    help="Length of flanking sequences")
PARSER.add_argument('--td', dest='thermodynamics', type=str, default=None,
                    help='''CSV file with thermodynamic features of all 10 mers''')
PARSER.add_argument('--di', dest='dinuc', type=str, default=None,
                    help='''TSV file with -stacking energies for all dinucleotides''')
PARSER.add_argument('-v', dest='vienna', type=str, default=None,
                    help='''Path to RNAVienna package''')
args = PARSER.parse_args()
sys.path.append(args.vienna)
import RNA


#start and end coordinates of potential RNA hairpin
start=11
stop=29
skewness_window=20
skewness_pos1=args.window-10-int(skewness_window/2)
skewness_pos2=args.window+10-int(skewness_window/2)

example_bed=pd.read_csv(args.outDir+'/exampleSequences.csv')

#calculate stacking energy
if args.dinuc:
    dinuc={}
    with open(args.dinuc,'r') as f:
        for ln in f:
            dinuc[ln.strip().split()[0]]=float(ln.strip().split()[1])
    example_bed["stacking_energy"]=example_bed.apply(lambda x: dinuc[x["stacking"]],1)

#calculate thermodynamical features
if args.thermodynamics:
    thermodynamics=pd.read_csv(args.thermodynamics)
    example_bed=example_bed.merge(thermodynamics,how='left',left_on='hyb',right_on='mers')
    example_bed=example_bed.drop(columns='mers')

#calculate GC-content
for col in ['up','hyb','down']:
    example_bed['GC_'+col+'_pct']=example_bed.apply(
            lambda x: (x[col].count('C')+x[col].count('G'))/len(x[col]),1)

#calculate RNA shape
example_bed["RNA_hairpin"]=example_bed.apply(lambda x: 
                                             RNA.fold(str(x["neighbourhood"])[args.window-stop+2:args.window-start+3])[1],1)

#calculate skewness
nuc=['A','T','G','C']

for i in range(len(nuc)):
    for j in range(i+1,len(nuc)):
        a=nuc[i]
        b=nuc[j]
        example_bed[a+b+'_skew']=example_bed.apply(lambda x:
                                                    Skew(a,b,x['neighbourhood'],skewness_window),1)
        example_bed[a+b+'_difference']=example_bed.apply(lambda x:
                                                          x[a+b+'_skew'][skewness_pos1]-x[a+b+'_skew'][skewness_pos2],1)
        example_bed=example_bed.drop(columns=[a+b+'_skew'])

example_bed.to_csv(args.outDir+'/generalFeatures.csv',index=None)

