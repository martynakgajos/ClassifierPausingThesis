#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import pandas as pd

# Fetch parser arguments
class Parser(object):
    def __init__(self, outDir, bed, region, pval, negativeDistance, maskY):
        self.outDir = outDir
        self.bed = bed
        self.pval = pval
        self.region = region
        self.negativeDistance = negativeDistance
        self.maskY = maskY
        
PARSER = argparse.ArgumentParser(description='Create Negative Set')
PARSER.add_argument('-o', dest='outDir', type=str, required=True,
                    help="Output directory.")
PARSER.add_argument('-b', dest='bed', type=str, required=True,
                    help='''Comma separated list of bed6 files with pausing sites.\n
                    Each bed6 file comes from biological replicate.\n
                    Fields: score - percentile, name - type of pausing site''')
PARSER.add_argument('-p', dest='pval', type=float, default=0.05,
                    help='''P-value cut-off''')
PARSER.add_argument('-r', dest='region', type=str, default=None,
                    help='''Comma separated list of abbrevation of regions
                    included in modeling.\n The name field of bed6 file''')
PARSER.add_argument('--maskY', dest='maskY', type=bool, default=False,
                    help="For female cell lines, mask sites at Y chromosome")
PARSER.add_argument('--negativeDistance', dest='negativeDistance', type=int, default=300,
                    help='''Examples for the negative set will be picked from
                    negativeDistance nucleotides upstream or downstream of a pausing site''')
args = PARSER.parse_args()

negativeDistance_noise=10

#Get list of pausing sites present in all replicates
print('Getting list of pausing sites present in all replicates')
beds=[]
for file_name in args.bed.split(','):
    df=pd.read_csv(file_name,delimiter='\t', header=None,
                   names=['chrom','chromStart','chromEnd','name','score','strand'])
    df['chrom'] = df['chrom'].astype(str)
    df=df[~df["chrom"].str.contains("chrM")]
    if args.maskY:
        df=df[~df["chrom"].str.contains("chrY")]
    if args.region:
        df=df[df["name"].isin(args.region.split(','))]
    df=df[df["score"]<=args.pval]
    df["score"]=1
    beds.append(df)
common_bed=beds[0]
for i in range(1,len(beds)):
    common_bed=common_bed.merge(beds[i],on=list(common_bed.columns))

#Create the negative set
negative_bed=common_bed.copy()
negative_bed["score"]=0
negative_bed["chromStart"]=negative_bed.apply(
        lambda x: x["chromStart"] + 
        int(np.sign(np.random.random()-0.5)*np.random.randint(args.negativeDistance,args.negativeDistance+negativeDistance_noise)),1)
negative_bed["chromEnd"]=negative_bed.apply(lambda x: x["chromStart"]+1,1)

#remove overlapping examples
dix=negative_bed[negative_bed[['chrom','chromStart','strand']].duplicated()].index
common_bed=common_bed.drop(index=dix)
negative_bed=negative_bed.drop(index=dix)

#remove overlapping positive and negative examples
example_bed=pd.concat((common_bed,negative_bed)).sort_values(by=["chrom","chromStart"]).reset_index(drop=True)
dix=example_bed[example_bed[['chrom','chromStart','strand']].duplicated(keep=False)].index
example_bed=example_bed.drop(index=dix)

example_bed.to_csv(args.outDir+'/exampleSet.csv',index=None)
print(sum(example_bed['score']==0)==sum(example_bed['score']==1))