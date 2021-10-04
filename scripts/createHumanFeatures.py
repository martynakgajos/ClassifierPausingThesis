#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import numpy as np
import pandas as pd
import MOODS.parsers
import MOODS.tools
import MOODS.scan
import rpy2.robjects.packages as rpackages
DNAshapeR = rpackages.importr('DNAshapeR')
import pyBigWig

PARSER = argparse.ArgumentParser(description='Script for modeling using RF.')
PARSER.add_argument('-o', dest='outDir', type=str, required=True,
                    help="Output directory.")
PARSER.add_argument('-w', dest='window', type=int, default=100,
                    help="Length of flanking sequences")
PARSER.add_argument('-m', dest='motifs', type=str, default=None,
                    help='''Directory with files with RBP's pfms''')
PARSER.add_argument('--tf', dest='matrices', type=str, default=None,
                    help='''Directory with files with TF's pfms''')
PARSER.add_argument('--shape', dest='shapes', type=str, default=None,
                    help='''Directory with bigWig files with DNA structures prediction''')

args = PARSER.parse_args()

pseudocount = 0.0001
pvalue = 0.01
shape_pos1=args.window-10
shape_pos2=args.window+15

example_bed=pd.read_csv(args.outDir+'/generalFeatures.csv')

#calculate DNAshape
with open(args.outDir+'sequences.fa','w') as f:
    for i in range(example_bed.shape[0]):
        f.write('>'+example_bed.iloc[i]['chrom']+':'+str(example_bed.iloc[i]['chromStart'])+example_bed.iloc[i]['strand']+'\n')
        f.write(str(example_bed.iloc[i]['neighbourhood'])+'\n')        
DNAshapeR.getShape(args.outDir+'sequences.fa')

measure=['min','max','mean','span','mean_diff']
DNAfeatures=['HelT','ProT','Roll','MGW','EP']
DNAlist=[]
for feat in DNAfeatures:
    with open(args.outDir+'sequences.fa.'+feat,'r') as fi:
        rows=[]
        l=[]
        for ln in fi:
            if ln.startswith('>'):
                if len(l):
                    x=np.array(l[shape_pos1:shape_pos2])
                    rows.append(head+[min(x),max(x),np.mean(x),max(x)-min(x),np.mean(abs(np.diff(x)))])
                l=[]
                head=[ln[1:].split(':')[0],int(ln[:-2].split(':')[1]),ln[-2]]
            else:
                for i in ln.strip().split(','):
                    try:
                        l.append(float(i))
                    except:
                        l.append(0)      
        x=np.array(l[shape_pos1:shape_pos2])
        rows.append(head+[min(x),max(x),np.mean(x),max(x)-min(x),np.mean(abs(np.diff(x)))])
    DNAlist.append(pd.DataFrame(rows,columns=['chrom','chromStart','strand']+[m+'_'+feat for m in measure]))
for i in range(len(DNAlist)): 
    example_bed=example_bed.merge(DNAlist[i],on=['chrom','chromStart','strand'],how='right')

#scan for TF motifs:
bg = MOODS.tools.flat_bg(4)
tfs=[f[:-4] for f in os.listdir(args.matrices)]
matrices = [MOODS.parsers.pfm_to_log_odds(args.matrices+f, bg, pseudocount) for f in os.listdir(args.matrices)]
thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]
for col in ['up','down','footprint']:
    temp=[]
    for index, row in example_bed.iterrows():
        results = MOODS.scan.scan_dna(str(row[col]), matrices, bg, thresholds, 7)
        temp.append([1*(len(r)>0) for r in results])
    example_bed[[col+'_'+tf for tf in tfs]] = pd.DataFrame(temp, index=example_bed.index)
pb=tfs.index('PauseButton')
temp=[]
for index, row in example_bed.iterrows():
    r=MOODS.scan.scan_dna(str(row['neighbourhood']), [matrices[pb]], bg, [thresholds[pb]], 7)[0]
    temp.append(1*(len(r)>0))
example_bed['PauseButton'] = pd.DataFrame(temp, index=example_bed.index)

#scan for RNA binding motifs
tfs=[f[:-4] for f in os.listdir(args.motifs)]
matrices = [MOODS.parsers.pfm_to_log_odds(args.motifs+f, bg, pseudocount) for f in os.listdir(args.motifs)]
thresholds = [MOODS.tools.threshold_from_p(m, bg, pvalue) for m in matrices]
for col in ['up']:
    temp=[]
    for index, row in example_bed.iterrows():
        results = MOODS.scan.scan_dna(str(row[col]), matrices, bg, thresholds, 7)
        temp.append([1*(len(r)>0) for r in results])
    example_bed[['RNAmotif_'+tf for tf in tfs]] = pd.DataFrame(temp, index=example_bed.index)

    example_bed['GC_'+col+'_pct']=example_bed.apply(
            lambda x: (x[col].count('C')+x[col].count('G'))/len(x[col]),1)

#DNA forms
formfiles=[i for i in os.listdir(args.shapes)]
for form in set([i.split('.')[0] for i in os.listdir(args.shapes)]):
    if form+'.bw' in formfiles:
        bwp=pyBigWig.open(args.shapes+form+'.bw')
        bwn=bwp
    else:
        bwp=pyBigWig.open(args.shapes+form+'.pos.bw')
        bwn=pyBigWig.open(args.shapes+form+'.neg.bw')
    example_bed[form+'_present_downstream']=example_bed.apply(lambda x:
                                                   np.nan_to_num(bwp.values(x['chrom'],x['chromStart']-100,x['chromStart'])).max() if x['strand']=='+'
                                                   else np.nan_to_num(bwn.values(x['chrom'],x['chromStart'],x['chromStart']+100)).max(),1)
    example_bed[form+'_mean_downstream']=example_bed.apply(lambda x:
                                                   np.nan_to_num(bwp.values(x['chrom'],x['chromStart']-100,x['chromStart'])).sum() if x['strand']=='+'
                                                   else np.nan_to_num(bwn.values(x['chrom'],x['chromStart'],x['chromStart']+100)).sum(),1)
    example_bed[form+'_present_upstream']=example_bed.apply(lambda x:
                                                   np.nan_to_num(bwp.values(x['chrom'],x['chromStart'],x['chromStart']+100)).max() if x['strand']=='+'
                                                   else np.nan_to_num(bwn.values(x['chrom'],x['chromStart']-100,x['chromStart'])).max(),1)
    example_bed[form+'_mean_upstream']=example_bed.apply(lambda x:
                                                   np.nan_to_num(bwp.values(x['chrom'],x['chromStart'],x['chromStart']+100)).sum() if x['strand']=='+'
                                                   else np.nan_to_num(bwn.values(x['chrom'],x['chromStart']-100,x['chromStart'])).sum(),1)

example_bed.to_csv(args.outDir+'/allFeatures.csv',index=None)

