#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pylab as plt
import rpy2.robjects.packages as rpackages
DNAshapeR = rpackages.importr('DNAshapeR')
import scipy.stats as st

def calculateMean(x,p):
    l=np.percentile(x,p)
    u=np.percentile(x,100-p)
    m=np.mean(x[(x>l)&(x<u)])
    return [l,u,m]

class Parser(object):
    def __init__(self, outDir):
        self.outDir = outDir
        
args = Parser('/project/owlmayerTemporary/Martyna/Thesis/Modeling/Results/human/GB/Data/')

#Pausing is between 101 and 102 of neighbourhood
shapes=['MGW', 'HelT', 'ProT', 'Roll', 'EP','Stretch', 'Tilt', 'Buckle', 'Shear', 'Opening','Rise', 'Shift', 'Stagger', 'Slide']
    
for shape in shapes:
    res = DNAshapeR.getShape(args.outDir+'sequences.fa',shapeType=shape)

#HelT, Rise, Roll, Shift, Slide, Tilt, Buckle, Opening, ProT, Shear, Stagger, Stretchor, EP
DNAfeatures=shapes
DNAlist={}
for feat in DNAfeatures:
    with open(args.outDir+'sequences.fa.'+feat,'r') as fi:
        rows={1:[],0:[]}
        l=[]
        for ln in fi:
            if ln.startswith('>'):
                if len(l):
                    x=np.array(l[88:105])
                    rows[head].append(x)
                l=[]
                head=int(ln[1])
            else:
                for i in ln.strip().split(','):
                    try:
                        l.append(float(i))
                    except:
                        l.append(0)      
        x=np.array(l[88:105])
        rows[head].append(x)
    DNAlist[feat]=rows

betweenf=['HelT', 'Rise', 'Roll', 'Shift', 'Slide', 'Tilt']

plt.subplots(len(betweenf),1,sharex=True,figsize=(6,8))
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,wspace=0,hspace=0)
for j in range(len(betweenf)):
    feat=betweenf[j]
    if feat in betweenf:
        pos=np.arange(-13.5,3.5)
    else:
        pos=np.arange(-14,3)
    plt.subplot(len(betweenf),1,j+1)
    plt.ylabel(feat)
    mean=np.array(DNAlist[feat][0]).mean(0)
    plt.plot(pos,mean,color='black',alpha=0.2)
    plt.boxplot(np.array(DNAlist[feat][1]),showfliers=False,positions=pos,usermedians=np.array(DNAlist[feat][1]).mean(0))
    plt.xlim([-14.5,2.5])
    if j!=(len(betweenf)-1):
        plt.xticks([],[])
    else:
        plt.xticks([-10,-1,0],['-10','-1','1'])
    print(feat)
    sig=[]
    pvals=[]
    for i in range(mean.shape[0]):
        p=st.ttest_ind(np.array(DNAlist[feat][1])[:,i],np.array(DNAlist[feat][0])[:,i],equal_var=False)[1]
        if p < 10**(-10):
            sig.append(pos[i])
            pvals.append(int(np.log10(p)))
    print(sig)
    #print(pvals)
plt.savefig(args.outDir+'/inbetween.pdf')

inner=['Buckle', 'MGW','Opening', 'ProT', 'Shear', 'Stagger', 'Stretch']
plt.subplots(len(inner),1,sharex=True,figsize=(6,8))
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,wspace=0,hspace=0)
for j in range(len(inner)):
    feat=inner[j]
    if feat not in inner:
        pos=np.arange(-13.5,3.5)
    else:
        pos=np.arange(-14,3)
    plt.subplot(len(inner),1,j+1)
    plt.ylabel(feat)
    mean=np.array(DNAlist[feat][0]).mean(0)
    plt.plot(pos,mean,color='black',alpha=0.2)
    plt.boxplot(np.array(DNAlist[feat][1]),showfliers=False,positions=pos,usermedians=np.array(DNAlist[feat][1]).mean(0))
    plt.xlim([-14.5,2.5])
    if j!=(len(inner)-1):
        plt.xticks([],[])
    else:
        plt.xticks([-10,-1,0],['-10','-1','1'])
    print(feat)
    sig=[]
    pvals=[]
    for i in range(mean.shape[0]):
        p=st.ttest_ind(np.array(DNAlist[feat][1])[:,i],np.array(DNAlist[feat][0])[:,i],equal_var=False)[1]
        if p < 10**(-10):
            sig.append(pos[i])
            pvals.append(int(np.log10(p)))
    print(sig)
    #print(pvals)
plt.savefig(args.outDir+'/inner.pdf')

