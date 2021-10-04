#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8) 

def feature_type(x):
    if x.startswith('stacking'):
        return 'stacking energy'
    elif x.startswith('last'):
        return 'last nucleotide added'
    elif x.startswith('next'):
        return 'next nucleotide added'
    elif x in ['dH','T','dS','Gibbs']:
        return 'hybrid features'
    elif x.endswith('HelT')|x.endswith('MGW')|x.endswith('ProT')|x.endswith('Roll')|x.endswith('EP'):
        return 'DNA shape'
    elif x.startswith('MR_')|x.startswith('Z_')|x.startswith('STR_')|x.startswith('DR')|x.startswith('IR')|x.startswith('GQ')|x.startswith('G4'):
        return 'DNA form'
    elif x.startswith('footprint_'):
        return 'footprint motif'
    elif x.startswith('hyb_'):
        return 'hybrid motif'
    elif x.startswith('down_'):
        return 'downstream motif'
    elif x.startswith('up_'):
        return 'upstream motif'
    elif x.startswith('MA'):
        return 'TF'
    elif x.startswith('-1_'):
        return '-1 nucleotide identity'
    elif x.startswith('-11_'):
        return '-11 nucleotide identity'
    elif x.startswith('-10_'):
        return '-10 nucleotide identity'
    elif x.startswith('+1_'):
        return '+1 nucleotide identity'
    elif x.startswith('-2'):
        return '-2 nucleotide identity'
    elif x.startswith('-3'):
        return '-3 nucleotide identity'
    elif x.startswith('RNAmotif'):
        return 'nascent RNA binding'
    else:
        return x

model_best=pickle.load(open(snakemake.input[0],'rb'))
X_test=np.load(snakemake.input[1])
y_test=np.load(snakemake.input[2])
feature_names=list(pd.read_csv(snakemake.input[3],header=None)[0])
k=snakemake.params[0]

pred=model_best.predict_proba(X_test)
scores=pred[:,1]
importances = permutation_importance(model_best, X_test,y_test, 'roc_auc', random_state =42,
                                n_repeats=k,n_jobs=k)

np.save(snakemake.output[0],importances['importances'])


permut_imp=pd.DataFrame(importances['importances'],index=feature_names)
permut_imp['feature']=permut_imp.apply(lambda x: feature_type(x.name),1)
permut_imp=permut_imp.groupby('feature').sum()

temp=permut_imp.std(1)
permut_imp['mean']=permut_imp.mean(1)
permut_imp['std']=temp

permut_imp=permut_imp[((permut_imp['mean']-3*permut_imp['std'])>0)&(permut_imp['mean']>0.005)]
permut_imp=permut_imp.sort_values('mean', ascending=True)

#k=permut_imp.shape[1]-1
n=permut_imp.shape[0]
ax=permut_imp[[i for i in np.arange(1,k)]].T.plot(kind='box',vert=False,showfliers=False, figsize=(6,0.5*n+1),
                                               color=dict(boxes='black', whiskers='black', medians='black', caps='black'))
plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0)
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.grid(False)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 
locs, _ = plt.xticks()
ax.spines['top'].set_bounds(0,locs[-1])
plt.yticks(range(1,n+1,1), permut_imp.index)
plt.title("Permutation feature importances [auc loss]", y=1.05)
plt.ylim(0.3,n+0.7)
plt.xlim(0,locs[-1])
plt.savefig(snakemake.output[1])