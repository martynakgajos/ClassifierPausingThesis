#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams.update({'font.size': 10})
plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8) 

lr=pickle.load(open(snakemake.input[0],'rb'))
rf=pickle.load(open(snakemake.input[1],'rb'))
bt=pickle.load(open(snakemake.input[2],'rb'))

lrdf=pd.DataFrame(lr.cv_results_)[['params', 'mean_test_score', 'mean_train_score','rank_test_score']]
lrdf['Method']='Logistic Regression'
btdf=pd.DataFrame(bt.cv_results_)[['params', 'mean_test_score', 'mean_train_score','rank_test_score']]
btdf['Method']='Gradient Boosting Classifier'
rfdf=pd.DataFrame(rf.cv_results_)[['params', 'mean_test_score', 'mean_train_score','rank_test_score']]
rfdf['Method']='Random Forest'

gs=pd.concat((lrdf,rfdf,btdf))
gs.to_csv(snakemake.output[0])

gs['Difference between train and test AUC']=gs['mean_train_score']-gs['mean_test_score']
gs['Test AUC']=gs['mean_test_score']
fig, ax = plt.subplots(figsize=(4.5,4))
sns.scatterplot( x='Difference between train and test AUC', y='Test AUC', data=gs, hue='Method',
           legend="brief",alpha=0.5,palette=["#25224c","#b4242f","#c7890f"])
plt.legend(loc='lower center')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig(snakemake.output[1])

