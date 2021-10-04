#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn import metrics
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8) 

X_train=np.load(snakemake.input[0])
y_train=np.load(snakemake.input[1])
X_test=np.load(snakemake.input[2])
y_test=np.load(snakemake.input[3])
models={'lr':['Logistic Regression',pickle.load(open(snakemake.input[4],'rb'))],
        'rf':['Random Forest',pickle.load(open(snakemake.input[5],'rb'))],
        'bt':['Gradient Boosting Classifier',pickle.load(open(snakemake.input[6],'rb'))]}

for model in models:
    model_best=models[model][1]
    scores=model_best.predict_proba(X_train)[:,1]   
    auc_train=metrics.roc_auc_score(y_train,scores)
    models[model].append(auc_train)
    scores=model_best.predict_proba(X_test)[:,1]   
    auc_test=metrics.roc_auc_score(y_test,scores)
    models[model].append(auc_test)
    
model_df=pd.DataFrame(models).T

model_df['Method']=model_df[0]
model_df['Test AUC']=model_df[3]
model_df['Difference between train and test AUC']=model_df[2]-model_df[3]

fig, ax = plt.subplots(figsize=(4.5,4))
sns.scatterplot( x='Difference between train and test AUC', y='Test AUC', data=model_df,
                hue='Method', legend="brief",alpha=0.8,palette=["#25224c","#b4242f","#c7890f"])
plt.legend(loc='lower center')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(left=0)
ax.set_ylim(model_df[3].min()-0.05,model_df[3].max()+0.05)
plt.tight_layout()
plt.savefig(snakemake.output[0])
