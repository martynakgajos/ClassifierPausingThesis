#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import ast
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8) 

def getCN(x):
    d=ast.literal_eval(x)
    return d['learning_rate']*d['n_estimators']

def getParam(x,p):
    d=ast.literal_eval(x)
    return d[p]

fn='/project/owlmayerTemporary/Martyna/Thesis/Modeling/Results/human/PP/Validation/Validation_scores.csv'
df=pd.read_csv(fn,index_col=[0])

fig, ax = plt.subplots(figsize=(4.6,4))
btdf=df[df['Method']=='Boosting Trees']
btdf['CN']=btdf.apply(lambda x: getCN(x['params']),1)
sns.swarmplot('CN', 'mean_test_score', data=btdf,size=4)
ax.set_ylim(bottom=0.48,top=0.8)
sns.despine()

fig, ax = plt.subplots(figsize=(4.6,4))
lrdf=df[df['Method']=='Logistic Regression']
lrdf['C']=lrdf.apply(lambda x: getParam(x['params'],'C'),1)
sns.swarmplot('C', 'mean_test_score', data=lrdf)
ax.set_ylim(bottom=0.48,top=0.8)
sns.despine()


p='max_depth'
fig, ax = plt.subplots(figsize=(4.6,4))
model_df=df[df['Method']=='Random Forest']
model_df[p]=model_df.apply(lambda x: getParam(x['params'],p),1)
sns.swarmplot('max_depth', 'mean_test_score', data=model_df)
ax.set_ylim(bottom=0.48,top=0.8)
sns.despine()

model_df['Test AUC']=model_df['mean_test_score']
model_df['Difference between train and test AUC']=model_df['mean_train_score']-model_df['mean_test_score']
sns.scatterplot( x='Difference between train and test AUC', y='Test AUC', data=model_df, alpha=0.8,hue=p)
plt.legend(loc='center left')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlim(left=0)
plt.tight_layout()

