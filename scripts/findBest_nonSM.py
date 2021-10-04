#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
import matplotlib.pyplot as plt
import seaborn as sns

dire='/project/owlmayerTemporary/Martyna/Thesis/Modeling/human/PP/Models/'
lr=pickle.load(open(dire+'lr_random.pickle','rb'))
rf=pickle.load(open(dire+'rf_random.pickle','rb'))
bt=pickle.load(open(dire+'bt_random.pickle','rb'))

lrdf=pd.DataFrame(lr.cv_results_)[['params', 'mean_test_score', 'mean_train_score','rank_test_score']]
lrdf['method']='Logistic Regression'
btdf=pd.DataFrame(bt.cv_results_)[['params', 'mean_test_score', 'mean_train_score','rank_test_score']]
btdf['method']='Boosting Trees'
rfdf=pd.DataFrame(rf.cv_results_)[['params', 'mean_test_score', 'mean_train_score','rank_test_score']]
rfdf['method']='Random Forest'

gs=pd.concat((lrdf,rfdf,btdf))

gs['Difference between train and test AUC']=gs['mean_train_score']-gs['mean_test_score']
gs['Test AUC']=gs['mean_test_score']
fig, ax = plt.subplots(figsize=(4.5,4))
sns.scatterplot( x='Difference between train and test AUC', y='Test AUC', data=gs, hue='method',
           legend="brief",alpha=0.5,palette=["#25224c","#b4242f","#c7890f"])
plt.legend(loc='lower center')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.tight_layout()

