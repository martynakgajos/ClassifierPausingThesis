#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import numpy as np
import pandas as pd
import json
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn import metrics
import matplotlib.pyplot as plt

model=pickle.load(open(snakemake.input[0],'rb'))
X_train=np.load(snakemake.input[1])
y_train=np.load(snakemake.input[2])
X_test=np.load(snakemake.input[3])
y_test=np.load(snakemake.input[4])

model_best=model.best_estimator_

with open(snakemake.output[3], 'w') as f:
    json.dump(model_best.get_params(), f)

model_best.fit(X_train,y_train)
pickle.dump(model_best,open(snakemake.output[0], 'wb'))

pred=model_best.predict_proba(X_test)
scores=pred[:,1]

auc=metrics.average_precision_score(y_test,scores)
fig, ax = plt.subplots(figsize=(4.5,4))
precision, recall, _ = metrics.precision_recall_curve(y_test,scores)
plt.plot(recall, precision, color='r',lw=1,alpha=0.5, label='AUC={0:0.2f}'.format(auc))
plt.plot([0, 1], [0.5, 0.5], linestyle='--', lw=2, color='grey', label='Baseline', alpha=.4)
plt.subplots_adjust(left=0.15, right=0.75, top=0.92, bottom=0.125)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.0])
plt.xlim([0.0, 1.0])
plt.gca().set_aspect('equal', adjustable='box')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.99, box.height])
ax.legend(loc='center left', bbox_to_anchor=(0.95, 0.5),frameon=False)
plt.title('Precision-Recall curve')
plt.savefig(snakemake.output[1])

fig, ax = plt.subplots(figsize=(4.5,4))
auc = metrics.roc_auc_score(y_test,scores)
av_auc=metrics.roc_auc_score(y_test,scores)
fpr, tpr,  threshold = metrics.roc_curve(y_test,scores)
plt.plot(fpr, tpr, color='r',lw=1,alpha=0.5, label='AUC={0:0.2f}'.format(auc))
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey', label='Baseline', alpha=.4)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.subplots_adjust(left=0.15, right=0.75, top=0.92, bottom=0.125)
plt.ylim([0.0, 1.0])
plt.xlim([0.0, 1.0])
plt.gca().set_aspect('equal', adjustable='box')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.99, box.height])
ax.legend(loc='center left', bbox_to_anchor=(0.95, 0.5),frameon=False)
plt.title('ROC curve')
plt.savefig(snakemake.output[2])