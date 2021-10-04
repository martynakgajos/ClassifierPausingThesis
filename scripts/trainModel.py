#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import numpy
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn import metrics

model=pickle.load(open(snakemake.input[0],'rb'))
X_train=np.load(snakemake.input[1])
y_train=np.load(snakemake.input[2])
X_test=np.load(snakemake.input[3])
y_test=np.load(snakemake.input[4])

model_best=model.best_estimator_
#model_best.set_params(presort=True)
model_best.fit(X_train,y_train)
pickle.dump(snakemake.output[0], 'wb'))

pred=model_best.predict_proba(X_test)
scores=pred[:,1]
result = permutation_importance(model_best, X_test,y_test, 'roc_auc', random_state =42,n_repeats=10),
                                n_repeats=snakemake.params[0],n_jobs=snakemake.params[0])
importances.append(result.importances_mean)

auc=metrics.average_precision_score(y_test,scores)
fig, ax = plt.subplots(figsize=(4.5,4))
precision, recall, _ = metrics.precision_recall_curve(y_test,scores)
plt.plot(recall, precision, color='r',lw=1,alpha=0.5, label='AUC={0:0.2f}'.format(auc))
plt.scatter(recall, precision,s=0.1, color='r')
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

fig, ax = plt.subplots(figsize=(4.5,4))
auc = metrics.roc_auc_score(y_test,scores)
av_auc=metrics.roc_auc_score(y_test,scores)
fpr, tpr,  threshold = metrics.roc_curve(y_test,scores)
plt.plot(fpr, tpr, color='r',lw=1,alpha=0.5, label='AUC={0:0.2f}'.format(auc))
plt.scatter(fpr, tpr, color='r',s=0.1)

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