#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV


X=np.load(snakemake.input[0])
y=np.load(snakemake.input[1])
jobs=snakemake.params[0]

#Logistic regression
penalty= ['elasticnet']
solver=['saga']
max_iter=[200]
C=[0.001,0.01,0.1,1,10,100,1000]
l1_ratio=[0,0.25,0.5,0.75,1]

lr = LogisticRegression()
lr_grid = {'penalty': penalty,'solver':solver,'C':C,'l1_ratio':l1_ratio,'max_iter':max_iter}
lr_random = GridSearchCV(lr,lr_grid, scoring ='roc_auc', return_train_score=True,
                         cv = 5, verbose=3, n_jobs = jobs)
lr_random.fit(X,y)
pickle.dump(lr_random,open(snakemake.output[0], 'wb'))

#Trees
max_features = ['auto']
n_estimators = [400,800,1500,3000,5000] # Number of trees in random forest
max_depth = [2,4,8,16,32,64,128] # Maximum number of levels in tree
min_samples_leaf = [2,5,10,15] # Minimum number of samples required to split a node
bootstrap = [False,True]  # Method of selecting samples for training each tree                 

#Random Forest                     
rf = RandomForestClassifier()
forest_grid = {'n_estimators': n_estimators,'max_features': max_features, 'max_depth': max_depth,
               'min_samples_leaf': min_samples_leaf, 'bootstrap': bootstrap}
rf_random = GridSearchCV(rf, forest_grid, scoring ='roc_auc',return_train_score=True,
                         cv = 5, verbose=3, n_jobs = jobs)
rf_random.fit(X,y)
pickle.dump(rf_random,open(snakemake.output[1], 'wb'))

#Boosting Trees
learning_rate=[0.001,0.005,0.01,0.05]

bt= GradientBoostingClassifier()
bt_grid = {'n_estimators': n_estimators,'max_features': max_features,'max_depth': max_depth,
           'min_samples_leaf': min_samples_leaf,'learning_rate': learning_rate}
bt_random = GridSearchCV(bt, bt_grid, scoring ='roc_auc',return_train_score=True,
                         cv = 5, verbose=3, n_jobs = jobs)
bt_random.fit(X,y)
pickle.dump(bt_random,open(snakemake.output[2], 'wb'))