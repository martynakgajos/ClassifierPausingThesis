#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

example_bed=pd.read_csv(snakemake.input[0])

metacol=["chrom","chromStart","chromEnd","name","strand","hyb","up","footprint","down","neighbourhood","stacking"]
labelcol=["score"]
X=example_bed[example_bed.columns[~example_bed.columns.isin(metacol)]]

X=pd.get_dummies(X)
X=X[X.columns[X.nunique()!=1]]
X=X[X.columns[~X.columns.isin(labelcol)]]

y = np.ravel(np.array(example_bed[labelcol]))
feature_list = list(X.columns)

scale=StandardScaler()
X=scale.fit_transform(X)

np.save(snakemake.output[0],X)
np.save(snakemake.output[1],y)
with open(snakemake.output[2], 'w') as f:
    for i in feature_list:
        f.write("%s\n" % i)