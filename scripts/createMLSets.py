#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from sklearn.model_selection import train_test_split

X=np.load(snakemake.input[0])
y=np.load(snakemake.input[1])

set_X, tuning_X, set_y, tuning_y = train_test_split(X, y, test_size = 0.2,
                                                    stratify=y, random_state=42)
train_X, test_X, train_y, test_y = train_test_split(set_X, set_y, test_size = 0.3,
                                                    stratify=set_y, random_state=42)

np.save(snakemake.output[0],tuning_X)
np.save(snakemake.output[1],train_X)
np.save(snakemake.output[2],test_X)

np.save(snakemake.output[3],tuning_y)
np.save(snakemake.output[4],train_y)
np.save(snakemake.output[5],test_y)