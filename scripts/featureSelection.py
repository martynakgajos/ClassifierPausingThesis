#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import seaborn as sns

def plot_collinear(df, correlation_threshold):
    corr_matrix = df.corr(method='spearman')
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k = 1).astype(np.bool))
    to_drop = [column for column in upper.columns if any(upper[column].abs() > correlation_threshold)]

    record_collinear = pd.DataFrame(columns = ['drop_feature', 'corr_feature', 'corr_value'])
    for column in to_drop:
        corr_features = list(upper.index[upper[column].abs() > correlation_threshold])
        corr_values = list(upper[column][upper[column].abs() > correlation_threshold])
        drop_features = [column for _ in range(len(corr_features))]
        temp_df = pd.DataFrame.from_dict({'drop_feature': drop_features,
                                         'corr_feature': corr_features,
                                         'corr_value': corr_values})
        record_collinear = record_collinear.append(temp_df, ignore_index = True)
    corr_matrix_plot = corr_matrix.loc[list(set(record_collinear['corr_feature'])),
                                       list(set(record_collinear['drop_feature']))]
    #corr_matrix_plot[abs(corr_matrix_plot)<correlation_threshold]=0
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    cg=sns.clustermap(corr_matrix_plot.values, cmap=cmap, center=0, linewidths=.01, dendrogram_ratio=0.1)
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)
    cg.ax_heatmap.set_yticks([x + 0.5 for x in list(range(corr_matrix_plot.shape[0]))])
    cg.ax_heatmap.set_yticklabels(list(corr_matrix_plot.index[cg.dendrogram_row.reordered_ind]), rotation=0, size = int(400 / corr_matrix_plot.shape[1]));
    cg.ax_heatmap.set_xticks([x + 0.5 for x in list(range(corr_matrix_plot.shape[1]))])
    cg.ax_heatmap.set_xticklabels(list(corr_matrix_plot.columns[cg.dendrogram_col.reordered_ind]), rotation=90, size = int(400 / corr_matrix_plot.shape[0]));
    

# X=np.load(snakemake.input[0])
# y=list(pd.read_csv(snakemake.input[1],header=None)[0])

X=np.load('/project/owlmayerTemporary/Martyna/Thesis/Modeling/Results/standard_human/PP/Features/normalizedFeatures.npy')
colnames=list(pd.read_csv('/project/owlmayerTemporary/Martyna/Thesis/Modeling/Results/standard_human/PP/Features/featureList.txt',header=None)[0])
features=pd.DataFrame(X, columns=colnames)

plot_collinear(features,0.7)