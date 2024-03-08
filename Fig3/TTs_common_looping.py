import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
from collections import Counter

pt_samples = ['PT1','PT2','PT3','PT4','PT5']
rt_samples = ['RT1','RT2','RT3','RT4','RT5']
# RTvsNT - (PTvsNT)
datdir = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis'
gain_df = pd.read_table('{}/individual_specific_loops/IndSp_Gained.txt'.format(datdir),sep='\t',header=0,index_col=0)
loss_df = pd.read_table('{}/individual_specific_loops/IndSp_Lost.txt'.format(datdir),sep='\t',header=0,index_col=0)
signal_df = pd.read_csv('{}/GC_loops/all_gene_medratio_mat.csv'.format(datdir),index_col=0)
# find PT's common genes with criterion at least 2 samples have them
common_cut = 2
gain_df['PTs'] = (gain_df[pt_samples].sum(axis=1) >= common_cut).astype(int)
loss_df['PTs'] = (loss_df[pt_samples].sum(axis=1) >= common_cut).astype(int)
summary_dict = {'CommonType':[],'Gained':[],'Lost':[]}
for i in range(1,11):
    summary_dict['CommonType'].append('>={}TTs'.format(i))
    summary_dict['Gained'].append(len(gain_df[gain_df.sum(axis=1)==i]))
    summary_dict['Lost'].append(len(loss_df[loss_df.sum(axis=1)==i]))

summary_df = pd.DataFrame(summary_dict)
summary_df.to_excel('{}/individual_specific_loops/TTsIndvsNTs_commons_summary.xlsx'.format(datdir),index=False)