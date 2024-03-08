import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly
import pybedtools as pybt
import plotly.graph_objects as go
import itertools
from scipy.stats import pearsonr

def bin2compartment(df,col):
    # find idx of flipped row first
    vals = df[col].values
    abs_sign_diff = np.abs(np.diff(np.sign(vals)))
    # idx of first row where the change is
    change_idx = np.flatnonzero(abs_sign_diff == 2)
    # +1 to get idx of second rows in the sign change too
    change_idx = np.stack((change_idx, change_idx + 1), axis=1)

    # now we have the locations where sign changes occur. We just need to extract
    # the `value` values at those locations to determine which of the two possibilities
    # to choose for each sign change (whichever has `value` closer to 0)
    min_idx = np.abs(vals[change_idx]).argmin(1)
    split_idx = change_idx[range(len(change_idx)), min_idx]
    split_df = np.split(df, split_idx)
    final_df_dict = {'chr':[],'start':[],'end':[],'compartment':[],'size_cat':[],'pcOri.avg':[]}
    for df in split_df:
        for chr, chr_df in df.groupby('chr'):
            final_df_dict['chr'].append(chr)
            final_df_dict['start'].append(chr_df['start'].min())
            final_df_dict['end'].append(chr_df['end'].max())
            size_mb = chr_df['end'].max()/1e6 - chr_df['start'].min()/1e6
            final_df_dict['pcOri.avg'].append(np.mean(chr_df[col]))
            if chr_df[col].mean()>=0:
                final_df_dict['compartment'].append('A')
            else:
                final_df_dict['compartment'].append('B')
            if size_mb <= 1:
                final_df_dict['size_cat'].append('<=1Mb')
            elif size_mb >= 5:
                final_df_dict['size_cat'].append('>=5Mb')
            else:
                final_df_dict['size_cat'].append('{}-{}Mb'.format(int(np.ceil(size_mb))-1,int(np.ceil(size_mb))))
    final_df = pd.DataFrame(final_df_dict)
    return final_df

def bin2compartment_avg(df,col):
    # find idx of flipped row first
    vals = df[col].values
    abs_sign_diff = np.abs(np.diff(np.sign(vals)))
    # idx of first row where the change is
    change_idx = np.flatnonzero(abs_sign_diff == 2)
    # +1 to get idx of second rows in the sign change too
    change_idx = np.stack((change_idx, change_idx + 1), axis=1)

    # now we have the locations where sign changes occur. We just need to extract
    # the `value` values at those locations to determine which of the two possibilities
    # to choose for each sign change (whichever has `value` closer to 0)
    min_idx = np.abs(vals[change_idx]).argmin(1)
    split_idx = change_idx[range(len(change_idx)), min_idx]
    split_df = np.split(df, split_idx)
    final_df_dict = {'chr':[],'start':[],'end':[],'compartment':[],'size_cat':[],
                     'pcQnm.avg':[],'padj.min':[],'maha.max':[]}
    for df in split_df:
        for chr, chr_df in df.groupby('chr'):
            final_df_dict['chr'].append(chr)
            final_df_dict['start'].append(chr_df['start'].min())
            final_df_dict['end'].append(chr_df['end'].max())
            size_mb = chr_df['end'].max()/1e6 - chr_df['start'].min()/1e6
            final_df_dict['pcQnm.avg'].append(np.mean(chr_df[col]))
            final_df_dict['maha.max'].append(np.max(chr_df['sample_maha']))
            final_df_dict['padj.min'].append(np.max(chr_df['padj']))
            if chr_df[col].mean()>=0:
                final_df_dict['compartment'].append('A')
            else:
                final_df_dict['compartment'].append('B')
            if size_mb <= 1:
                final_df_dict['size_cat'].append('<=1Mb')
            elif size_mb >= 5:
                final_df_dict['size_cat'].append('>=5Mb')
            else:
                final_df_dict['size_cat'].append('{}-{}Mb'.format(int(np.ceil(size_mb))-1,int(np.ceil(size_mb))))
    final_df = pd.DataFrame(final_df_dict)
    return final_df

def plot_sankeyplot(sankey_dict, node_x, node_y, outname):
    layout = go.Layout(
        autosize=False,
        width=500,
        height=500
    )
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=['<b>' + label + '</b>' for label in sankey_dict['label']],
            color=sankey_dict['color_node'],
            x=node_x,
            y=node_y
        ),
        link=dict(
            source=sankey_dict['source'],  # indices correspond to labels, eg A1, A2, A1, B1, ...
            target=sankey_dict['target'],
            value=sankey_dict['value'],
            color=sankey_dict['color_link']
        ),
        textfont=dict(color="rgba(0,0,0,0)", size=1)
    )], layout=layout)
    fig.update_layout(title_text="{}-{} CompartmentA/B Transition".format(sankey_dict['label'][0].split('_')[0],
                                                                          sankey_dict['label'][2].split('_')[0]),
                      font_size=10)
    fig.write_image(outname)

datdir = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/dcHic/100k'
nt_tt_pcqnm = 'NT_TT_differential.intra_sample_combined.pcQnm.bedGraph'
pt_rt_pcqnm = 'PT_RT_differential.intra_sample_combined.pcQnm.bedGraph'
nt_pt_rt_pcori = 'NT_PT_RT_combined.pcOri.srt.bedGraph'

fdr_cut = 0.1
res = 100000
diff_color_dict = {'Diff':'red','noDiff':'black'}
samples = ['NT1', 'NT2', 'PT1', 'PT2', 'PT3', 'PT4', 'PT5', 'RT1', 'RT2', 'RT3', 'RT4', 'RT5']
size_cate_color = {'<=1Mb':'navy',
                   '1-2Mb':'steelblue',
                   '2-3Mb':'lightseagreen',
                   '3-4Mb':'forestgreen',
                   '4-5Mb':'darkorange',
                   '>=5Mb':'brown'
                   }
nt_pt_rt_pcori_df = pd.read_table('{}/{}'.format(datdir, nt_pt_rt_pcori))
nt_tt_pcqnm_df = pd.read_table('{}/{}'.format(datdir, nt_tt_pcqnm))
pt_rt_pcqnm_df = pd.read_table('{}/{}'.format(datdir, pt_rt_pcqnm))

# count number of compartment A/B in size category
compartment_dict = {}
compartment_stats = {'Compartment A':{cat:[] for cat in size_cate_color},
                     'Compartment B':{cat:[] for cat in size_cate_color}
                     }
for sample in samples:
    print("Process {}".format(sample))
    current_pcori = nt_pt_rt_pcori_df[['chr','start','end','{}_{}'.format(sample,res)]]
    current_comp = bin2compartment(current_pcori,'{}_{}'.format(sample,res))
    compartment_dict[sample] = current_comp
    for cate in size_cate_color:
        compartment_stats['Compartment A'][cate].append(len(current_comp.loc[(current_comp['compartment']=='A')&
                                                                        (current_comp['size_cat']==cate)]))
        compartment_stats['Compartment B'][cate].append(len(current_comp.loc[(current_comp['compartment'] == 'B') &
                                                                        (current_comp['size_cat'] == cate)]))

# save process compartment result
writer = pd.ExcelWriter('{}/compartments_all_samples.xlsx'.format(datdir), engine='xlsxwriter')
for sample in samples:
    compartment_dict[sample].to_excel(writer, sheet_name=sample, index=False)
writer.save()

# Fig.1B
fig, axes = plt.subplots(2, 1, sharex=True)
for idx, compartment in enumerate(compartment_stats.keys()):
    bottom_height = np.zeros(len(samples))
    for cate in size_cate_color:
        axes[idx].bar(samples,
                      compartment_stats[compartment][cate],
                      bottom = bottom_height,
                      color = size_cate_color[cate],
                      label = cate)
        bottom_height += np.array(compartment_stats[compartment][cate])
        axes[idx].tick_params(axis='y',labelsize=12)
plt.legend(bbox_to_anchor=(1.02, 1.5))
plt.savefig('{}/results/CompartmentNum_Barplot.png'.format(datdir),dpi=300, bbox_inches='tight')
plt.close()

# Density plot for NTvsTT, PTvsRT
nt_tt_pcqnm_df['mark'] = np.where(nt_tt_pcqnm_df['padj'] <= fdr_cut,'Diff','noDiff')
nt_tt_pcqnm_df['-log10padj'] = -np.log10(nt_tt_pcqnm_df['padj'])
pt_rt_pcqnm_df['mark'] = np.where(pt_rt_pcqnm_df['padj'] <= fdr_cut,'Diff','noDiff')
pt_rt_pcqnm_df['-log10padj'] = -np.log10(pt_rt_pcqnm_df['padj'])
# Fig.1C
fig, axes = plt.subplots(2,1,figsize=(6,12))
sns.scatterplot(data=nt_tt_pcqnm_df, x='NT', y='TT', hue='mark', linewidth=0, s=2.5,
                palette=diff_color_dict, ax=axes[0], legend=False)
xpoints = ypoints = axes[0].get_xlim()
axes[0].plot(xpoints, ypoints, linestyle='--', color='k', lw=1, scalex=False, scaley=False)

sns.scatterplot(data=pt_rt_pcqnm_df, x='PT', y='RT', hue='mark', linewidth=0, s=2.5, palette=diff_color_dict,ax=axes[1])
xpoints = ypoints = axes[1].get_xlim()
axes[1].plot(xpoints, ypoints, linestyle='--', color='k', lw=1, scalex=False, scaley=False)
# axes[1].legend(bbox_to_anchor=(1.25, 0.55))
plt.savefig('{}/results/pcQnm_Scatter_082623.png'.format(datdir),dpi=300,bbox_inches='tight')
plt.close()

hg19_size = 3099734149 - 57264655 - 4485509 # hg19 all scaffold - chrY - chrUn accroding https://www.ncbi.nlm.nih.gov/grc/human/data
nt_tt_pcqnm_df_db = nt_tt_pcqnm_df[nt_tt_pcqnm_df['padj']<=fdr_cut]
pt_rt_pcqnm_df_db = pt_rt_pcqnm_df[pt_rt_pcqnm_df['padj']<=fdr_cut]
# percent of all genome
print("nt_tt db {}%".format(100*len(nt_tt_pcqnm_df_db)*res/hg19_size))
print("pt_rt db {}%".format(100*len(pt_rt_pcqnm_df_db)*res/hg19_size))

nt_tt_pcqnm_df_db['NT_comp'] = np.where(nt_tt_pcqnm_df_db['NT']>0,'A','B')
nt_tt_pcqnm_df_db['TT_comp'] = np.where(nt_tt_pcqnm_df_db['TT']>0,'A','B')
nt_tt_pcqnm_df_db['SwitchType'] = nt_tt_pcqnm_df_db['NT_comp'] + '-' + nt_tt_pcqnm_df_db['TT_comp']

nt_tt_sankey_dict = {
    'source': [0, 0, 1, 1],
    'target': [2, 3, 2, 3],
    'value': [len(nt_tt_pcqnm_df_db[nt_tt_pcqnm_df_db['SwitchType']=='A-A']),
              len(nt_tt_pcqnm_df_db[nt_tt_pcqnm_df_db['SwitchType']=='A-B']),
              len(nt_tt_pcqnm_df_db[nt_tt_pcqnm_df_db['SwitchType']=='B-A']),
              len(nt_tt_pcqnm_df_db[nt_tt_pcqnm_df_db['SwitchType']=='B-B'])],
    'label': ['NT_A', 'NT_B', 'TT_A', 'TT_B'],
    'color_link': ['darkorange', 'darkviolet', 'darkorange', 'darkviolet'],
    'color_node': ['darkorange', 'darkviolet', 'darkorange', 'darkviolet']
}

plot_sankeyplot(nt_tt_sankey_dict,[0.1, 0.1, 0.9, 0.9],[0.2, 0.8, 0.135, 0.74],'{}/NT_TT_CompTrans_080623.png'.format(datdir))

pt_rt_pcqnm_df_db['PT_comp'] = np.where(pt_rt_pcqnm_df_db['PT']>0,'A','B')
pt_rt_pcqnm_df_db['RT_comp'] = np.where(pt_rt_pcqnm_df_db['RT']>0,'A','B')
pt_rt_pcqnm_df_db['SwitchType'] = pt_rt_pcqnm_df_db['PT_comp'] + '-' + pt_rt_pcqnm_df_db['RT_comp']

pt_rt_sankey_dict = {
    'source': [0, 0, 1, 1],
    'target': [2, 3, 2, 3],
    'value': [len(pt_rt_pcqnm_df_db[pt_rt_pcqnm_df_db['SwitchType']=='A-A']),
              len(pt_rt_pcqnm_df_db[pt_rt_pcqnm_df_db['SwitchType']=='A-B']),
              len(pt_rt_pcqnm_df_db[pt_rt_pcqnm_df_db['SwitchType']=='B-A']),
              len(pt_rt_pcqnm_df_db[pt_rt_pcqnm_df_db['SwitchType']=='B-B'])],
    'label': ['PT_A', 'PT_B', 'RT_A', 'RT_B'],
    'color_link': ['darkorange', 'darkviolet', 'darkorange', 'darkviolet'],
    'color_node': ['darkorange', 'darkviolet', 'darkorange', 'darkviolet']
}

plot_sankeyplot(pt_rt_sankey_dict,[0.1, 0.1, 0.9, 0.9],[0.2, 0.8, 0.4, 0.9999],'{}/PT_RT_CompTrans_080623.png'.format(datdir))

# boxplot show compartment scores
nt_tt_pcqnm_boxplot = pd.melt(nt_tt_pcqnm_df_db[['SwitchType','NT','TT']],
                              id_vars=['SwitchType'], value_vars=['NT', 'TT'])
nt_tt_pcqnm_boxplot['variable2'] = np.where(nt_tt_pcqnm_boxplot['value']>0,'A','B')
nt_tt_pcqnm_boxplot['types'] = nt_tt_pcqnm_boxplot['variable'] + '_' + nt_tt_pcqnm_boxplot['variable2']

plt.figure()
sns.boxplot(nt_tt_pcqnm_boxplot,x='types',y='value',order=['NT_A','TT_B','NT_B','NT_A','NT_A','TT_A','NT_B','TT_B'],
            palette={'NT_A':'darkorange','NT_B':'darkviolet','TT_A':'darkorange','TT_B':'darkviolet'})
plt.savefig('{}/NT_TT_compartment_score.png'.format(datdir),dpi = 400)
plt.close()

pt_rt_pcqnm_boxplot = pd.melt(pt_rt_pcqnm_df_db[['SwitchType','PT','RT']],
                              id_vars=['SwitchType'], value_vars=['PT', 'RT'])
pt_rt_pcqnm_boxplot['variable2'] = np.where(pt_rt_pcqnm_boxplot['value']>0,'A','B')
pt_rt_pcqnm_boxplot['types'] = pt_rt_pcqnm_boxplot['variable'] + '_' + pt_rt_pcqnm_boxplot['variable2']

plt.figure()
sns.boxplot(pt_rt_pcqnm_boxplot,x='types',y='value',order=['PT_A','RT_B','PT_B','RT_A','PT_A','RT_A','PT_B','RT_B'],
            palette={'PT_A':'darkorange','PT_B':'darkviolet','RT_A':'darkorange','RT_B':'darkviolet'})
plt.savefig('{}/PT_RT_compartment_score.png'.format(datdir),dpi = 400)
plt.close()

# save
writer2 = pd.ExcelWriter('{}/compartments_dcHiC.xlsx'.format(datdir), engine='xlsxwriter')
nt_tt_pcqnm_df.to_excel(writer2, sheet_name='NT_TT_dcHiC', index=False)
nt_tt_pcqnm_df_db.to_excel(writer2, sheet_name='NT_TT_diff', index=False)
pt_rt_pcqnm_df.to_excel(writer2, sheet_name='PT_RT_dcHiC', index=False)
pt_rt_pcqnm_df_db.to_excel(writer2, sheet_name='PT_RT_diff', index=False)
writer2.save()

# Fig.1E
# find pearson corr among samples
samples_pc_corr_mat = pd.DataFrame(0,columns=samples,index=samples)
for subset in itertools.combinations(samples, 2):
    sample1 = subset[0]
    sample2 = subset[1]
    pearsonr_result = pearsonr(nt_pt_rt_pcori_df['{}_{}'.format(sample1, res)],
                               nt_pt_rt_pcori_df['{}_{}'.format(sample2, res)])
    samples_pc_corr_mat.loc[sample1,sample2] = pearsonr_result[0]
    samples_pc_corr_mat.loc[sample2,sample1] = pearsonr_result[0]

def get_lower_tri_heatmap(df, show_text=False, output="cooc_matrix.png"):
    mask = np.triu(np.ones_like(df))
    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns_plot = sns.heatmap(df, cmap=cmap, vmin=0.5 ,vmax=1, center=0.75, annot=df, mask=mask,
                           square=True, linewidths=.5, cbar_kws={"shrink": .5},annot_kws={"fontsize":16,"color":"black"})
    for text in ax.texts:
        text.set_visible(show_text)
    plt.savefig(output)
    plt.close()

def get_tri_heatmap(df, loc='lower', fontsize=16, show_text=False, output="cooc_matrix.png"):
    if loc == 'lower':
        df_tri = df.where(np.tril(np.ones(df.shape)).astype(bool))
    elif loc == 'upper':
        df_tri = df.where(np.triu(np.ones(df.shape)).astype(bool))
    else:
        print("Wrong loc {} params, use lower".format(loc))
        df_tri = df.where(np.tril(np.ones(df.shape)).astype(bool))
    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns_plot = sns.heatmap(df_tri, cmap=cmap, vmin=0.5, vmax=1, center=0.75, annot=df,
                           square=True, linewidths=.5, cbar_kws={"shrink": .5},
                           annot_kws={"fontsize": fontsize, "color": "black"})
    if loc == 'upper':
        ax.xaxis.tick_top()  # x axis on top
        ax.xaxis.set_label_position('top')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')

    for text in ax.texts:
        text.set_visible(show_text)
    plt.savefig(output)
    plt.close()

get_tri_heatmap(samples_pc_corr_mat, loc='lower', show_text=False, output='{}/samples_pearsonr_heatmap.png'.format(datdir))
samples_pc_corr_mat = samples_pc_corr_mat.replace(0,1)
tissueType = ['NT','PT','RT']
sample_pc_merge_corr_mat = pd.DataFrame(0, columns=tissueType, index=tissueType)
sample_pc_merge_corr_mat.loc['NT','NT'] = samples_pc_corr_mat.loc['NT1','NT2']
sample_pc_merge_corr_mat.loc['NT','PT'] = samples_pc_corr_mat.loc[['NT1','NT2'],['PT{}'.format(i) for i in range(1,6)]].mean().mean()
sample_pc_merge_corr_mat.loc['NT','RT'] = samples_pc_corr_mat.loc[['NT1','NT2'],['RT{}'.format(i) for i in range(1,6)]].mean().mean()
sample_pc_merge_corr_mat.loc['PT','PT'] = samples_pc_corr_mat.loc[['PT{}'.format(i) for i in range(1,6)],
                                                                  ['PT{}'.format(i) for i in range(1,6)]].mean().mean()
sample_pc_merge_corr_mat.loc['PT','RT'] = samples_pc_corr_mat.loc[['PT{}'.format(i) for i in range(1,6)],
                                                                  ['RT{}'.format(i) for i in range(1,6)]].mean().mean()
sample_pc_merge_corr_mat.loc['RT','RT'] = samples_pc_corr_mat.loc[['RT{}'.format(i) for i in range(1,6)],
                                                                  ['RT{}'.format(i) for i in range(1,6)]].mean().mean()
get_tri_heatmap(sample_pc_merge_corr_mat, loc='upper', fontsize=30, show_text=True, output='{}/samples_mean_pearsonr_heatmap.png'.format(datdir))

# write sample pcori into igv bedGraph
# header: track name=" TT  PC" description="BedGraph format" visibility=full color= 234,100,0  altColor= 122,16,180
# priority=20 plotType="points"
for sample in samples:
    with open('{}/{}.pcori.bedGraph'.format(datdir,sample), 'w') as file:
        file.write('track name=\"{}\" description=\"BedGraph format\" '
                   'visibility=full color=234,100,0  altColor=122,16,180 priority=20 plotType=\"points\"\n'.format(sample))
        nt_pt_rt_pcori_df[['chr','start','end','{}_100000'.format(sample)]].to_csv(file, header=False, index=False, sep='\t',mode='a')


# def find_switch(compartment_df1,compartment_df2,col1,col2):
#     comp_bt1 = pybt.BedTool.from_dataframe(compartment_df1[['chr', 'start', 'end', 'compartment']])
#     comp_bt2 = pybt.BedTool.from_dataframe(compartment_df2[['chr', 'start', 'end', 'compartment']])
#     inter_df = pybt.BedTool.intersect(comp_bt1, comp_bt2, wa=True, wb=True).to_dataframe(disable_auto_names=True,
#                                                                                          names=['chr1',
#                                                                                                 'start1',
#                                                                                                 'end1',
#                                                                                                 '{}_comp'.format(col1),
#                                                                                                 'chr2',
#                                                                                                 'start2',
#                                                                                                 'end2',
#                                                                                                 '{}_comp'.format(col2)])
#     # remove short overlap fragment, this step might lose some small compartment in comp_bt2
#     inter_df['interSizeRatio'] = inter_df.apply(lambda x:
#                                                 (min(x['end1'], x['end2']) - max(x['start1'], x['start2'])) /
#                                                 (x['end1'] - x['start1']),
#                                                 axis=1)
#     inter_df['SwitchType'] = inter_df['{}_comp'.format(col1)] + '-' + inter_df['{}_comp'.format(col2)]
#     return inter_df
#
# # profiling differential information. As we need to compare, we use pcQnm
# compartment_NT_pcqnm = bin2compartment_avg(nt_tt_pcqnm_df[['chr','start','end','NT','padj','sample_maha']],'NT')
# compartment_TT_pcqnm = bin2compartment_avg(nt_tt_pcqnm_df[['chr','start','end','TT','padj','sample_maha']],'TT')
# compartment_PT_pcqnm = bin2compartment_avg(pt_rt_pcqnm_df[['chr','start','end','PT','padj','sample_maha']],'PT')
# compartment_RT_pcqnm = bin2compartment_avg(pt_rt_pcqnm_df[['chr','start','end','RT','padj','sample_maha']],'RT')
# # save
# writer2 = pd.ExcelWriter('{}/compartments_combined.xlsx'.format(datdir), engine='xlsxwriter')
# compartment_NT_pcqnm.to_excel(writer2, sheet_name='NT', index=False)
# compartment_TT_pcqnm.to_excel(writer2, sheet_name='TT', index=False)
# compartment_PT_pcqnm.to_excel(writer2, sheet_name='PT', index=False)
# compartment_RT_pcqnm.to_excel(writer2, sheet_name='RT', index=False)
# writer2.save()
#
# compartment_NT_pcqnm['key'] = compartment_NT_pcqnm['chr'] + '_' + compartment_NT_pcqnm['start'].astype(str) + '_' \
#                               + compartment_NT_pcqnm['end'].astype(str)
# compartment_TT_pcqnm['key'] = compartment_TT_pcqnm['chr'] + '_' + compartment_TT_pcqnm['start'].astype(str) + '_' \
#                               + compartment_TT_pcqnm['end'].astype(str)
# compartment_PT_pcqnm['key'] = compartment_PT_pcqnm['chr'] + '_' + compartment_PT_pcqnm['start'].astype(str) + '_' \
#                               + compartment_PT_pcqnm['end'].astype(str)
# compartment_RT_pcqnm['key'] = compartment_RT_pcqnm['chr'] + '_' + compartment_RT_pcqnm['start'].astype(str) + '_' \
#                               + compartment_RT_pcqnm['end'].astype(str)
# NT_pcqnm_dict = dict(zip(compartment_NT_pcqnm['key'], compartment_NT_pcqnm['pcQnm.avg']))
# TT_pcqnm_dict = dict(zip(compartment_TT_pcqnm['key'], compartment_TT_pcqnm['pcQnm.avg']))
# PT_pcqnm_dict = dict(zip(compartment_PT_pcqnm['key'], compartment_PT_pcqnm['pcQnm.avg']))
# RT_pcqnm_dict = dict(zip(compartment_RT_pcqnm['key'], compartment_RT_pcqnm['pcQnm.avg']))
# # find switch
# NT_TT_switch = find_switch(compartment_NT_pcqnm,compartment_TT_pcqnm,'NT','TT')
# PT_RT_switch = find_switch(compartment_PT_pcqnm,compartment_RT_pcqnm,'PT','RT')
#
# # assign compartment score back
# NT_TT_switch['key.NT'] = NT_TT_switch['chr1'] + '_' + NT_TT_switch['start1'].astype(str) + '_' + NT_TT_switch['end1'].astype(str)
# NT_TT_switch['NT.pcQnm'] = NT_TT_switch['key.NT'].map(NT_pcqnm_dict)
# NT_TT_switch['key.TT'] = NT_TT_switch['chr2'] + '_' + NT_TT_switch['start2'].astype(str) + '_' + NT_TT_switch['end2'].astype(str)
# NT_TT_switch['TT.pcQnm'] = NT_TT_switch['key.TT'].map(TT_pcqnm_dict)
# PT_RT_switch['key.PT'] = PT_RT_switch['chr1'] + '_' + PT_RT_switch['start1'].astype(str) + '_' + PT_RT_switch['end1'].astype(str)
# PT_RT_switch['PT.pcQnm'] = PT_RT_switch['key.PT'].map(PT_pcqnm_dict)
# PT_RT_switch['key.RT'] = PT_RT_switch['chr2'] + '_' + PT_RT_switch['start2'].astype(str) + '_' + PT_RT_switch['end2'].astype(str)
# PT_RT_switch['RT.pcQnm'] = PT_RT_switch['key.RT'].map(RT_pcqnm_dict)
# NT_TT_switch.drop(['key.NT','key.TT'],inplace=True,axis=1)
# PT_RT_switch.drop(['key.PT','key.RT'],inplace=True,axis=1)
# # remove key
# writer3 = pd.ExcelWriter('{}/compartments_switch.xlsx'.format(datdir), engine='xlsxwriter')
# NT_TT_switch.to_excel(writer3, sheet_name='NT2TT', index=False)
# PT_RT_switch.to_excel(writer3, sheet_name='PT2RT', index=False)
# writer3.save()
#
# # In terms of the number of the switch
# def build_sankey(switch_df,groups):
#     comp1_count = switch_df[['chr1','start1','end1','{}_comp'.format(groups[0])]].drop_duplicates()['{}_comp'.format(groups[0])].value_counts()
#     comp2_count = switch_df[['chr1','start1','end1','{}_comp'.format(groups[0])]].drop_duplicates()['{}_comp'.format(groups[0])].value_counts()
#     # initiate sankey dict
#     # darkorange:A; darkviolet:B
#     sankey_dict = {
#         'source': [0, 0, 1, 1],
#         'target': [2, 3, 2, 3],
#         'value': [comp1_count['A'], comp1_count['B'], comp2_count['A'], comp2_count['B']],
#         'label': ['{}_A'.format(groups[0]), '{}_B'.format(groups[0]), '{}_A'.format(groups[1]), '{}_B'.format(groups[1])],
#         'color_link': ['darkorange', 'darkviolet', 'darkorange', 'darkviolet'],
#         'color_node': ['darkorange', 'darkviolet', 'darkorange', 'darkviolet']
#     }
#
#     a2b_count = len(switch_df[switch_df['SwitchType'] == 'A-B'])
#     b2a_count = len(switch_df[switch_df['SwitchType'] == 'B-A'])
#     sankey_dict['value'][0] = sankey_dict['value'][0] - a2b_count
#     sankey_dict['value'][1] = a2b_count
#     sankey_dict['value'][2] = b2a_count
#     sankey_dict['value'][3] = sankey_dict['value'][1] - b2a_count
#     return sankey_dict
#
# def plot_sankeyplot(sankey_dict,outname):
#     layout = go.Layout(
#         autosize=False,
#         width=500,
#         height=500
#     )
#     fig = go.Figure(data=[go.Sankey(
#         node=dict(
#             pad=15,
#             thickness=20,
#             line=dict(color="black", width=0.5),
#             label=['<b>' + label + '</b>' for label in sankey_dict['label']],
#             color=sankey_dict['color_node']
#         ),
#         link=dict(
#             source=sankey_dict['source'],  # indices correspond to labels, eg A1, A2, A1, B1, ...
#             target=sankey_dict['target'],
#             value=sankey_dict['value'],
#             color=sankey_dict['color_link']
#         ),
#         textfont=dict(size=16)
#     )], layout=layout)
#     fig.update_layout(title_text="{}-{} CompartmentA/B Transition".format(sankey_dict['label'][0].split('_')[0],
#                                                                           sankey_dict['label'][2].split('_')[0]),
#                       font_size=10)
#     fig.write_image(outname)
#
# nt_tt_sankey_dict = build_sankey(NT_TT_switch,['NT','TT'])
# pt_rt_sankey_dict = build_sankey(PT_RT_switch,['PT','RT'])
#
# # Fig. 1D
# plot_sankeyplot(nt_tt_sankey_dict,'{}/NT_TT_CompTrans_080623.png'.format(datdir))
# plot_sankeyplot(pt_rt_sankey_dict,'{}/PT_RT_CompTrans_080623.png'.format(datdir))

