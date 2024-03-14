import pandas as pd
import numpy as np
import pybedtools as pybt
import matplotlib.pyplot as plt
import seaborn as sns
import logomaker

# v4 change terminology
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
####### load files
hg19_f = '/Users/kfang/Documents/lab/Jin_lab/refseq/hg19/hg19.refseq.select.annot.noalt.srt.052621.bed'
module_f = '/Users/kfang/Documents/lab/Jin_lab/scRNA_scATAC/scRNA/annotation/Modules.xlsx'
outdir = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/ComptTADsLoops'

# load seven modules
ER_mod = pd.read_excel(module_f,sheet_name='ESR1',header=None)
HER2_mod = pd.read_excel(module_f,sheet_name='ERBB2',header=None)
proliferation_mod = pd.read_excel(module_f,sheet_name='AURKA',header=None)
invasion_mod = pd.read_excel(module_f,sheet_name='PLAU',header=None)
immune_mod = pd.read_excel(module_f,sheet_name='STAT1',header=None)
angiogenesis_mod = pd.read_excel(module_f,sheet_name='VEGF',header=None)
homeostasis_mod = pd.read_excel(module_f,sheet_name='Homeostasis',header=None)
module_list = ['ER','HER2','Proliferation','Invasion','ImmuneResponse','Angiogenesis','Homeostasis']
modules_dict = {'ER':ER_mod, 'HER2':HER2_mod,
                'Proliferation':proliferation_mod,
                'Invasion':invasion_mod,
                'ImmuneResponse':immune_mod,
                'Angiogenesis':angiogenesis_mod,
                'Homeostasis':homeostasis_mod}

# load hg19 genes
hg19_genes = pd.read_table(hg19_f,header=None)

# load compartments
comp_f = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/dcHic/100k/compartments_dcHiC.xlsx'
ptrt_compt_diff = pd.read_excel(comp_f,sheet_name='PT_RT_diff')

# load tads
tads_f = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/tads_calling/diff_tads/win5_40k/GISTA_test2/tads_sub_array_marks_gcut0.7_icut0.7.xlsx'
ptrt_tads_diff = pd.read_excel(tads_f,sheet_name='RTvsPT_GVH')

# acquire genes for compt and tads
def acquire_genes(genes_bt,input_bt):
    input_genes = pybt.BedTool.intersect(input_bt, genes_bt, wa=True, wb=True).to_dataframe(
        disable_auto_names=True,
        names=['chr1', 'start1', 'end1',
               'chr2', 'start2', 'end2', 'gene'])
    return input_genes

genes_bt = pybt.BedTool.from_dataframe(hg19_genes[[0,1,2,4]])
ptrt_compt_diff_bt = pybt.BedTool.from_dataframe(ptrt_compt_diff[['chr','start','end']])
ptrt_tads_diff_bt = pybt.BedTool.from_dataframe(ptrt_tads_diff[['chr','start(min)','end(max)']])

ptrt_compt_diff_genes = acquire_genes(genes_bt,ptrt_compt_diff_bt)
ptrt_tads_diff_genes = acquire_genes(genes_bt,ptrt_tads_diff_bt)

# assign back diff information
ptrt_compt_diff_genes = pd.merge(ptrt_compt_diff_genes, ptrt_compt_diff[['chr','start','end','SwitchType']],
                                 left_on=['chr1','start1','end1'], right_on=['chr','start','end'])
ptrt_tads_diff_genes = pd.merge(ptrt_tads_diff_genes,
                                ptrt_tads_diff[['chr','start(min)','end(max)','GroupChange','ChangeTypes','IndividualChange']],
                                left_on=['chr1','start1','end1'],right_on=['chr','start(min)','end(max)'])

name2name = {'Indi-Var-H':'IVH','Indi-Var-M&C':'IMC'}
ptrt_tads_diff_genes['SubTypes'] = ptrt_tads_diff_genes['IndividualChange'].map(name2name) + '-' + ptrt_tads_diff_genes['ChangeTypes']
ptrt_tads_diff_genes[['chr1','start1','end1','chr2','start2','end2','gene']].to_csv('{}/RTPT_tads_gene_annot.csv'.format(outdir))
# one gene only one type of change in three organization level
def find_across_genes(diff_genes_df,gene_col,type_col):
    gene_count = diff_genes_df[[gene_col, type_col]].drop_duplicates()[gene_col].value_counts()
    return gene_count[gene_count>1].index.values

ptrt_compt_across_genes = find_across_genes(ptrt_compt_diff_genes, 'gene', 'SwitchType')
ptrt_tads_across_genes = find_across_genes(ptrt_tads_diff_genes,'gene','SubTypes')

def remove_across(diff_genes_df, acrossed_genes, gene_col):
    df = diff_genes_df.copy()
    df_remained = df[~df[gene_col].isin(acrossed_genes)]
    df_removed = df[df[gene_col].isin(acrossed_genes)]
    df_removed['intersect'] = df_removed[['end1','end2']].min(axis=1) - df_removed[['start1','start2']].max(axis=1)
    df_removed = df_removed.sort_values('intersect', ascending=False).drop_duplicates(['gene'])
    df_removed.drop(columns='intersect',inplace=True)
    df_final = pd.concat([df_remained,df_removed]).sort_index().reset_index(drop=True)
    return df_final

ptrt_compt_diff_genes_filt = remove_across(ptrt_compt_diff_genes, ptrt_compt_across_genes, 'gene')

ptrt_tads_diff_genes_filt = remove_across(ptrt_tads_diff_genes, ptrt_tads_across_genes, 'gene')

# build gene feature dataframe
allsubtypes = ['A-A','A-B','B-A','B-B','IVH-ND','IVH-SF','IVH-Mixed','IMC-ND','IMC-SF','IMC-Mixed','Gained','Lost']
# RTvsPT
ptrt_genes_features_df = pd.DataFrame(0,index=hg19_genes[4],columns=allsubtypes)
for idx, row in ptrt_compt_diff_genes_filt[['gene','SwitchType']].drop_duplicates().iterrows():
    ptrt_genes_features_df.loc[row['gene'],row['SwitchType']] += 1

for idx, row in ptrt_tads_diff_genes_filt[['gene','SubTypes']].drop_duplicates().iterrows():
    ptrt_genes_features_df.loc[row['gene'], row['SubTypes']] += 1


ptrt_genes_features_c2_df = ptrt_genes_features_df.copy()
loop_f = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/Individual_specific_loops/RTsIndvsPT_commons.xlsx'
ptrt_loop_all_gained = pd.read_excel(loop_f,sheet_name='>=1_RTs_Gained')
ptrt_loop_all_gained.columns = ['genes'] + list(ptrt_loop_all_gained.columns.values[1:])
for idx,row in ptrt_loop_all_gained.iterrows():
    ptrt_genes_features_df.loc[row['genes'],'Gained'] += 1

ptrt_loop_all_lost = pd.read_excel(loop_f,sheet_name='>=1_RTs_Lost')
ptrt_loop_all_lost.columns = ['genes'] + list(ptrt_loop_all_lost.columns.values[1:])
for idx,row in ptrt_loop_all_lost.iterrows():
    ptrt_genes_features_df.loc[row['genes'],'Lost'] += 1

ptrt_genes_features_df = ptrt_genes_features_df[ptrt_genes_features_df.sum(axis=1)!=0]
ptrt_genes_features_df= ptrt_genes_features_df[ptrt_genes_features_df[['Gained','Lost']].sum(axis=1)!=2]
ptrt_genes_features_df.to_csv('{}/RTPT_gene_feature.csv'.format(outdir))

# categorize looping gene into only in 1RTs_specific,>2RTs
ptrt_loop_c2_gained = pd.read_excel(loop_f,sheet_name='>=2_RTs_Gained')
ptrt_loop_c2_gained.columns = ['genes'] + list(ptrt_loop_c2_gained.columns.values[1:])
for idx,row in ptrt_loop_c2_gained.iterrows():
    ptrt_genes_features_c2_df.loc[row['genes'],'Gained'] += 1

ptrt_loop_c2_lost = pd.read_excel(loop_f,sheet_name='>=2_RTs_Lost')
ptrt_loop_c2_lost.columns = ['genes'] + list(ptrt_loop_c2_lost.columns.values[1:])
for idx,row in ptrt_loop_c2_lost.iterrows():
    ptrt_genes_features_c2_df.loc[row['genes'],'Lost'] += 1

ptrt_genes_features_c2_df = ptrt_genes_features_c2_df[ptrt_genes_features_c2_df.sum(axis=1)!=0]
c2_genes = list(set(ptrt_genes_features_c2_df.index.values) & set(ptrt_genes_features_df.index.values))
ptrt_genes_features_c2_df = ptrt_genes_features_c2_df.loc[c2_genes,]
# heatmap visualize
tads_subtypes = ['IVH-ND','IVH-SF','IVH-Mixed','IMC-ND','IMC-SF','IMC-Mixed']
loop_subtypes = ['Gained','Lost']
tads_loop_c1_df = pd.DataFrame(0,index=loop_subtypes,columns=tads_subtypes)
tads_loop_c2_df = pd.DataFrame(0,index=loop_subtypes,columns=tads_subtypes)
with pd.ExcelWriter('{}/tads_loop_genes.xlsx'.format(outdir)) as writer:
    for idx, type_x in enumerate(loop_subtypes):
        for idy,type_y in enumerate(tads_subtypes):
            current_df = ptrt_genes_features_df[ptrt_genes_features_df[[type_x, type_y]].sum(axis=1) == 2]
            tads_loop_c1_df.loc[type_x,type_y] = len(current_df)
            current_df.to_excel(writer, sheet_name='>=1RTs-{}_{}'.format(type_x,type_y))

    for idx, type_x in enumerate(loop_subtypes):
        for idy, type_y in enumerate(tads_subtypes):
            current_df = ptrt_genes_features_c2_df[ptrt_genes_features_c2_df[[type_x, type_y]].sum(axis=1) == 2]
            tads_loop_c2_df.loc[type_x, type_y] = len(current_df)
            current_df.to_excel(writer, sheet_name='>=2RTs-{}_{}'.format(type_x, type_y))

fig, axs = plt.subplots(2,1,sharex=True,figsize=(10,8))
sns.heatmap(tads_loop_c1_df,square=True,cmap='rocket_r',edgecolor='white',
            linewidths=1,annot=True,cbar=False,ax=axs[0],annot_kws={'fontsize':14,'weight':'bold'},fmt='d')
axs[0].set_yticklabels(axs[0].get_yticklabels(),fontsize=16)
sns.heatmap(tads_loop_c2_df,square=True,cmap='rocket_r',edgecolor='white',
            linewidths=1,annot=True,cbar=False,ax=axs[1],annot_kws={'fontsize':14,'weight':'bold'},fmt='d')
axs[1].set_yticklabels(axs[1].get_yticklabels(),fontsize=16)
axs[1].set_xticklabels(['HS-ND','HS-SF','HS-Mixed','LS-ND','LS-SF','LS-Mixed'],fontsize=16)
plt.savefig('{}/2type_looping_heatmap.png'.format(outdir),dpi=500)
plt.close()


subtype2org = {'A-A':'Compartment','A-B':'Compartment','B-A':'Compartment','B-B':'Compartment',
               'IVH-ND':'TAD','IVH-SF':'TAD','IVH-Mixed':'TAD','IMC-ND':'TAD','IMC-SF':'TAD','IMC-Mixed':'TAD',
               'Gained':'Loop','Lost':'Loop'}

tab10 = sns.color_palette()

# {'A-A':'darkorange','A-B':'mediumorchid','B-A':'lightsalmon','B-B':'darkviolet'},
# {'IVH-ND':tab10[9],'IVH-SF':tab10[8],'IVH-Mixed':tab10[7],'IMC-ND':tab10[6],'IMC-SF':tab10[5],'IMC-Mixed':tab10[4]}
# {'Gained':tab10[0],'Lost':tab10[1]}

# plot circle bar plot Fig4A
# Load the dataset
data = ptrt_genes_features_df.copy()
# rename
data = data.rename(columns={'IVH-ND': 'HS-ND', 'IVH-SF': 'HS-SF', 'IVH-Mixed': 'HS-Mixed',
                            'IMC-ND': 'LS-ND', 'IMC-SF': 'LS-SF', 'IMC-Mixed': 'LS-Mixed'})
# List of dataframes and categories
dataframes = [data[['A-A', 'A-B', 'B-A', 'B-B']],
              data[['HS-ND', 'HS-SF', 'HS-Mixed', 'LS-ND', 'LS-SF', 'LS-Mixed']],
              data[['Gained', 'Lost']]]
categories_list = [['A-A', 'A-B', 'B-A', 'B-B'],
                   ['HS-ND', 'HS-SF', 'HS-Mixed', 'LS-ND', 'LS-SF', 'LS-Mixed'],
                   ['Gained', 'Lost']]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
          "#2ea02c", "#d62758"]

def draw_custom_yaxis(ax, position, color, max_val, label):
    """Draw a custom y-axis on a polar plot."""
    y_ticks = np.linspace(0, max_val, 5)
    for y in y_ticks:
        ax.plot([position, position], [0, y], color=color, lw=0.5)
        ax.text(position, y, str(int(y)), color=color, ha='center', va='center')
    ax.text(position, max_val + 0.05 * max_val, label, color=color, ha='center', va='center')


def combined_circular_barplot_log_annotated(dataframes, original_dataframes, categories_list, colors, ax, log_scale,
                                            max_radius=1.5, offset=1.5):
    """Version with log10 scale and raw value annotations."""
    all_radii = []
    if log_scale:
        for df in dataframes:
            all_radii.extend(np.log10(df.sum().values))
    else:
        for df in dataframes:
            all_radii.extend(df.sum().values)

    all_original_values = []
    for df in original_dataframes:
        all_original_values.extend(df.sum().values)

    all_categories = []
    for cat_list in categories_list:
        all_categories.extend(cat_list)

    total_bars = len(all_categories)
    all_theta = np.linspace(0, 2 * np.pi, total_bars, endpoint=False)
    width = 2 * np.pi / total_bars

    bars = ax.bar(all_theta, all_radii, width=width, align='center', color=colors, edgecolor="k", alpha=0.6)

    # Annotate bars with their original values
    for bar, value in zip(bars, all_original_values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height()-0.5,
                str(int(value)), ha='center', va='bottom',fontsize=16)

    ax.set_xticks(all_theta)
    # Increase font size and offset for xticklabels
    for theta, label in zip(all_theta, all_categories):
        ax.text(theta, max_radius * offset, label, ha='center', va='center', fontsize=18)
    ax.set_xticklabels([])
    ax.yaxis.grid(False)
    ax.set_rlabel_position(90)
    ax.set_rticks([])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    return bars


# Create the figure
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})

# Create combined circular bar plot with log10 scale and raw value annotations
combined_circular_barplot_log_annotated(dataframes, dataframes, categories_list, colors, ax, log_scale=True,
                                        max_radius=3, offset=1.55)

# # Draw custom y-axes
# draw_custom_yaxis(ax, np.pi, colors[0], 3, "Compartment (log10)")
# draw_custom_yaxis(ax, np.pi / 2, colors[len(categories_list[0])], 3, "TAD (log10)")
# draw_custom_yaxis(ax, -np.pi / 2, colors[len(categories_list[0]) + len(categories_list[1])], 3, "Loop (log10)")
# plt.title("Combined Gene Features Distribution (log10 scale with raw annotations)", y=1.1)
plt.tight_layout()
plt.savefig('{}/genecount_circle_barplot'.format(outdir),dpi=500)
plt.close()

# correlated with module for 12 categories
# Initialize a dictionary to store the count of common genes between each column of RTPT_gene_feature and each subsheet of Modules.xlsx
common_genes_count = {}
module_sizes = {mod:len(modules_dict[mod]) for mod in module_list}
# Iterate over each column (excluding the identifier column) in RTPT_gene_feature.csv
for column in ptrt_genes_features_df.columns:
    common_genes_count[column] = []

    # Extract genes where the feature value is 1
    genes_from_rtpt = ptrt_genes_features_df[ptrt_genes_features_df[column] == 1].index.values.tolist()

    # Iterate over each subsheet in Modules.xlsx
    for mod in module_list:
        genes_from_module = modules_dict[mod].iloc[:, 0].tolist()

        # Count common genes
        common_count = len(set(genes_from_rtpt).intersection(set(genes_from_module)))
        common_genes_count[column].append(common_count)

# Convert the dictionary to a DataFrame for easier visualization
common_genes_df = pd.DataFrame(common_genes_count, index=module_list)
normalized_common_genes_df = common_genes_df.divide(pd.Series(module_sizes), axis=0)*100
normalized_common_genes_df_melt = pd.melt(normalized_common_genes_df.reset_index(), id_vars='index')

def heatmap(x, y, size, color, outname=None, x_order=None, y_order=None, marker='s',figsize=(16,9),point_size=50):
    plt.figure(figsize=figsize)
    def _value_to_color(val):
        val_position = float((val - color_min)) / (color_max - color_min)
        ind = int(val_position * (n_colors - 1))  # target index in the color palette
        return palette[ind]

    plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1)  # Setup a 1x15 grid
    ax = plt.subplot(plot_grid[:, :-1])  # Use the leftmost  columns of the grid for the main plot

    # Mapping from column names to integer coordinates
    if x_order is None:
        x_labels = [v for v in sorted(x.unique())]
    else:
        x_labels = x_order

    if y_order is None:
        y_labels = [v for v in sorted(y.unique())]
    else:
        y_labels = y_order
    x_to_num = {p[1]: p[0] for p in enumerate(x_labels)}
    y_to_num = {p[1]: p[0] for p in enumerate(y_labels)}

    n_colors = 256  # Use 256 colors for the diverging color palette
    palette = sns.diverging_palette(220, 20, n=n_colors)  # Create the palette
    # Range of values that will be mapped to the palette, i.e. min and max possible correlation
    color_min, color_max = [np.min(color), np.max(color)]

    # position of value in the input range, relative to the length of the input range
    size_scale = point_size
    ax.scatter(
        x=x.map(x_to_num),  # Use mapping for x
        y=y.map(y_to_num),  # Use mapping for y
        s=size * size_scale,  # Vector of square sizes, proportional to size parameter
        marker=marker,  # Use square as scatterplot marker
        c=color.apply(_value_to_color),  # Vector of square color values, mapped to color palette
    )

    # Show column labels on the axes
    ax.set_xticks([x_to_num[v] for v in x_labels])
    ax.set_xticklabels(x_labels, rotation=45, horizontalalignment='right',fontsize=20, fontweight='bold')
    ax.set_yticks([y_to_num[v] for v in y_labels])
    ax.set_yticklabels(y_labels,fontsize=20, fontweight='bold')
    # put int grid center
    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)
    # add color bar
    # Add color legend on the right side of the plot
    ax = plt.subplot(plot_grid[:, -1])  # Use the rightmost column of the plot

    col_x = [0] * len(palette)  # Fixed x coordinate for the bars
    bar_y = np.linspace(color_min, color_max, n_colors)  # y coordinates for each of the n_colors bars

    bar_height = bar_y[1] - bar_y[0]
    ax.barh(
        y=bar_y,
        width=[5] * len(palette),  # Make bars 5 units wide
        left=col_x,  # Make bars start at 0
        height=bar_height,
        color=palette,
        linewidth=0
    )
    ax.set_xlim(1, 2)  # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
    ax.grid(False)  # Hide grid
    ax.set_facecolor('white')  # Make background white
    ax.set_xticks([])  # Remove horizontal ticks
    ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3))  # Show vertical ticks for min, middle and max
    ax.yaxis.tick_right()  # Show vertical ticks on the right
    if outname is None:
        plt.show()
    else:
        plt.savefig(outname,bbox_inches='tight',dpi=500)
        plt.close()

heatmap(
    x=normalized_common_genes_df_melt['variable'],
    y=normalized_common_genes_df_melt['index'],
    size=normalized_common_genes_df_melt['value'].abs(),
    color= normalized_common_genes_df_melt['value'],
    marker='h',
    x_order = normalized_common_genes_df.columns,
    y_order = np.flip(module_list),
    figsize=(16,11),
    outname='{}/module_12types.png'.format(outdir)
)

# Fig4C is done by network_visual


# module for top 5
top5_edges = [['IVH-Mixed','Gained'],['IVH-SF','Gained'],['IVH-Mixed','Lost'],['IVH-SF','Lost'],['IMC-Mixed','Gained']]
top5_common_genes_count = {'_'.join(edge):[] for edge in top5_edges}
for edge in top5_edges:
    # Extract genes where the feature value is 1
    edge_genes = ptrt_genes_features_df[ptrt_genes_features_df[edge].sum(axis=1)==2].index.values

    # Iterate over each subsheet in Modules.xlsx
    for mod in module_list:
        genes_from_module = modules_dict[mod].iloc[:, 0].tolist()

        # Count common genes
        common_count = len(set(edge_genes).intersection(set(genes_from_module)))
        top5_common_genes_count['_'.join(edge)].append(common_count)

# Convert the dictionary to a DataFrame for easier visualization
top5_common_genes_count_df = pd.DataFrame(top5_common_genes_count, index=module_list)
normalized_top5_common_genes_count_df = top5_common_genes_count_df.divide(pd.Series(module_sizes), axis=0)*100
normalized_top5_common_genes_count_df_melt = pd.melt(normalized_top5_common_genes_count_df.reset_index(), id_vars='index')

heatmap(
    x=normalized_top5_common_genes_count_df_melt['variable'],
    y=normalized_top5_common_genes_count_df_melt['index'],
    size=normalized_top5_common_genes_count_df_melt['value'].abs(),
    color= normalized_top5_common_genes_count_df_melt['value'],
    marker='h',
    x_order = normalized_top5_common_genes_count_df.columns,
    y_order = np.flip(module_list),
    figsize=(9,8),
    outname='{}/module_top5_connection.png'.format(outdir)
)

# module for assinged edge
assign_edges = [['IVH-Mixed','Gained'],['IVH-SF','Gained'],['IVH-ND','Gained'],['IVH-SF','Lost'],['IVH-Mixed','Lost']]
edge_xls = pd.ExcelFile('{}/tads_loop_genes.xlsx'.format(outdir))
edge_genes = {}
for edge in assign_edges:
    edge_genes['>=1RTs-{}_{}'.format(edge[1], edge[0])] = pd.read_excel(edge_xls,
                                                                    sheet_name='>=1RTs-{}_{}'.format(edge[1], edge[0]),
                                                                    index_col=0).index.values
    edge_genes['>=2RTs-{}_{}'.format(edge[1], edge[0])] = pd.read_excel(edge_xls,
                                                                    sheet_name='>=2RTs-{}_{}'.format(edge[1], edge[0]),
                                                                    index_col=0).index.values

assign_edges_module = {'>=1RTs':{'{}_{}'.format(edge[0],edge[1]):[] for edge in assign_edges},
                       '>=2RTs':{'{}_{}'.format(edge[0],edge[1]):[] for edge in assign_edges}}
for edge in assign_edges:
    for cond in assign_edges_module:
        current_genes = edge_genes['{}-{}_{}'.format(cond,edge[1],edge[0])]
        # Iterate over each subsheet in Modules.xlsx
        for mod in module_list:
            genes_from_module = modules_dict[mod].iloc[:, 0].tolist()
            # Count common genes
            common_count = len(set(current_genes).intersection(set(genes_from_module)))
            assign_edges_module[cond]['{}_{}'.format(edge[0],edge[1])].append(common_count)

# Convert the dictionary to a DataFrame for easier visualization
common_genes_count_1rts = pd.DataFrame(assign_edges_module['>=1RTs'], index=module_list)
common_genes_count_2rts = pd.DataFrame(assign_edges_module['>=2RTs'], index=module_list)
normalized_common_genes_count_1rts = common_genes_count_1rts.divide(pd.Series(module_sizes), axis=0)*100
normalized_common_genes_count_1rts_melt = pd.melt(normalized_common_genes_count_1rts.reset_index(), id_vars='index')
normalized_common_genes_count_2rts = common_genes_count_2rts.divide(pd.Series(module_sizes), axis=0)*100
normalized_common_genes_count_2rts_melt = pd.melt(normalized_common_genes_count_2rts.reset_index(), id_vars='index')
common_genes_count_1rts_melt = pd.melt(common_genes_count_1rts.reset_index(), id_vars='index')
common_genes_count_2rts_melt = pd.melt(common_genes_count_2rts.reset_index(), id_vars='index')

heatmap(
    x=common_genes_count_1rts_melt['variable'],
    y=common_genes_count_1rts_melt['index'],
    size=common_genes_count_1rts_melt['value'].abs(),
    color= common_genes_count_1rts_melt['value'],
    marker='h',
    x_order = ['HS-ND_Gained','HS-SF_Gained','HS-Mixed_Gained','HS-SF_Lost','HS-Mixed_Lost'],
    y_order = np.flip(module_list),
    figsize=(6,3),
    point_size=20,
    outname='{}/module_assigned_connection_1rts_heatmap.png'.format(outdir)
)

heatmap(
    x=common_genes_count_2rts_melt['variable'],
    y=common_genes_count_2rts_melt['index'],
    size=common_genes_count_2rts_melt['value'].abs(),
    color= common_genes_count_2rts_melt['value'],
    marker='h',
    x_order = ['HS-ND_Gained','HS-SF_Gained','HS-Mixed_Gained','HS-SF_Lost','HS-Mixed_Lost'],
    y_order = np.flip(module_list),
    figsize=(6,3),
    point_size=20,
    outname='{}/module_assigned_connection_2rts_heatmap.png'.format(outdir)
)