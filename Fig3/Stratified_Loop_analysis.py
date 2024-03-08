import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
from functools import reduce
from collections import Counter
from scipy.stats import zscore
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

pd.options.mode.chained_assignment = None

res = 20000
cutoff = 0.98
tissue_outdir2 = '/data/kun/Lava/SLG-DEGs/fithic/tissue/{}/'.format(res)
nt_samples = ['NT1','NT2']
pt_samples = ['PT1','PT2','PT3','PT4','PT5']
rt_samples = ['RT1','RT2','RT3','RT4','RT5']
tissue_samples = nt_samples + pt_samples + rt_samples
cell_outdir2 = '/data/kun/Lava/SLG-DEGs/fithic/cell/{}/'.format(res)
cell_samples = ['MCF7','MCF7TR_TCC','T47D','T47DTR','ZR75-1','ZR75-1TR']
geneanno = '/Users/KunFang/Documents/lab/Jin_lab/refseq/hg19.refseq.select.annot.noalt.srt.052621.bed'
# part1 outdir
outdir1 = '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/filt_loops'
# part2 outdir
outdir2 = '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/annot_loops'
# part3 outdir
outdir3 = '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/GC_loops'

################### part 1. apply qvalue_cutoff to filter some insignificant interactions ######################
for t_sample in tissue_samples:
    print(t_sample)
    current_df = pd.read_csv(tissue_outdir2+"/{}/{}_fithic.spline_pass1.res{}.significances.txt.gz".format(t_sample,t_sample,res),sep='\t')
    contact_q3 = current_df['contactCount'].quantile(cutoff)
    current_df_q3 = current_df[current_df['contactCount'] >= contact_q3]
    current_df_q3.to_csv('{}/{}_loops_q{}_filt.txt'.format(outdir1,t_sample,cutoff),sep='\t',index=False)

for c_sample in cell_samples:
    print(c_sample)
    current_df = pd.read_csv(cell_outdir2+"/{}/{}_fithic.spline_pass1.res{}.significances.txt.gz".format(c_sample,c_sample,res),sep='\t')
    contact_q3 = current_df['contactCount'].quantile(cutoff)
    current_df_q3 = current_df[current_df['contactCount']>=contact_q3]
    current_df_q3.to_csv('{}/{}_loops_q{}_filt.txt'.format(outdir1,t_sample,cutoff),sep='\t',index=False)

############### part 2. annotate filtered interactions and extract P1D1 loops from them ####################
q_val = cutoff
file_suffix = '_loops_q{}_filt.txt'.format(q_val)
# unit 10M valid pairs
# unit 10M valid pairs
seq_depth = {'NT1':15.33,'NT2':10.47,'PT1':14.88,'PT2':15.03,'PT3':7.57,'PT4':11.58,'PT5':9.10,
             'RT1':10.08,'RT2':10.85,'RT3':7.79,'RT4':10.18,'RT5':11.83,'MCF7':18.62,'MCF7TR':5.51,'T47D':14.16,
             'T47DTR':16.21,'ZR75-1':14.01,'ZR75-1TR':11.77}
resolution = 20000
distal_up_up = 200000
distal_up_down = 10000
distal_down_up = -10000
distal_down_down = -200000
promoter_up = 4000
promoter_down = -2000
D_overlap_cutoff = 1
P_overlap_cutoff = promoter_up/(2*resolution)

# half length
merge_overlap_cutoff = 0.5
only_DP = True

def region_df(up_b,down_b,gene_anno_df,anno_str):
    # divide gene list df to plus and minus strand
    plus_gene_anno_df = gene_anno_df[gene_anno_df[3]=='+'][[0,1,1,3,4,5]]
    minus_gene_anno_df = gene_anno_df[gene_anno_df[3]=='-'][[0,2,2,3,4,5]]
    plus_gene_anno_df.columns =[0,1,2,3,4,5]
    minus_gene_anno_df.columns =[0,1,2,3,4,5]
    # make the distal or promoter region annotation file
    plus_gene_region_df = plus_gene_anno_df[[1,2]].sub([up_b,down_b],axis='columns')
    minus_gene_region_df = minus_gene_anno_df[[1,2]].add([down_b,up_b],axis='columns')
    plus_gene_region_df = plus_gene_region_df.clip(lower=0)
    minus_gene_region_df = minus_gene_region_df.clip(lower=0)
    plus_gene_region_df[0] = gene_anno_df[0]
    plus_gene_region_df[3] = gene_anno_df[4]
    plus_gene_region_df[4] = anno_str
    plus_gene_region_df = plus_gene_region_df[[0,1,2,3,4]]
    minus_gene_region_df[0] = gene_anno_df[0]
    minus_gene_region_df[3] = gene_anno_df[4]
    minus_gene_region_df[4] = anno_str
    minus_gene_region_df = minus_gene_region_df[[0,1,2,3,4]]
    gene_region_df = pd.concat([plus_gene_region_df,minus_gene_region_df]).sort_index()
    return gene_region_df

def load_loop_f(samples,file_path,file_suffix,seq_depth):
    samples_df_dict = {}
    select_col = ['chr1','start1','end1','chr2','start2','end2','contactCount']
    for sample in samples:
        print(sample)
        tmp_file = "{}/{}{}".format(file_path,sample,file_suffix)
        sif_df = pd.read_table(tmp_file,sep='\t')
        sif_df = sif_df[select_col]
        # VPPM
        sif_df['NormC'] = np.round(sif_df['contactCount']/seq_depth[sample],2)
        samples_df_dict[sample] = sif_df
    return samples_df_dict

class Loop_Anno(object):
    def __init__(self,loop_df,gene_upD_df,gene_downD_df,gene_P_df,D_overlap_cutoff,P_overlap_cutoff):
        # loop_df should have at least 7 cols: chr1/start1/end1/chr2/start2/end2/VPPM
        # make sure all coordinates are greate than 0
        loop_df_copy = loop_df.copy()
        loop_df_copy[['start1','end1','start2','end2']][loop_df_copy[['start1','end1','start2','end2']]<0] = 0
        # make coordinates are all int type
        loop_df_copy[['start1','end1','start2','end2']] = loop_df_copy[['start1','end1','start2','end2']].astype(int)
        loop_df_copy['loop_id'] = np.arange(len(loop_df_copy))
        self.loop_df = loop_df_copy
        self.gene_upD_df = gene_upD_df
        self.gene_downD_df = gene_downD_df
        self.gene_P_df = gene_P_df
        self.D_overlap_cutoff = D_overlap_cutoff
        self.P_overlap_cutoff = P_overlap_cutoff

    def _intersect_files_df(self,sif_end_df,region_df,region_overlap_cutoff):
        # intersect the sif-end df and region annotation df and then transfer the intersection result into a dictionary.
        intersect_file1 = BedTool.from_dataframe(sif_end_df)
        intersect_file2 = BedTool.from_dataframe(region_df)
        result_file = BedTool.intersect(intersect_file1,intersect_file2,wa=True,wb=True,f=region_overlap_cutoff)
        result_inter_df = result_file.to_dataframe()
        result_inter_df.columns = [['chr1','start1','end1','loop_id','chr2','start2','end2','gene','region']]
        result_inter_df = result_inter_df[['chr1','start1','end1','loop_id','gene','region']]
        return result_inter_df

    def annot(self):
        # make sure all coordinates are greate than 0
        sif_end1_df, sif_end2_df = self.loop_df[['chr1','start1','end1','loop_id']],self.loop_df[['chr2','start2','end2','loop_id']]
        # sif end1 intersection
        print("Perform sif end1 D/P intersction..")
        sif_end1_upD_df = self._intersect_files_df(sif_end1_df,self.gene_upD_df,D_overlap_cutoff)
        sif_end1_downD_df = self._intersect_files_df(sif_end1_df,self.gene_downD_df,D_overlap_cutoff)
        sif_end1_D_df = pd.concat([sif_end1_upD_df,sif_end1_downD_df])
        sif_end1_P_df = self._intersect_files_df(sif_end1_df,self.gene_P_df,P_overlap_cutoff)

        # sif end2 intersection
        print("Perform sif end2 D/P intersction..")
        sif_end2_upD_df = self._intersect_files_df(sif_end2_df,self.gene_upD_df,D_overlap_cutoff)
        sif_end2_downD_df = self._intersect_files_df(sif_end2_df,self.gene_downD_df,D_overlap_cutoff)
        sif_end2_D_df = pd.concat([sif_end2_upD_df,sif_end2_downD_df])
        sif_end2_P_df = self._intersect_files_df(sif_end2_df,self.gene_P_df,P_overlap_cutoff)

        # rename cols
        sif_end1_D_df.columns = ['chr','start','end','loop_id','gene','region']
        sif_end1_P_df.columns = ['chr','start','end','loop_id','gene','region']
        sif_end2_D_df.columns = ['chr','start','end','loop_id','gene','region']
        sif_end2_P_df.columns = ['chr','start','end','loop_id','gene','region']

        print("Finding D-P and P-D matched end..")
        # merge loop end based on loop id: we focus on sif_end1_D_sif_end2_P or sif_end1_P_sif_end2_D
        end1_D_end2_P = pd.merge(sif_end1_D_df,sif_end2_P_df,how="inner",left_on=['gene','loop_id'],right_on=['gene','loop_id'])
        end1_P_end2_D = pd.merge(sif_end1_P_df,sif_end2_D_df,how="inner",left_on=['gene','loop_id'],right_on=['gene','loop_id'])

        all_annot = pd.concat([end1_D_end2_P,end1_P_end2_D])
        id_contactcount = dict(zip(self.loop_df['loop_id'],self.loop_df['contactCount']))
        id_VPPM = dict(zip(self.loop_df['loop_id'],self.loop_df['NormC']))
        all_annot['contactCount'] = all_annot['loop_id'].map(id_contactcount)
        all_annot['NormC'] = all_annot['loop_id'].map(id_VPPM)
        all_annot['region'] = all_annot['gene'] + '_' + all_annot['region_x'] + '_' + all_annot['region_y']
        all_annot.sort_values(by='loop_id',ascending=True,inplace=True)
        all_annot = all_annot[['chr_x','start_x','end_x','chr_y','start_y','end_y','region','contactCount','NormC']]
        return all_annot


print('Loading files..')
gene_anno_df = pd.read_table(geneanno,sep='\t',header=None)
gene_up_distal_df = region_df(distal_up_up,distal_up_down,gene_anno_df,'UpD')
gene_down_distal_df = region_df(distal_down_up,distal_down_down,gene_anno_df,'DownD')
gene_promoter_df = region_df(promoter_up,promoter_down,gene_anno_df,'P')
tissue_samples = nt_samples + pt_samples + rt_samples
tissue_samples_df_dict = load_loop_f(tissue_samples,outdir1,file_suffix,seq_depth)
cell_samples_df_dict = load_loop_f(cell_samples,outdir1,file_suffix,seq_depth)

for t_sample in tissue_samples:
    current_df_t = tissue_samples_df_dict[t_sample]
    loop_annot = Loop_Anno(current_df_t,gene_up_distal_df,gene_down_distal_df,gene_promoter_df,D_overlap_cutoff,
                           P_overlap_cutoff)
    annot_df = loop_annot.annot()
    annot_df = annot_df.groupby(['chr_x','start_x','end_x','chr_y','start_y','end_y','contactCount','NormC']).agg(lambda x:','.join(x)).reset_index()
    annot_df.to_csv("{}/{}_loops_q{}_filt_annot.txt".format(outdir2,t_sample,q_val),sep='\t',index=False)
    annot_df_simple = annot_df[['chr_x','start_x','end_x','chr_y','start_y','end_y']].drop_duplicates()
    print("Original {} loops' number : {}; P1D1 loops' abs number: {}".format(t_sample,len(current_df_t),
                                                                              len(annot_df_simple)))

for c_sample in cell_samples:
    current_df_c = cell_samples_df_dict[c_sample]
    loop_annot = Loop_Anno(current_df_c,gene_up_distal_df,gene_down_distal_df,gene_promoter_df,D_overlap_cutoff,
                           P_overlap_cutoff)
    annot_df = loop_annot.annot()
    annot_df = annot_df.groupby(['chr_x','start_x','end_x','chr_y','start_y','end_y','contactCount','NormC']).agg(lambda x:','.join(x)).reset_index()
    annot_df.to_csv("{}/{}_loops_q{}_filt_annot.txt".format(outdir2,c_sample,q_val),sep='\t',index=False)
    annot_df_simple = annot_df[['chr_x','start_x','end_x','chr_y','start_y','end_y']].drop_duplicates()
    print("Original {} loops' number : {}; P1D1 loops' abs number: {}".format(c_sample,len(current_df_c),
                                                                              len(annot_df_simple)))

##################### part3. transfer P1D1 loops to the gene-centric loop ###########################
file_suffix = '_loops_q{}_filt_annot.txt'.format(cutoff)
tissue_samples = nt_samples + pt_samples + rt_samples
tissue_gc_loops_dict = {}
gc_all_info = {'Sample':[],'The number of Genes':[]}
all_gc_normc_list = []
all_gc_count_list = []
for t_sample in tissue_samples:
    print("Processing {}".format(t_sample))
    with open("{}/{}{}".format(outdir2,t_sample,file_suffix),'r') as current_f:
        # dictionary: 'Gene':{'UpDP_Signal':[],'DownDP_Signal':[],'Total_Signal':[]}
        current_dict = {}
        # Skip header
        header_line = current_f.readline()
        for line in current_f:
            line_info = line.strip().split()
            line_regions = line_info[-1].split(',')
            line_normC = float(line_info[-2])
            line_count = int(line_info[-3])
            for region in line_regions:
                current_gene = region.split('_')[0]
                updown = region.split('_')[1:]
                if current_gene in current_dict:
                    current_dict[current_gene]['Loops_Count'] += 1
                    if 'UpD' in updown:
                        current_dict[current_gene]['UpDP_NormC'] += line_normC
                        current_dict[current_gene]['UpDP_Count'] += line_count
                    elif 'DownD' in updown:
                        current_dict[current_gene]['DownDP_NormC'] += line_normC
                        current_dict[current_gene]['DownDP_Count'] += line_count
                    else:
                        print("Abnormal region {}".format(region))
                else:
                    current_dict[current_gene] = {'UpDP_NormC':0,'DownDP_NormC':0,'UpDP_Count':0,'DownDP_Count':0,'Loops_Count':1}
                    if 'UpD' in updown:
                        current_dict[current_gene]['UpDP_NormC'] += line_normC
                        current_dict[current_gene]['UpDP_Count'] += line_count
                    elif 'DownD' in updown:
                        current_dict[current_gene]['DownDP_NormC'] += line_normC
                        current_dict[current_gene]['DownDP_Count'] += line_count
                    else:
                        print("Abnormal region {}".format(region))

        df_dict = {'Gene':[],'UpDP_NormC':[],'DownDP_NormC':[],'UpDP_Count':[],'DownDP_Count':[],'Loops_Count':[]}
        for gene in current_dict:
            df_dict['Gene'].append(gene)
            df_dict['UpDP_NormC'].append(current_dict[gene]['UpDP_NormC'])
            df_dict['DownDP_NormC'].append(current_dict[gene]['DownDP_NormC'])
            df_dict['UpDP_Count'].append(current_dict[gene]['UpDP_Count'])
            df_dict['DownDP_Count'].append(current_dict[gene]['DownDP_Count'])
            df_dict['Loops_Count'].append(current_dict[gene]['Loops_Count'])
        current_df = pd.DataFrame(df_dict)
        current_df['Total_NormC'] = current_df['UpDP_NormC'] + current_df['DownDP_NormC']
        current_df['Total_Count'] = current_df['UpDP_Count'] + current_df['DownDP_Count']
        tmp_df = current_df[['Gene','Total_NormC']]
        tmp_df2 = current_df[['Gene','Total_Count']]
        tmp_df.columns = ['Gene',t_sample]
        tmp_df2.columns = ['Gene',t_sample]
        all_gc_normc_list.append(tmp_df)
        all_gc_count_list.append(tmp_df2)
        current_df.to_csv("{}/{}_GC_loops.txt".format(outdir2,t_sample),sep='\t',index=False)
        gc_all_info['Sample'].append(t_sample)
        gc_all_info['The number of Genes'].append(len(current_df))
        # plot GC information
        tissue_gc_loops_dict[t_sample] = current_df
    current_f.close()

gene_normc_mat = reduce(lambda left,right: pd.merge(left,right,on=['Gene'],how='outer'), all_gc_normc_list).fillna(0)
gene_count_mat = reduce(lambda left,right: pd.merge(left,right,on=['Gene'],how='outer'), all_gc_count_list).fillna(0)
gene_normc_mat.to_csv("{}/all_gene_normc_mat.txt".format(outdir3),sep='\t',index=False)
gene_count_mat.to_csv("{}/all_gene_count_mat.txt".format(outdir3),sep='\t',index=False)
gc_all_info_df = pd.DataFrame(gc_all_info)
gc_all_info_df.to_excel("{}/Tissue_GC_loops_summary.xlsx".format(outdir3))
plt.figure()
sns.barplot(data=gc_all_info_df,x='Sample',y='The number of Genes',linewidth=2.5, color='blue',edgecolor=".2")
sns.despine()
plt.savefig("{}/All_sample_genes.png".format('/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/'
                                             'BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/analysis_results'))
plt.close()

# Loops Count distribution
fig,axs = plt.subplots(3,4,figsize=(16,12))
for idx,t_sample in enumerate(tissue_samples):
    ax_id = axs[idx//4,idx%4]
    g = sns.histplot(data=tissue_gc_loops_dict[t_sample],x='Loops_Count',ax=ax_id,kde=True)
    g.set(ylim=(0,3500))
    ax_id.title.set_text(t_sample)
plt.tight_layout()
plt.savefig("{}/All_sample_loops_counts.png".format('/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/'
                                                    'BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/analysis_results'))
plt.close()

fig,axs = plt.subplots(3,4,figsize=(16,12))
for idx,t_sample in enumerate(tissue_samples):
    ax_id = axs[idx//4,idx%4]
    g = sns.histplot(data=tissue_gc_loops_dict[t_sample],x='Total_NormC',ax=ax_id,kde=True)
    g.set(ylim=(0,3500))
    ax_id.title.set_text(t_sample)
plt.tight_layout()
plt.savefig("{}/All_sample_total_signal.png".format('/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/'
                                                    'BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/analysis_results'))


################## part4 build median_ratio matrix with R code ####################
# library(DESeq2)
#
# loop_mat_path = '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/Individual_specific_loops/'
# setwd(loop_mat_path)
# data_raw = read.csv('/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/GC_loops/all_gene_count_mat.txt',header=T,sep='\t',row.names = 'Gene')
# meta = data.frame(condition=factor(c(rep("nt",2),rep("pt",5),rep("rt",5))))
# meta$color = c(rep("gray",2),rep("yellow",5),rep("orange",5))
# rownames(meta) = colnames(data_raw)
# meta$condition <- relevel(meta$condition, ref = "nt")
# dds <- DESeqDataSetFromMatrix(countData = data_raw, colData = meta, design = ~ condition)
# dds <- estimateSizeFactors(dds)
# sizeFactors(dds)
# data_deseq2 <- counts(dds, normalized=TRUE)
# write.csv(data_deseq2,'/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/GC_loops/all_gene_medratio_mat.csv')

################## part5 Individual specific analysis ####################
# step1. row zscore for the matrix
# step2. calculate distance and angle between (NT1, NT2) and (RT1, RT1)...
# step3. histogram check the distribution of the sign(angle)_distance and find cutoff for the differentiated genes
# step4. based on the overall distribution fit a model and calculate the p-value
# output: genes list of PT/RT vs NT individual specific
def cal_dist(df, ctrl_col, trt_sample):
    ctrl_df = df[ctrl_col]
    trt_df = df[[trt_sample]*len(ctrl_col)]
    return np.sign(trt_df.mean(axis=1)-ctrl_df.mean(axis=1)) * np.linalg.norm(trt_df.values - ctrl_df.values,axis=1)

def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def cal_x_bigauss(y1,y2,params):
    x1 = params[0] - np.sqrt(2*params[1]**2*np.log(params[2]/y1))
    x2 = params[3] + np.sqrt(2*params[4]**2*np.log(params[5]/y2))
    return x1,x2

analy_path = '/Users/KunFang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/SLG-DLG/fithic/results_res20k_083121/analysis/'
dat_mat = pd.read_csv('{}/GC_loops/all_gene_medratio_mat.csv'.format(analy_path),index_col=0)
dat_mat_col = dat_mat.columns
dat_mat_zscore = dat_mat.apply(zscore,axis=1).apply(pd.Series)
dat_mat_zscore.columns = dat_mat_col

DE_mat_list = []
fig_row = 2
fig_col = 5
fig,ax = plt.subplots(fig_row,fig_col,figsize=(35,8))
for idx,sample in enumerate(pt_samples+rt_samples):
    print(sample)
    new_col = '{}vsNT'.format(sample)
    dat_mat_zscore[new_col] = cal_dist(dat_mat_zscore,nt_samples,sample)
    row_idx = idx//fig_col
    col_idx = idx%fig_col
    g = sns.histplot(data=dat_mat_zscore,x=new_col,ax=ax[row_idx,col_idx])
    hist_heights = [patch.get_height() for patch in g.patches]
    hist_x = [patch.get_x() for patch in g.patches]
    peaks,_=find_peaks(hist_heights,prominence=100)
    A1,A2 = np.array(hist_heights)[peaks][0],np.array(hist_heights)[peaks][-1]
    mu1,mu2 = np.array(hist_x)[peaks][0],np.array(hist_x)[peaks][-1]
    p0 = (mu1,.2,A1,mu2,.2,A2)
    print("Initial Params: {}".format(p0))
    params,cov=curve_fit(bimodal,hist_x,hist_heights,p0)
    sigma=np.sqrt(np.diag(cov))
    hist_new_height = bimodal(np.array(hist_x),*params)
    ax[row_idx,col_idx].plot(np.array(hist_x),hist_new_height,color='red',lw=2,label='model')
    height_cut1,height_cut2 = A1/2,A2/2
    x1,x2 = cal_x_bigauss(height_cut1,height_cut2,params)
    print("DownR dist cut: {}; UpR dist cut: {}".format(x1,x2))
    ax[row_idx,col_idx].axvline(x1,lw=1,ls='--',c='black')
    ax[row_idx,col_idx].axvline(x2,lw=1,ls='--',c='black')
    ax[row_idx,col_idx].title.set_text(new_col)
    ax[row_idx,col_idx].set(xlabel="Distance to NTs")
    bool_val = (dat_mat_zscore[new_col] >= x2).astype(int) - (dat_mat_zscore[new_col] <= x1).astype(int)
    DE_mat_list.append(bool_val)
# ax[2,2].set_visible(False)
# ax[2,3].set_visible(False)
plt.tight_layout()
plt.savefig('{}/analysis_results/individual_specific_hist.png'.format(analy_path))
plt.close()
DE_mat = pd.concat(DE_mat_list,axis=1)
DE_mat.columns = pt_samples+rt_samples
DE_mat[DE_mat==0]=0
UpR_mat = DE_mat.copy()
DownR_mat = DE_mat.copy()
UpR_mat[UpR_mat==-1] =0
DownR_mat[DownR_mat==1] = 0
DownR_mat[DownR_mat==-1] = 1
UpR_mat.to_csv('{}/individual_specific_loops/IndSp_Gained.txt'.format(analy_path),sep='\t')
DownR_mat.to_csv('{}/individual_specific_loops/IndSp_Lost.txt'.format(analy_path),sep='\t')

DE_summary = pd.DataFrame(pd.concat([UpR_mat.sum().T,DownR_mat.sum().T]))
DE_summary['Type'] = np.array(['Gained']*10+['Lost']*10)
DE_summary.columns = ['Counts','Type']
DE_summary['Sample'] = DE_summary.index.values
# Fig.3C
plt.figure()
g= sns.barplot(data=DE_summary,x='Sample',y='Counts',hue='Type',linewidth=2.5,edgecolor=".1")
g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
sns.despine()
plt.savefig("{}/analysis_results/TTvsNT_individual_specific.png".format(analy_path),bbox_inches='tight')
plt.close()

# plot heatmap
fig,ax = plt.subplots(5,2,figsize=(10,35))
for idx,sample in enumerate(pt_samples):
    row_idx = idx
    col_idx = 0
    current_genes = DE_mat.loc[DE_mat[sample]!=-0,sample].sort_values(ascending=False).index.values
    current_df = dat_mat.loc[current_genes,["NT1","NT2",sample]]
    current_df_zscore = current_df.apply(zscore,axis=1).apply(pd.Series)
    current_df_zscore.columns = ["NT1","NT2",sample]
    sns.heatmap(current_df_zscore,cmap='coolwarm',ax=ax[row_idx,col_idx],yticklabels=False)
for idx,sample in enumerate(rt_samples):
    row_idx = idx
    col_idx = 1
    current_genes = DE_mat.loc[DE_mat[sample]!=-0,sample].sort_values(ascending=False).index.values
    current_df = dat_mat.loc[current_genes,["NT1","NT2",sample]]
    current_df_zscore = current_df.apply(zscore,axis=1).apply(pd.Series)
    current_df_zscore.columns = ["NT1","NT2",sample]
    sns.heatmap(current_df_zscore,cmap='coolwarm',ax=ax[row_idx,col_idx],yticklabels=False)
plt.tight_layout()
# Fig.3C
plt.savefig('{}/analysis_results/individual_specific_heatmap.png'.format(analy_path))
plt.close()

# RTvsNT - (PTvsNT)
gain_df = pd.read_table('{}/individual_specific_loops/IndSp_Gained.txt'.format(analy_path), sep='\t', header=0, index_col=0)
loss_df = pd.read_table('{}/individual_specific_loops/IndSp_Lost.txt'.format(analy_path), sep='\t', header=0, index_col=0)
signal_df = pd.read_csv('{}/GC_loops/all_gene_medratio_mat.csv'.format(analy_path), index_col=0)
# find PT's common genes with criterion at least 2 samples have them
common_cut = 2
gain_df['PTs'] = (gain_df[pt_samples].sum(axis=1) >= common_cut).astype(int)
loss_df['PTs'] = (loss_df[pt_samples].sum(axis=1) >= common_cut).astype(int)

# simple RT_ind - PTs
bar_dict = {'Sample': [], 'Type': [], 'The Number of Genes': []}
vsPT_cols = []
for sample in rt_samples:
    current_col = "{}vsPTs".format(sample)
    vsPT_cols.append(current_col)
    gain_df[current_col] = gain_df[sample] - gain_df["PTs"]
    gain_df.loc[gain_df[current_col] <= 0, current_col] = 0
    bar_dict['Sample'].append(sample)
    bar_dict['Type'].append('Gained')
    bar_dict['The Number of Genes'].append(gain_df[current_col].sum())
    loss_df[current_col] = loss_df[sample] - loss_df["PTs"]
    loss_df.loc[loss_df[current_col] <= 0, current_col] = 0
    bar_dict['Sample'].append(sample)
    bar_dict['Type'].append('Lost')
    bar_dict['The Number of Genes'].append(loss_df[current_col].sum())

# visual
bar_df = pd.DataFrame(bar_dict)
# Fig.3D
bar_df.to_csv("{}/Individual_specific_loops/RTIndvsPTs_summary.csv".format(analy_path))
plt.figure()
g = sns.barplot(data=bar_df, x='Sample', y='The Number of Genes', hue='Type', linewidth=2.5, edgecolor=".1")
g.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("{}/analysis_results/RTIndvsPTs_individual_specific.png".format(analy_path), bbox_inches='tight')
plt.close()

###### flowerplot is plotted by R with the following code
# library(plotrix)
# flower_plot < - function(sample, value, start, a, b,
#                          ellipse_col,
#                          circle_col,
#                          circle_text_cex=1
#                          )
# {
#     par(bty="n", ann=F, xaxt="n", yaxt="n", mar=c(1, 1, 1, 1))
# plot(c(0, 20), c(0, 20), type="n")
# n < - length(sample)
# deg < - 360 / n
# res < - lapply(1: n, function(t)
# {
#     draw.ellipse(x=10 + (2 + b[t]) * cos((start + deg * (t - 1)) * pi / 180),
#                  y=10 + (2 + b[t]) * sin((start + deg * (t - 1)) * pi / 180),
#                  col=ellipse_col[t],
#                  border="black",
#                  a=a, b=2 + b[t], angle=deg * (t - 1))
# text(x=10 + (3 + b[t]) * cos((start + deg * (t - 1)) * pi / 180),
#      y=10 + + (3 + b[t]) * sin((start + deg * (t - 1)) * pi / 180),
#      value[t]
#      )
# })
# draw.circle(x=10, y=10, r=1.0, col=circle_col, border=circle_col)
# }
#
# df = data.frame("Samples" = c("RT1_Gained", "RT1_Lost", "RT2_Gained", "RT2_Lost",
#                               "RT3_Gained", "RT3_Lost", "RT4_Gained", "RT4_Lost", "RT5_Gained", "RT5_Lost"))
# df$Values = c(1581, 559, 1189, 500, 1233, 563, 1668, 607, 1428, 427)
# # col2rgb("#3274A1") col2rgb("#E1812C")
# myblue = rgb(50, 116, 161, max=255, alpha=220, names="blue50")
# myorange = rgb(225, 129, 44, max=255, alpha=220, names="orange50")
# mycircle = rgb(255, 255, 255, max=255, alpha=200, names="circle")
# df$epi_color = c(myblue, myorange, myblue, myorange,
#                  myblue, myorange, myblue, myorange, myblue, myorange)
#
# df$epi_size = 1 * df$Values / as.integer(max(df$Values))
# flower_plot(df$Samples, df$Values, 90, 0.5, df$epi_size, df$epi_color, mycircle)

# extract loss gain signal for RT Individual
gain_signal_list = []
loss_signal_list = []
writer = pd.ExcelWriter('{}/individual_specific_loops/RTsIndvsPT.xlsx'.format(analy_path))
for sample in rt_samples:
    current_col = "{}vsPTs".format(sample)
    current_gain_signal = signal_df.loc[gain_df[gain_df[current_col] == 1].index.values, pt_samples + rt_samples]
    current_gain_df = gain_df.loc[gain_df[current_col] == 1, pt_samples + rt_samples]
    current_gain_df.to_excel(writer, sheet_name="{}_Gained".format(sample))
    gain_signal_list.append(current_gain_signal)
    current_loss_signal = signal_df.loc[loss_df[loss_df[current_col] == 1].index.values, pt_samples + rt_samples]
    current_loss_df = loss_df.loc[loss_df[current_col] == 1, pt_samples + rt_samples]
    current_loss_df.to_excel(writer, sheet_name="{}_Lost".format(sample))
    loss_signal_list.append(current_loss_signal)
writer.save()

dl_df = pd.concat(gain_signal_list + loss_signal_list)
plt.figure()
dl_df_zscore = dl_df.apply(zscore, axis=1).apply(pd.Series)
dl_df_zscore.columns = pt_samples + rt_samples
sns.heatmap(dl_df_zscore, cmap='coolwarm', yticklabels=False)
plt.tight_layout()
plt.savefig('{}/analysis_results/RTIndvsPTs_specific_heatmap.png'.format(analy_path))
plt.close()

# RTs Ind's common
RTsInd_gain = gain_df.loc[gain_df[vsPT_cols].sum(axis=1) >= 1, vsPT_cols]
RTsInd_loss = loss_df.loc[loss_df[vsPT_cols].sum(axis=1) >= 1, vsPT_cols]

RTsInd_gain['Commons'] = RTsInd_gain.sum(axis=1)
RTsInd_loss['Commons'] = RTsInd_loss.sum(axis=1)

# Fig.3E visualized in excel
writer_2 = pd.ExcelWriter('{}/individual_specific_loops/RTsIndvsPT_commons.xlsx'.format(analy_path))
common_dict = {'Common Type': [], "Gained": [], "Lost": []}
gain_common_dict = Counter(RTsInd_gain['Commons'])
loss_common_dict = Counter(RTsInd_loss['Commons'])
for key in sorted(gain_common_dict.keys()):
    current_type = ">={}-RTs".format(key)
    common_dict['Common Type'].append(current_type)
    # current_gain_signal = signal_df.loc[RTsInd_gain[RTsInd_gain['Commons'] >= key].index.values,pt_samples+rt_samples]
    current_gain_signal = RTsInd_gain[RTsInd_gain['Commons'] >= key]
    current_gain_signal.to_excel(writer_2, sheet_name=">={}_RTs_Gained".format(key))
    common_dict['Gained'].append(len(current_gain_signal))

for key in sorted(loss_common_dict.keys()):
    current_type = ">={}-RTs".format(key)
    # current_loss_signal = signal_df.loc[RTsInd_loss[RTsInd_loss['Commons'] >= key].index.values,pt_samples+rt_samples]
    current_loss_signal = RTsInd_loss[RTsInd_loss['Commons'] >= key]
    current_loss_signal.to_excel(writer_2, sheet_name=">={}_RTs_Lost".format(key))
    if current_type in common_dict['Common Type']:
        common_dict['Lost'].append(len(current_loss_signal))
    else:
        common_dict['Common Type'].append(current_type)
        common_dict['Gained'].append(0)
        common_dict['Lost'].append(len(current_loss_signal))
writer_2.save()