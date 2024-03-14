import sys
import time
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from collections import Counter, defaultdict
from scipy import stats
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from pathlib import Path
from scipy.stats import pearsonr
import logomaker
from functools import reduce

# based on TopDom result
# v2 implements a new split-fuse method
# v3 reduce to C tads array and combine cluster method
# v4 change terminology: MV: inter-tumor moderately variable; SV: inter-tumor significant variable; C: conserved; HS: individual high-specific; LS:individual low-specific
# v5 re-implement soft border by new method
# unfinished! v6 assign conservation score and p-value/fdr
# v7 output class table and the summary for outputs. unfinished! identify individual specific TADs array with distance
class TadsArray(object):
    def __init__(self, tads_files, samples, binsize, winsize):
        # it is important that tads_files and samples in the same order
        # flexibility represents how many bins that allow to consider two ids are closed. e.g. if flexibilty = 1, then 17 and 18
        # are considered as closed but not for 16 and 18
        self.tads_files = tads_files
        self.samples = samples
        self.binsize = binsize
        self.winsize = winsize
    def loading_files(self):
        print("Loading files")
        data_tads = {sample:0 for sample in samples}
        tad_ids = []
        for idx,file in enumerate(self.tads_files):
            current_df = pd.read_csv(file,sep='\t')
            # drop row that has same from.coord and to.coord, annoying topdom bugs
            current_df = current_df.drop(current_df.loc[current_df.loc[:,'from.coord']==current_df.loc[:,'to.coord'],:].index.values)
            current_df.reset_index(inplace=True)
            # make sure same data type of "to.id"
            current_df["to.id"] = pd.to_numeric(current_df["to.id"], downcast="float")
            combined_mark = list(current_df.loc[:,'chr']+'_'+current_df.loc[:,'to.id'].astype(str))
            tad_ids.append(combined_mark)
            current_df.loc[:,'mark'] = combined_mark
            # assign the vector value for each region
            current_df.loc[:,'vec'] = current_df.loc[:,'size']/self.binsize
            current_df.loc[current_df.loc[:,'tag']=='gap','vec'] = -current_df.loc[current_df.loc[:,'tag']=='gap','vec']
            current_df.loc[current_df.loc[:,'tag']=='boundary','vec'] = -current_df.loc[current_df.loc[:,'tag']=='boundary','vec']
            data_tads[self.samples[idx]] = current_df
        return data_tads, tad_ids
    def tads_basic_stats(self, data_tads):
        data_stat_tads = {'Sample':[],'domain':[],'boundary':[],'gap':[],'mean.tads.size':[]}
        data_dist_tads = {'Sample':[],'size':[]}
        for sample in self.samples:
            current_df = data_tads[sample]
            # get some statistics
            current_size = current_df.loc[current_df.loc[:,'tag']=='domain','size']
            data_dist_tads['Sample'] += [sample] * len(current_size)
            data_dist_tads['size'] += list(current_size)
            data_stat_tads['Sample'].append(sample)
            data_stat_tads['mean.tads.size'].append(current_size.mean())
            for (k,v) in Counter(current_df['tag']).items():
                data_stat_tads[k].append(v)
        data_dist_tads_df = pd.DataFrame.from_dict(data_dist_tads)
        data_stat_tads_df = pd.DataFrame.from_dict(data_stat_tads)
        return data_dist_tads_df,data_stat_tads_df
    def _split_by_common(self, data_tads, common_coor):
        # return data_tads_split (data_tads_split['NT1']['chr1'][0])
        total_num_tads_array = {sample:0 for sample in samples}
        data_tads_split ={sample:{} for sample in samples}
        for sample in self.samples:
            current_df = data_tads[sample]
            # group by chromosomes
            current_gb = current_df.groupby('chr')
            for chrom in current_gb.groups:
                current_chr_df = current_df.loc[current_gb.groups[chrom],:]
                # +1 is essential for the correct splits
                idx_list = current_chr_df[current_chr_df.loc[:,'mark'].isin(common_coor)].index.values + 1
                idx_list_sort = sorted(idx_list)
                idx_list_sort = idx_list_sort - current_chr_df.index.values[0]
                # avoid generating empty array
                if idx_list_sort[-1] == len(current_chr_df):
                    idx_list_sort = idx_list_sort[:-1]
                # merge second last tads array and last tads array for each chromosome (empirical)
                idx_list_sort = idx_list_sort[:-1]
                data_tads_split[sample][chrom] = np.split(current_chr_df,idx_list_sort)
                total_num_tads_array[sample] += len(data_tads_split[sample][chrom])
                # make sure each chromosome in data_tads_split[sample] has same number of tads array
        test_total = list(total_num_tads_array.values())[0]
        for sample in self.samples:
            if total_num_tads_array[sample] != test_total:
                print("Encounter different number of tads array, exit")
                exit(1)
        return data_tads_split
    def _tads_matrix(self,data_tads_split,chroms):
        vector_matrices = {}
        # vector_annot['tads_id'] = [chr,start,end,{'NT1':"d12,b10,g1,d20"}]
        vector_annot = defaultdict(list)
        tads_id_list = []
        for chrom in chroms:
            for i in range(len(data_tads_split[self.samples[0]][chrom])):
                vec_list = []
                try:
                    tads_id = chrom+'_'+str(i)
                    tads_id_list.append(tads_id)
                    # chr
                    vector_annot[tads_id].append(data_tads_split[self.samples[0]][chrom][i].iloc[0,1])
                    # start
                    vector_annot[tads_id].append(data_tads_split[self.samples[0]][chrom][i].iloc[0,3])
                    # end
                    vector_annot[tads_id].append(data_tads_split[self.samples[0]][chrom][i].iloc[-1,5])
                    current_sample_annot = {}
                    for sample in self.samples:
                        current_annot_list = list(data_tads_split[sample][chrom][i].loc[:,'tag'].str[0] +
                                                  (data_tads_split[sample][chrom][i].loc[:,'size']/binsize).astype(str))
                        current_sample_annot[sample]= ','.join(current_annot_list)
                        current_vec = data_tads_split[sample][chrom][i].loc[:,'vec']
                        current_vec.index.values[:] = np.arange(len(current_vec))
                        vec_list.append(data_tads_split[sample][chrom][i].loc[:,'vec'])
                    vector_annot[tads_id].append(current_sample_annot)
                    current_vec_matrix = pd.concat(vec_list,axis=1)
                    current_vec_matrix.columns = self.samples
                    vector_matrices[chrom+'_'+str(i)] = current_vec_matrix.T
                except IndexError:
                    # last element is blank
                    print("Skip blank elements")
                    continue
        return vector_matrices, vector_annot, tads_id_list
    def _matrix2array(self,matrix):
        matrix_fill = matrix.fillna(0)
        matrix_dict = matrix_fill.T.to_dict('list')
        final_dict = {}
        for key in matrix_dict:
            final_dict[key] = [i for i in matrix_dict[key] if i !=0]
        return final_dict
    def _softseg(self,trt_list,ctrl_list):
        # flexibility is the parameter define the degree of the softness of the border
        flexibility = 2
        trt_abs = np.abs(trt_list)
        ctrl_abs = np.abs(ctrl_list)
        trt_cumsum = np.cumsum(trt_abs)
        ctrl_cumsum = np.cumsum(ctrl_abs)
        pile_trt = np.tile(trt_cumsum,(len(ctrl_cumsum),1))
        index_trt = pile_trt - np.asarray(ctrl_cumsum)[:,None]
        index_trt_abs = np.abs(index_trt)
        # find the index of the soft border
        ctrl_idx,trt_idx = np.where(index_trt_abs <= flexibility)
        # print(trt_idx,ctrl_idx)
        # avoid consecutive soft border
        new_ctrl_idx = []
        new_trt_idx = []
        for idx,value in enumerate(ctrl_idx):
            if value not in new_ctrl_idx and trt_idx[idx] not in new_trt_idx:
                new_ctrl_idx.append(value)
                new_trt_idx.append(trt_idx[idx])
            else:
                continue
        # print(new_trt_idx,new_ctrl_idx)
        new_ctrl_idx = np.array(new_ctrl_idx)
        new_trt_idx = np.array(new_trt_idx)
        new_ctrl_idx += 1
        new_trt_idx += 1
        new_ctrl_idx = np.concatenate([[0],new_ctrl_idx])
        new_trt_idx = np.concatenate([[0],new_trt_idx])
        ctrl_tads_num = np.diff(new_ctrl_idx)
        trt_tads_num = np.diff(new_trt_idx)
        ctrl_tads_split_idx = np.cumsum(ctrl_tads_num)
        trt_tads_split_idx = np.cumsum(trt_tads_num)
        return trt_tads_split_idx[:-1],ctrl_tads_split_idx[:-1]
    def new_idx(self,split_idx_dict,minlen_sample):
        final_split_idx = {sample:0 for sample in samples}
        seg_idx = [split_idx_dict[minlen_sample][sample][0] for sample in split_idx_dict[minlen_sample]]
        # find the common seg for all pairwised comparison
        common_seg_idx = sorted(list(set(seg_idx[0]).intersection(*seg_idx)))
        final_split_idx[minlen_sample] = common_seg_idx
        for sample in split_idx_dict[minlen_sample]:
            target_seg_idx = split_idx_dict[minlen_sample][sample][1]
            current_seg_idx =  split_idx_dict[minlen_sample][sample][0]
            new_idx = [current_seg_idx.tolist().index(idx) for idx in common_seg_idx]
            final_split_idx[sample] = [target_seg_idx[idx] for idx in new_idx]
        return final_split_idx
    def find_subcommon(self,input_dict,samples):
        split_idx_dict = {sample:{} for sample in samples}
        min_len = len(min(input_dict.items(),key=lambda x:len(x[1]))[1])
        minlen_samples = [sample for sample in samples if len(input_dict[sample])==min_len]
        final_seg = {}
        # store all pairwise segmentation index
        for trt_ctrl_sample in itertools.combinations(samples,2):
            trt_sample = trt_ctrl_sample[0]
            ctrl_sample = trt_ctrl_sample[1]
            trt_split_idx, ctrl_split_idx = self._softseg(input_dict[trt_sample],input_dict[ctrl_sample])
            if len(trt_split_idx)==0 or len(ctrl_split_idx) ==0:
                return {sample:[np.array(input_dict[sample])] for sample in input_dict}
            else:
                split_idx_dict[trt_sample][ctrl_sample] = [trt_split_idx,ctrl_split_idx]
                split_idx_dict[ctrl_sample][trt_sample] = [ctrl_split_idx,trt_split_idx]
        if len(minlen_samples) ==1:
            # find the shortest segment
            final_split_idx = self.new_idx(split_idx_dict,minlen_samples[0])
        else:
            idx_len = []
            tmp_store = {}
            for idx,minlen_sample in enumerate(minlen_samples):
                split_idx = self.new_idx(split_idx_dict,minlen_sample)
                idx_len.append(len(split_idx[minlen_sample]))
                tmp_store[idx] = split_idx
            minidx = int(np.argmin(idx_len))
            final_split_idx = tmp_store[minidx]
        for sample in samples:
            final_seg[sample] = np.split(input_dict[sample],final_split_idx[sample])
        return final_seg
    def tads_array_chunk(self, data_tads, tad_ids):
        # find common coordinates id for all samples
        common_coor = list(set(tad_ids[0]).intersection(*tad_ids))
        print("The number of common coordinates(0 flexiblity): {}".format(len(common_coor)))
        chroms = [i for i in data_tads[self.samples[0]].groupby('chr').groups]
        print("First round segmentation")
        data_tads_split = self._split_by_common(data_tads,common_coor)
        # construct tad vector matrix, each row represent the tads array for a sample (length might be varied)
        print("Construct vector matrices")
        vector_matrices, vector_annot, tads_id_list = self._tads_matrix(data_tads_split,chroms)
        print("Second round segmentation")
        tads_sub_array_dict = {}
        tads_sub_annot_dict = {}
        tads_sub_id_list =[]
        for array_mat_id in tqdm(vector_matrices):
            current_df = vector_matrices[array_mat_id]
            current_dict = self._matrix2array(current_df)
            current_split_seg = self.find_subcommon(current_dict,samples)
            chunk_num = len(current_split_seg[self.samples[0]])
            for idx in range(chunk_num):
                tads_sub_id = "{}_{}".format(array_mat_id,idx)
                tads_sub_id_list.append(tads_sub_id)
                tads_sub_array_dict[tads_sub_id] = {sample:current_split_seg[sample][idx] for sample in self.samples}
                if idx != 0:
                    tads_sub_annot_dict[tads_sub_id] = {sample:[vector_annot[array_mat_id][0],
                                                                vector_annot[array_mat_id][1] + self.binsize * np.sum(np.abs(np.concatenate(current_split_seg[sample][:idx]))),
                                                                vector_annot[array_mat_id][1] + self.binsize * np.sum(np.abs(np.concatenate(current_split_seg[sample][:idx+1])))]
                                                        for sample in self.samples}
                else:
                    tads_sub_annot_dict[tads_sub_id] = {sample:[vector_annot[array_mat_id][0],
                                                                vector_annot[array_mat_id][1],
                                                                vector_annot[array_mat_id][1] + self.binsize * np.sum(np.abs(current_split_seg[sample][0]))]
                                                        for sample in self.samples}
        return tads_sub_array_dict, tads_sub_annot_dict, tads_sub_id_list

class TadsMatrix(object):
    def __init__(self,tads_array_dict,samples,tads_id_list):
        self.tads_array_dict = tads_array_dict
        self.samples = samples
        self.tads_id_list = tads_id_list
        self.flexibility =2
    def _neo_del_score(self,input_dict):
        # neo positive; del negative
        sum_series = pd.DataFrame.from_dict(input_dict,orient='index').sum(axis=1)
        # row ctrl, column trt, mat_row_expand is dup sum_series in rows
        mat_row_expand = np.tile(sum_series.tolist(),(len(self.samples),1))
        mat_col_expand = mat_row_expand.T
        # find diff between any two samples
        nd_mat = pd.DataFrame((mat_row_expand - mat_col_expand)/2,index=self.samples,columns=self.samples)
        # if neo_del less than flexibility, make it 0. Might change to some relative conserve information in the future
        nd_mat[nd_mat.abs()<=self.flexibility] = 0
        return nd_mat
    def _split_fuse_score(self,trt_list,ctrl_list):
        if len(trt_list) > len(ctrl_list):
            split_fuse_score = sum(np.sort(np.abs(trt_list))[::-1][len(ctrl_list):])
        elif len(trt_list) < len(ctrl_list):
            split_fuse_score = -sum(np.sort(np.abs(ctrl_list))[::-1][len(trt_list):])
        else:
            split_fuse_score = 0.0
        if np.abs(split_fuse_score) <= self.flexibility:
            split_fuse_score = 0.0
        return split_fuse_score
    def matrix_build(self,input_dict):
        # return 'skew-symmetric' matrix with shape sample x sample. The 'rows' used as the ctrl, for exmaple: the element in ['NT1','NT2'] means
        # use 'NT1' as ctrl, the tads array change in 'NT2'; each element in the matrix is in two dimensional [neo/del,split/fuse].
        sf_df_raw = pd.DataFrame(0,index=self.samples,columns=self.samples)
        nd_df_raw = self._neo_del_score(input_dict)
        nd_df = nd_df_raw.abs()
        for sample in itertools.combinations(self.samples,2):
            trt_sample = sample[0]
            ctrl_sample = sample[1]
            trt_array = input_dict[trt_sample]
            ctrl_array = input_dict[ctrl_sample]
            split_fuse_score = self._split_fuse_score(trt_array,ctrl_array)
            sf_df_raw.loc[trt_sample,ctrl_sample] = -split_fuse_score
        sf_df_raw = pd.DataFrame(np.triu(sf_df_raw) - np.triu(sf_df_raw,1).T,index=self.samples,columns=self.samples)
        sf_df = sf_df_raw.abs()
        detail_info = {'nd_raw':nd_df_raw,'sf_raw':sf_df_raw}
        return nd_df, sf_df, detail_info
    def fit(self,comparison_dict):
        # out_matrix_dict = {}
        # stats_dict is feature dict
        stats_dict = {'tads_id':self.tads_id_list}
        for compare_type in comparison_dict:
            stats_type_nd = compare_type + '_nd'
            stats_type_sf = compare_type + '_sf'
            stats_dict[stats_type_nd] = []
            stats_dict[stats_type_sf] = []
        print("Building featured matrix")
        total_info = {'total_detail':{},'total_nd':{},'total_sf':{}}
        for tads_array_id in tqdm(self.tads_id_list):
            current_dict = self.tads_array_dict[tads_array_id]
            current_nd, current_sf, detail_info = self.matrix_build(current_dict)
            total_info['total_detail'][tads_array_id] = detail_info
            total_info['total_nd'][tads_array_id] = current_nd
            total_info['total_sf'][tads_array_id] = current_sf
            # fill the stats_dict
            for compare_type in comparison_dict:
                stats_type_nd = compare_type + '_nd'
                stats_type_sf = compare_type + '_sf'
                trt_samples = comparison_dict[compare_type][0]
                ctrl_samples = comparison_dict[compare_type][1]
                # inter and intra tissue type comparison
                current_part_nd = current_nd.loc[ctrl_samples,trt_samples]
                current_part_sf = current_sf.loc[ctrl_samples,trt_samples]
                if trt_samples != ctrl_samples:
                    # inter
                    stats_dict[stats_type_nd].append(np.round(current_part_nd.stack().mean(),3))
                    stats_dict[stats_type_sf].append(np.round(current_part_sf.stack().mean(),3))
                else:
                    tri_nd = np.array(current_part_nd)[np.triu_indices(len(trt_samples),k=1)]
                    tri_sf = np.array(current_part_sf)[np.triu_indices(len(trt_samples),k=1)]
                    stats_dict[stats_type_nd].append(np.round(np.mean(tri_nd),3))
                    stats_dict[stats_type_sf].append(np.round(np.mean(tri_sf),3))
        stats_df = pd.DataFrame(stats_dict)
        return total_info, stats_df

def permutation_test(df, gcomp, lvl1_coef=0.7, lvl2_coef=0.3, n_permutations=1000):
    """Performs a Monte Carlo permutation test on a row of data and returns the p-value."""
    data = df.copy()
    # Set the seed for the random number generator
    np.random.seed(42)
    lvl1 = gcomp.split('.')[2]
    lvl2 = gcomp.split('.')[0]
    # Q99 of TT.vs.NT_nd
    Qnd_group = data['{}_nd'.format(gcomp)].quantile(0.99)
    Qsf_group = data['{}_sf'.format(gcomp)].quantile(0.99)
    # Q99 value of Intra-NT_nd
    Qnd_lvl1 = data['Intra-{}_nd'.format(lvl1)].quantile(0.99)
    Qsf_lvl1 = data['Intra-{}_sf'.format(lvl1)].quantile(0.99)
    Qnd_lvl2 = data['Intra-{}_nd'.format(lvl2)].quantile(0.99)
    Qsf_lvl2 = data['Intra-{}_sf'.format(lvl2)].quantile(0.99)
    # print('Qnd_group:{},Qnd_lvl1:{},Qnd_lvl2:{}'.format(Qnd_group,Qnd_lvl1,Qnd_lvl2))
    # print('Qsf_group:{},Qsf_lvl1:{},Qsf_lvl2:{}'.format(Qsf_group, Qsf_lvl1, Qsf_lvl2))
    # Calculate the observed test statistic
    observed_statistic = (data['{}_nd'.format(gcomp)]/Qnd_group + data['{}_sf'.format(gcomp)]/Qsf_group) - \
                         lvl1_coef * (data['Intra-{}_nd'.format(lvl1)]/Qnd_lvl1 + data['Intra-{}_sf'.format(lvl1)]/Qsf_lvl1) - \
                         lvl2_coef * (data['Intra-{}_nd'.format(lvl2)]/Qnd_lvl2 + data['Intra-{}_sf'.format(lvl2)]/Qsf_lvl2)

    # Generate a distribution of test statistics under the null hypothesis
    permuted_obs = [pd.Series(np.random.permutation(observed_statistic)) for i in range(n_permutations)]
    permuted_obs_df = pd.concat(permuted_obs,axis=1)

    # Calculate the p-value as the proportion of permuted test statistics that are as or more extreme than the observed test statistic
    p_values = (permuted_obs_df >= observed_statistic.values[:, None]).mean(axis=1)
    return p_values

def cal_pval(feature_df, gcomp):
    max_d_number = 0
    final_pval = 0
    # optimize lvl1 and lvl2 coefficient to get the maximum number of dchange
    only_inter_pval = 0
    for i in range(10):
        for j in range(10):
            lvl1_coef = i / 10
            lvl2_coef = j / 10
            p_values = permutation_test(df=feature_df, lvl1_coef=lvl1_coef, lvl2_coef=lvl2_coef,
                                                    gcomp=gcomp)
            sig_num = len(p_values[p_values < 0.05])
            if i == 0 and j == 0:
                only_inter_pval = permutation_test(df=feature_df, lvl1_coef=lvl1_coef, lvl2_coef=lvl2_coef,
                                                    gcomp=gcomp)
            if max_d_number < sig_num:
                max_d_number = sig_num
                final_pval = p_values
            else:
                continue
    return [only_inter_pval,final_pval]

def class_process(input_df,group_level,individual_level1,individual_level2,group_diff_cut,individual_diff_cut):
    # group_level TT.vs.NT_nd
    df = input_df.copy()
    group_cut_up = input_df.stack()[input_df.stack()!=0].quantile(group_diff_cut)
    indi_cut_up = input_df.stack()[input_df.stack()!=0].quantile(individual_diff_cut)
    group_cut_down = input_df.stack()[input_df.stack()!=0].quantile(0.1)
    indi_cut_down = group_cut_down
    print("up_cut_val: {}; down_up_val: {}".format(group_cut_up,group_cut_down))
    df.loc[df.loc[:,group_level]>group_cut_down,group_level+'_class'] = 'medium'
    df.loc[df.loc[:,individual_level1]>indi_cut_down,individual_level1+'_class'] = 'medium'
    df.loc[df.loc[:,individual_level2]>indi_cut_down,individual_level2+'_class'] = 'medium'
    df.loc[df.loc[:,group_level]>group_cut_up,group_level+'_class'] = 'high'
    df.loc[df.loc[:,individual_level1]>indi_cut_up,individual_level1+'_class'] = 'high'
    df.loc[df.loc[:,individual_level2]>indi_cut_up,individual_level2+'_class'] = 'high'
    df = df.fillna('low')
    return df, [group_cut_up,group_cut_down]

def combine_nd_sf_annot(tads_id,nd_marked_df,sf_marked_df,group_level,individual_level1,individual_level2):
    # group_level TT.vs.NT
    # diff_cut 0.7 means we define top 30% as the drastically change
    combine_dict = {'high_high':'high','high_medium':'high','high_low':'high','medium_high':'high','medium_medium':'medium','medium_low':'medium',
                    'low_high':'high','low_medium':'medium','low_low':'low'}
    types_dict = {'high_high':'Mixed','high_medium':'Mixed','high_low':'ND','medium_high':'Mixed','medium_medium':'NotAvail','medium_low':'NotAvail',
                  'low_high':'SF','low_medium':'NotAvail','low_low':'NotAvail'}
    final_df = pd.DataFrame(tads_id)
    final_df[group_level] = (nd_marked_df['{}_nd_class'.format(group_level)] + '_' +
                             sf_marked_df['{}_sf_class'.format(group_level)]).map(combine_dict)
    final_df[individual_level1] = (nd_marked_df['{}_nd_class'.format(individual_level1)] + '_' +
                                   sf_marked_df['{}_sf_class'.format(individual_level1)]).map(combine_dict)
    final_df[individual_level2] = (nd_marked_df['{}_nd_class'.format(individual_level2)] + '_' +
                                   sf_marked_df['{}_sf_class'.format(individual_level2)]).map(combine_dict)
    final_df['GroupChange'] = 'MV'
    final_df.loc[final_df[group_level]=='high','GroupChange'] = 'SV'
    final_df.loc[(final_df[group_level]=='low') & (final_df[individual_level1]!='high') & (final_df[individual_level2]!='high'),'GroupChange'] = 'C'
    # add ND, SF and ND_SF_Mixed
    final_df['ChangeTypes'] = (nd_marked_df['{}_nd_class'.format(group_level)] + '_' +
                               sf_marked_df['{}_sf_class'.format(group_level)]).map(types_dict)
    # add Individual-high, Individual-M & C
    final_df['IndividualChange'] = 'NotAvail'
    final_df.loc[(final_df[group_level] == 'high') & (final_df[individual_level2] == 'high'),'IndividualChange'] = 'HS'
    final_df.loc[
        (final_df[group_level] == 'high') & (final_df[individual_level2] != 'high'), 'IndividualChange'] = 'LS'
    # final_df.loc[(final_df[group_level]=='low') & (final_df[individual_level1]='low') & (final_df[individual_level2]=='low'),'bio_mark'] = 'C'
    # final_df.loc[(final_df[group_level]=='low') & (final_df[individual_level1]=='meidum') & (final_df[individual_level2]=='low'),'bio_mark'] = 'C'
    return final_df

def coor_annot(features_df_marked, tads_sub_id_list,tads_sub_annot_dict):
    cols = list(features_df_marked.columns)
    tads_subarray_start = []
    tads_subarray_end = []
    tads_subarray_chr = []
    for tads_sub_id in tads_sub_id_list:
        # to avoid overlapping, pick the largest start and smallest end coodinates
        min_start = min(tads_sub_annot_dict[tads_sub_id].items(),key=lambda x:x[1][1])[1][1]
        max_end = max(tads_sub_annot_dict[tads_sub_id].items(),key=lambda x:x[1][2])[1][2]
        tads_subarray_start.append(int(min_start))
        tads_subarray_end.append(int(max_end))
        tads_subarray_chr.append(tads_sub_annot_dict[tads_sub_id][samples[0]][0])
    features_df_marked.loc[:,'start(min)'] = tads_subarray_start
    features_df_marked.loc[:,'end(max)'] = tads_subarray_end
    features_df_marked.loc[:,'chr'] = tads_subarray_chr
    cols = ['chr','start(min)','end(max)']+cols
    features_df_marked = features_df_marked[cols]
    return features_df_marked

def custom_discrete_cmap(numcate,colors,interval):
    newcolors = plt.get_cmap('viridis',numcate).colors
    # assign new colors of interest
    for i in range(0,numcate,interval):
        newcolors[i:i+interval, :] = matplotlib.colors.to_rgba(colors[i])
    # create the customized color map
    cmap = matplotlib.colors.ListedColormap(newcolors)
    return cmap

def plot_heatmap(df,cmap,title,outdir,outname):
    plt.figure()
    ax = sns.heatmap(df,yticklabels=False,cmap=cmap,cbar=False)
    ax.set_title(title)
    ax.vlines([1, 2], colors='black',*ax.get_ylim())
    plt.savefig('{}/{}_{}.png'.format(outdir,outname,title))
    plt.close()

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

def get_tri_heatmap(df, loc='lower', vmin=0.5, vmax=1, fontsize=16, show_text=False, output="cooc_matrix.png"):
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
    sns_plot = sns.heatmap(df_tri, cmap=cmap, vmin=vmin, vmax=vmax, center=(vmin+vmax)/2, annot=df,
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

def summarize_changes(array_dict,changetype,ctrl_samples,trt_samples):
    summarize_df = pd.DataFrame(0, index=trt_samples, columns=['N','D','S','F'])
    for cs in ctrl_samples:
        cs_array = array_dict[cs]
        for ts in trt_samples:
            ts_array = array_dict[ts]
            if changetype == 'ND':
                if len(ts_array[ts_array<0]) > len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'D'] += 1
                elif len(ts_array[ts_array<0]) < len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'N'] += 1
                else:
                    summarize_df.loc[ts, 'N'] += 0
            elif changetype == 'SF':
                if len(ts_array) > len(cs_array):
                    summarize_df.loc[ts, 'S'] += 1
                elif len(ts_array) < len(cs_array):
                    summarize_df.loc[ts, 'F'] += 1
                else:
                    summarize_df.loc[ts, 'F'] += 0
            elif changetype == 'Mixed':
                if len(ts_array[ts_array<0]) > len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'D'] += 1
                elif len(ts_array[ts_array<0]) < len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'N'] += 1
                else:
                    summarize_df.loc[ts, 'N'] += 0
                if len(ts_array) > len(cs_array):
                    summarize_df.loc[ts, 'S'] += 1
                elif len(ts_array) < len(cs_array):
                    summarize_df.loc[ts, 'F'] += 1
                else:
                    summarize_df.loc[ts, 'F'] += 0
    return summarize_df

def summarize_change_all(dchange_df,tads_sub_array_dict,ctrl_samples,trt_samples):
    # return 1.individual_sample_df; 2. as group_df; 3. detailed change for each row
    summarize_list_dict = {'HS':{'ND':[],'SF':[],'Mixed':[]},
                           'LS':{'ND':[],'SF':[],'Mixed':[]}}
    summarize_dict = {}
    sample_changes = {sample:[] for sample in trt_samples}
    sample_changes['tads_id'] = []
    for idx,row in dchange_df.iterrows():
        current_summarize = summarize_changes(tads_sub_array_dict[row['tads_id']],row['ChangeTypes'],ctrl_samples,trt_samples)
        sample_changes['tads_id'].append(row['tads_id'])
        for sample in trt_samples:
            sample_changes[sample].append("{}N;{}D;{}S;{}F".format(current_summarize.loc[sample, 'N'],
                                                                  current_summarize.loc[sample, 'D'],
                                                                  current_summarize.loc[sample, 'S'],
                                                                  current_summarize.loc[sample, 'F']
                                                                  ),
                                                 )

        summarize_list_dict[row['IndividualChange']][row['ChangeTypes']].append(current_summarize)
        # detailed change
        changes = '{}N;{}D;{}S;{}F'.format(current_summarize['N'].sum(), current_summarize['D'].sum(),
                                           current_summarize['S'].sum(), current_summarize['F'].sum())
        summarize_dict[row['tads_id']] = changes
    summarize_all_dict = {'HS':{'ND':reduce(lambda x, y: x.add(y, fill_value=0), summarize_list_dict['HS']['ND']),
                                        'SF':reduce(lambda x, y: x.add(y, fill_value=0), summarize_list_dict['HS']['SF']),
                                        'Mixed':reduce(lambda x, y: x.add(y, fill_value=0), summarize_list_dict['HS']['Mixed'])
                                        },
                          'LS': {'ND': reduce(lambda x, y: x.add(y, fill_value=0), summarize_list_dict['LS']['ND']),
                                           'SF': reduce(lambda x, y: x.add(y, fill_value=0), summarize_list_dict['LS']['SF']),
                                            'Mixed': reduce(lambda x, y: x.add(y, fill_value=0),
                                                            summarize_list_dict['LS']['Mixed'])
                                         }
                          }
    return summarize_dict, summarize_all_dict,pd.DataFrame(sample_changes)

def plotlogo(summarize_df,outname):
    # visualize summarize_df
    color_scheme = {
        'N': [0, .5, 0],
        'D': [1, 0, 0],
        'S': [1, .65, 0],
        'F': [0, 0, 1]
    }
    tick_labels = summarize_df.index.values
    summarize_df_norm = summarize_df/summarize_df.sum(axis=1).max()
    summarize_df_norm.reset_index(drop=True,inplace=True)
    logo = logomaker.Logo(summarize_df_norm, font_name='Arial Rounded MT Bold',color_scheme=color_scheme)
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left'], visible=True, bounds=[0, 2])
    logo.style_spines(spines=['bottom'], visible=True)
    logo.ax.set_ylabel(r"Scaled"+"\n"+ r"Frequency",fontsize=12,weight='bold')
    logo.ax.set_ylim([0 , 1])
    logo.ax.set_yticks([0, 1])
    logo.ax.set_xticks([])
    tick_locations = range(len(tick_labels))
    logo.ax.set_xticks(tick_locations)
    logo.ax.set_xticklabels(tick_labels)
    logo.ax.tick_params(axis='both', which='major', labelsize=16, labelcolor='black')
    for label in logo.ax.get_xticklabels() + logo.ax.get_yticklabels():
        label.set_weight("bold")
    plt.savefig(outname,dpi=500)
    plt.close()
    return summarize_df_norm

if __name__ == "__main__":
    normal_samples = ['NT1','NT2']
    primary_samples = ['PT1','PT2','PT3','PT4','PT5']
    recurrent_samples = ['RT1','RT2','RT3','RT4','RT5']
    tumor_samples = primary_samples + recurrent_samples
    samples = normal_samples + primary_samples +recurrent_samples
    binsize = 40000
    winsize = 5
    group_cutoff = 0.7
    individual_cutoff = 0.7
    datdir = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/tads_calling/diff_tads/input_data'
    outdir = '/Users/kfang/Documents/lab/Jin_lab/Labmates/lava/BRCA_tissue_Analysis/tads_calling/diff_tads/win5_40k/GISTA_013124'
    Path(outdir).mkdir(parents=True, exist_ok=True)

    tads_files = ["{}/{}_{}_win{}_total.txt".format(datdir,sample,binsize,winsize) for sample in samples]
    Tads_array = TadsArray(tads_files,samples,binsize,winsize)
    data_tads, tad_ids = Tads_array.loading_files()
    data_dist_tads_df,data_stat_tads_df = Tads_array.tads_basic_stats(data_tads)
    data_stat_tads_df.to_csv('{}/tads_basic_information.csv'.format(outdir),index=False)
    # the format of the key in tads_array_dict: 'chr1_1' means the first tad array in chr1
    # key in sample_df_dict: chr1_1; value: chr/start/end/detail_info
    tads_sub_array_dict, tads_sub_annot_dict, tads_sub_id_list = Tads_array.tads_array_chunk(data_tads,tad_ids)
    # tads_array_size
    tads_sub_array_size = []
    for tads_id in tads_sub_id_list:
        current_size = 0
        for sample in samples:
            current_size += (tads_sub_annot_dict[tads_id][sample][2] - tads_sub_annot_dict[tads_id][sample][1])/1000000
        current_size /= len(samples)
        tads_sub_array_size.append(current_size)

    plt.figure()
    g = sns.displot(data=pd.DataFrame(tads_sub_array_size), x=0, log_scale=10)
    plt.xlabel('TADs array size (Mb)')
    plt.savefig('{}/TADs_sub_array_size_dist.png'.format(outdir))
    plt.tight_layout()
    plt.close()

    Tads_matrix = TadsMatrix(tads_sub_array_dict,samples,tads_sub_id_list)
    # trt type at 0 position and ctrl type at the 1 position
    comparison_dict ={'TT.vs.NT':[tumor_samples,normal_samples],'RT.vs.PT':[recurrent_samples,primary_samples],
                      'Intra-TT':[tumor_samples,tumor_samples],'Intra-NT':[normal_samples,normal_samples],
                      'Intra-PT':[primary_samples,primary_samples],'Intra-RT':[recurrent_samples,recurrent_samples]}
    # total_info with key total_detail, total_sf, total_nd, nd:Neo-Del, sf: Split-Fuse
    total_info, features_df = Tads_matrix.fit(comparison_dict)
    features_df.to_csv("{}/Tads_sub_array_feature_matrix.csv".format(outdir),index=False)
    TTvsNT_type = ['TT.vs.NT','Intra-NT','Intra-TT']
    RTvsPT_type = ['RT.vs.PT','Intra-PT','Intra-RT']
    nd_TTvsNT_type =[comp_type+'_nd' for comp_type in TTvsNT_type]
    sf_TTvsNT_type =[comp_type+'_sf' for comp_type in TTvsNT_type]
    nd_RTvsPT_type = [comp_type+'_nd' for comp_type in RTvsPT_type]
    sf_RTvsPT_type = [comp_type+'_sf' for comp_type in RTvsPT_type]
    features_nd_TTvsNT = features_df.loc[:,nd_TTvsNT_type]
    features_sf_TTvsNT = features_df.loc[:,sf_TTvsNT_type]
    features_nd_RTvsPT = features_df.loc[:,nd_RTvsPT_type]
    features_sf_RTvsPT = features_df.loc[:,sf_RTvsPT_type]
    # calculate pval for SV
    print("Permutation Test")
    TTvsNT_pvals = cal_pval(pd.concat([features_nd_TTvsNT,features_sf_TTvsNT],axis=1),'TT.vs.NT')
    RTvsPT_pvals = cal_pval(pd.concat([features_nd_RTvsPT,features_sf_RTvsPT],axis=1),'RT.vs.PT')

    # give mnemonic
    features_nd_TTvsNT_marked, TTvsNT_nd_cuts = class_process(features_nd_TTvsNT,'TT.vs.NT_nd','Intra-NT_nd','Intra-TT_nd',group_cutoff,individual_cutoff)
    features_sf_TTvsNT_marked, TTvsNT_sf_cuts = class_process(features_sf_TTvsNT,'TT.vs.NT_sf','Intra-NT_sf','Intra-TT_sf',group_cutoff,individual_cutoff)
    features_nd_RTvsPT_marked, RTvsPT_nd_cuts = class_process(features_nd_RTvsPT,'RT.vs.PT_nd','Intra-PT_nd','Intra-RT_nd',group_cutoff,individual_cutoff)
    features_sf_RTvsPT_marked, RTvsPT_sf_cuts = class_process(features_sf_RTvsPT,'RT.vs.PT_sf','Intra-PT_sf','Intra-RT_sf',group_cutoff,individual_cutoff)
    # nd sf overall distribution
    TTNT_nd_data = features_nd_TTvsNT.stack()
    TTNT_sf_data = features_sf_TTvsNT.stack()
    RTPT_nd_data = features_nd_RTvsPT.stack()
    RTPT_sf_data = features_sf_RTvsPT.stack()
    # visualize
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    f,ax =plt.subplots(2,2,sharey=True,sharex=True)
    ax1 = ax[0][0]
    sns.kdeplot(data=np.log2(TTNT_nd_data+0.1), ax=ax1,color='black',fill=True)
    ax1.axvline(np.log2(TTvsNT_nd_cuts[0]+0.1), linewidth=1, color='r',ls='--')
    ax1.axvline(np.log2(TTvsNT_nd_cuts[1]+0.1), linewidth=1, color='g',ls='--')
    ax1.set_ylabel('TT.vs.NT\nDensity')
    ax2 = ax[0][1]
    sns.kdeplot(data=np.log2(TTNT_sf_data+0.1), ax=ax2,color='purple',fill=True)
    ax2.axvline(np.log2(TTvsNT_sf_cuts[0]+0.1), linewidth=1, color='r', ls='--')
    ax2.axvline(np.log2(TTvsNT_sf_cuts[1]+0.1), linewidth=1, color='g', ls ='--')
    ax3 = ax[1][0]
    sns.kdeplot(data=np.log2(RTPT_nd_data+0.1), ax=ax3,color='black',fill=True)
    ax3.axvline(np.log2(TTvsNT_nd_cuts[0]+0.1),  linewidth=1, color='r',ls='--')
    ax3.axvline(np.log2(TTvsNT_nd_cuts[1]+0.1),  linewidth=1, color='g',ls='--')
    ax3.set_ylabel('RT.vs.PT\nDensity')
    ax3.set_xlabel('Log2(Neo-Del Score+0.1)')
    ax4 = ax[1][1]
    sns.kdeplot(data=np.log2(RTPT_sf_data+0.1), ax=ax4,color='purple',fill=True)
    ax4.axvline(np.log2(RTvsPT_sf_cuts[0]+0.1), linewidth=1, color='r', ls='--')
    ax4.axvline(np.log2(RTvsPT_sf_cuts[1]+0.1), linewidth=1, color='g', ls ='--')
    ax4.set_xlabel('Log2(Split-Fuse Score+0.1)')
    plt.tight_layout()
    plt.savefig(f'{outdir}/ND_SF_toCate.png',dpi=300)
    plt.close()

    features_TTvsNT_df_marked = combine_nd_sf_annot(features_df['tads_id'],features_nd_TTvsNT_marked,features_sf_TTvsNT_marked,
                                             'TT.vs.NT','Intra-NT','Intra-TT')
    features_RTvsPT_df_marked = combine_nd_sf_annot(features_df['tads_id'],features_nd_RTvsPT_marked,features_sf_RTvsPT_marked,
                                                    'RT.vs.PT','Intra-PT','Intra-RT')
    # assign coordinate back
    features_TTvsNT_df_marked = coor_annot(features_TTvsNT_df_marked,tads_sub_id_list,tads_sub_annot_dict)
    features_RTvsPT_df_marked = coor_annot(features_RTvsPT_df_marked,tads_sub_id_list,tads_sub_annot_dict)
    # specific change type
    features_TTvsNT_dchange = features_TTvsNT_df_marked.loc[
                                 features_TTvsNT_df_marked['GroupChange'] == 'SV', :]
    features_RTvsPT_dchange = features_RTvsPT_df_marked.loc[
                                 features_RTvsPT_df_marked['GroupChange'] == 'SV', :]
    TTvsNT_summarize_dict, TTvsNT_summarize_list_dict,TTvsNT_sample_changes = summarize_change_all(features_TTvsNT_dchange,
                                                                      tads_sub_array_dict,normal_samples,
                                                                      primary_samples+recurrent_samples)
    # the maximum changecode's number should be less than len(ctrl_samples)*len(trt_samples)
    features_TTvsNT_df_marked['ChangeCodes'] = features_TTvsNT_df_marked['tads_id'].map(TTvsNT_summarize_dict)
    features_TTvsNT_df_marked = pd.merge(features_TTvsNT_df_marked,TTvsNT_sample_changes,on='tads_id',how='outer')
    features_TTvsNT_df_marked = features_TTvsNT_df_marked.fillna('NotAvail')
    RTvsPT_summarize_dict, RTvsPT_summarize_list_dict, RTvsPT_sample_changes = summarize_change_all(features_RTvsPT_dchange,
                                                                      tads_sub_array_dict,primary_samples,recurrent_samples)
    features_RTvsPT_df_marked['ChangeCodes'] = features_RTvsPT_df_marked['tads_id'].map(RTvsPT_summarize_dict)
    features_RTvsPT_df_marked = pd.merge(features_RTvsPT_df_marked, RTvsPT_sample_changes, on='tads_id', how='outer')
    features_RTvsPT_df_marked = features_RTvsPT_df_marked.fillna('NotAvail')
    # plot sample change logo
    TTvsNT_ivh_nd_norm = plotlogo(TTvsNT_summarize_list_dict['HS']['ND'],'{}/TTNT_IVH_ND_logo.png'.format(outdir))
    TTvsNT_ivh_sf_norm = plotlogo(TTvsNT_summarize_list_dict['HS']['SF'],'{}/TTNT_IVH_SF_logo.png'.format(outdir))
    TTvsNT_ivh_mix_norm = plotlogo(TTvsNT_summarize_list_dict['HS']['Mixed'],'{}/TTNT_IVH_Mixed_logo.png'.format(outdir))

    TTvsNT_imc_nd_norm = plotlogo(TTvsNT_summarize_list_dict['LS']['ND'],'{}/TTNT_IMC_ND_logo.png'.format(outdir))
    TTvsNT_imc_sf_norm = plotlogo(TTvsNT_summarize_list_dict['LS']['SF'],'{}/TTNT_IMC_SF_logo.png'.format(outdir))
    TTvsNT_imc_mix_norm = plotlogo(TTvsNT_summarize_list_dict['LS']['Mixed'],'{}/TTNT_IMC_Mixed_logo.png'.format(outdir))

    TTvsNT_nd_norm = plotlogo(TTvsNT_summarize_list_dict['HS']['ND'] + TTvsNT_summarize_list_dict['LS']['ND'],
                              '{}/TTNT_ND_logo.png'.format(outdir))
    TTvsNT_sf_norm = plotlogo(TTvsNT_summarize_list_dict['HS']['SF'] + TTvsNT_summarize_list_dict['LS']['SF'],
                              '{}/TTNT_SF_logo.png'.format(outdir))
    TTvsNT_mix_norm = plotlogo(TTvsNT_summarize_list_dict['HS']['Mixed'] + TTvsNT_summarize_list_dict['LS']['Mixed'],
                               '{}/TTNT_Mixed_logo.png'.format(outdir))

    RTvsPT_ivh_nd_norm = plotlogo(RTvsPT_summarize_list_dict['HS']['ND'],
                                  '{}/RTPT_IVH_ND_logo.png'.format(outdir))
    RTvsPT_ivh_sf_norm = plotlogo(RTvsPT_summarize_list_dict['HS']['SF'],
                                  '{}/RTPT_IVH_SF_logo.png'.format(outdir))
    RTvsPT_ivh_mix_norm = plotlogo(RTvsPT_summarize_list_dict['HS']['Mixed'],
                                   '{}/RTPT_IVH_Mixed_logo.png'.format(outdir))

    RTvsPT_imc_nd_norm = plotlogo(RTvsPT_summarize_list_dict['LS']['ND'],
                                  '{}/RTPT_IMC_ND_logo.png'.format(outdir))
    RTvsPT_imc_sf_norm = plotlogo(RTvsPT_summarize_list_dict['LS']['SF'],
                                  '{}/RTPT_IMC_SF_logo.png'.format(outdir))
    RTvsPT_imc_mix_norm = plotlogo(RTvsPT_summarize_list_dict['LS']['Mixed'],
                                   '{}/RTPT_IMC_Mixed_logo.png'.format(outdir))

    RTvsPT_nd_norm = plotlogo(
        RTvsPT_summarize_list_dict['HS']['ND'] + RTvsPT_summarize_list_dict['LS']['ND'],
        '{}/RTPT_ND_logo.png'.format(outdir))
    RTvsPT_sf_norm = plotlogo(
        RTvsPT_summarize_list_dict['HS']['SF'] + RTvsPT_summarize_list_dict['LS']['SF'],
        '{}/RTPT_SF_logo.png'.format(outdir))
    RTvsPT_mix_norm = plotlogo(
        RTvsPT_summarize_list_dict['HS']['Mixed'] + RTvsPT_summarize_list_dict['LS']['Mixed'],
        '{}/RTPT_Mixed_logo.png'.format(outdir))

    # finalize the pval
    TTvsNT_mask = np.full(len(features_TTvsNT_df_marked), False)
    TTvsNT_mask[features_TTvsNT_df_marked[features_TTvsNT_df_marked['GroupChange']=="SV"].index.values] = True
    features_TTvsNT_df_marked['pval_SV'] = np.where(TTvsNT_mask, np.minimum(TTvsNT_pvals[0], TTvsNT_pvals[1]),
                                                     TTvsNT_pvals[1])
    RTvsPT_mask = np.full(len(features_RTvsPT_df_marked), False)
    RTvsPT_mask[features_RTvsPT_df_marked[features_RTvsPT_df_marked['GroupChange'] == "SV"].index.values] = True
    features_RTvsPT_df_marked['pval_SV'] = np.where(RTvsPT_mask, np.minimum(RTvsPT_pvals[0], RTvsPT_pvals[1]),
                                                     RTvsPT_pvals[1])
    features_TTvsNT_df_dchange = features_TTvsNT_df_marked.loc[features_TTvsNT_df_marked['GroupChange']=='SV',:]
    features_RTvsPT_df_dchange = features_RTvsPT_df_marked.loc[features_RTvsPT_df_marked['GroupChange']=='SV',:]
    with pd.ExcelWriter("{}/tads_sub_array_marks_gcut{}_icut{}.xlsx".format(outdir,group_cutoff,individual_cutoff),
                        engine="xlsxwriter") as writer:
        features_TTvsNT_df_marked.to_excel(writer, index=False, sheet_name="TTvsNT_total_info")
        features_TTvsNT_df_dchange.to_excel(writer, index=False, sheet_name="TTvsNT_SV")
        features_RTvsPT_df_marked.to_excel(writer, index=False, sheet_name="RTvsPT_total_info")
        features_RTvsPT_df_dchange.to_excel(writer, index=False, sheet_name="RTvsPT_SV")

    # summary
    TTvsNT_bio_mark = Counter(features_TTvsNT_df_marked['GroupChange'])
    RTvsPT_bio_mark = Counter(features_RTvsPT_df_marked['GroupChange'])
    summary_dict = {'GroupChange':[],'Type':[],'Count':[]}
    for key in TTvsNT_bio_mark:
        summary_dict['GroupChange'].append(key)
        summary_dict['Type'].append('TT.vs.NT')
        summary_dict['Count'].append(TTvsNT_bio_mark[key])

    for key in RTvsPT_bio_mark:
        summary_dict['GroupChange'].append(key)
        summary_dict['Type'].append('RT.vs.PT')
        summary_dict['Count'].append(RTvsPT_bio_mark[key])
    summary_df = pd.DataFrame(summary_dict)
    plt.figure()
    sns.barplot(data=summary_df,x='GroupChange',y='Count',hue='Type',order=['C','MV','SV'])
    plt.savefig('{}/TADs_array_GroupChange_summary.png'.format(outdir))
    plt.close()

    # individual specific number
    # indi_spe_list = ['high_high_high','high_medium_high','high_low_high']
    # plot clustermap
    value2int = {'high':3.0,'medium':2.0,'low':1.0}
    colors_list = ['limegreen','goldenrod','red']
    # high: palevioletred, low:tan, con:forestgreen
    features_TTvsNT_df_marked['GroupChange'] = pd.Categorical(features_TTvsNT_df_marked['GroupChange'],
                                               categories=['C','MV','SV'],ordered=True)
    features_RTvsPT_df_marked['GroupChange'] = pd.Categorical(features_RTvsPT_df_marked['GroupChange'],
                                                           categories=['C','MV','SV'],ordered=True)
    features_TTvsNT_list = []
    features_TTvsNT_len_list = []
    for group,df in features_TTvsNT_df_marked.groupby('GroupChange'):
        print(group)
        if group in ['C','MV']:
            current_df = df[['TT.vs.NT','Intra-NT','Intra-TT']].replace(value2int)
            current_df.sort_values(['TT.vs.NT','Intra-NT','Intra-TT'],inplace=True)
            features_TTvsNT_list.append(current_df)
            features_TTvsNT_list.append(pd.DataFrame(0,index=np.arange(30),columns=list(current_df.columns)))
            features_TTvsNT_len_list.append(len(current_df))
        else:
            current_indi = df[(df['TT.vs.NT']=='high')&(df['Intra-TT']=='high')][['TT.vs.NT','Intra-NT','Intra-TT']].replace(value2int)
            current_indi.sort_values(['TT.vs.NT','Intra-NT','Intra-TT'],inplace=True)
            features_TTvsNT_list.append(current_indi)
            features_TTvsNT_list.append(pd.DataFrame(0,index=np.arange(30),columns=list(current_indi.columns)))
            features_TTvsNT_len_list.append(len(current_indi))
            current_noindi = df[(df['TT.vs.NT']=='high')&(df['Intra-TT']!='high')][['TT.vs.NT','Intra-NT','Intra-TT']].replace(value2int)
            current_noindi.sort_values(['TT.vs.NT','Intra-NT','Intra-TT'],inplace=True)
            features_TTvsNT_list.append(current_noindi)
            features_TTvsNT_len_list.append(len(current_noindi))

    features_TTvsNT_final_df = pd.concat(features_TTvsNT_list).T
    cmap = custom_discrete_cmap(4,['white','limegreen','goldenrod','red'],1)
    plt.figure(figsize=(15,6))
    ax = sns.heatmap(features_TTvsNT_final_df, xticklabels=False, cmap=cmap)
    ax.set_title('TT.vs.NT')
    ax.hlines([1, 2], colors='black', *ax.get_xlim())
    # ax.set(xticklabels=[])
    # ax.set(xlabel=None)
    # ax.tick_params(bottom=False)
    plt.savefig(f'{outdir}/TTvsNT_heatmap.png')
    plt.close()

    features_RTvsPT_list = []
    features_RTvsPT_len_list = []
    for group,df in features_RTvsPT_df_marked.groupby('GroupChange'):
        print(group)
        if group in ['C','MV']:
            current_df = df[['RT.vs.PT','Intra-PT','Intra-RT']].replace(value2int)
            current_df.sort_values(['RT.vs.PT','Intra-PT','Intra-RT'],inplace=True)
            features_RTvsPT_list.append(current_df)
            features_RTvsPT_list.append(pd.DataFrame(0,index=np.arange(30),columns=list(current_df.columns)))
            features_RTvsPT_len_list.append(len(current_df))
        else:
            current_indi = df[(df['RT.vs.PT']=='high')&(df['Intra-RT']=='high')][['RT.vs.PT','Intra-PT','Intra-RT']].replace(value2int)
            current_indi.sort_values(['RT.vs.PT','Intra-PT','Intra-RT'],inplace=True)
            features_RTvsPT_list.append(current_indi)
            features_RTvsPT_list.append(pd.DataFrame(0,index=np.arange(30),columns=list(current_indi.columns)))
            features_RTvsPT_len_list.append(len(current_indi))
            current_noindi = df[(df['RT.vs.PT']=='high')&(df['Intra-RT']!='high')][['RT.vs.PT','Intra-PT','Intra-RT']].replace(value2int)
            current_noindi.sort_values(['RT.vs.PT','Intra-PT','Intra-RT'],inplace=True)
            features_RTvsPT_list.append(current_noindi)
            features_RTvsPT_len_list.append(len(current_noindi))

    features_RTvsPT_final_df = pd.concat(features_RTvsPT_list).T
    cmap = custom_discrete_cmap(4,['white','limegreen','goldenrod','red'],1)
    plt.figure(figsize=(15,6))
    ax = sns.heatmap(features_RTvsPT_final_df, xticklabels=False, cmap=cmap)
    ax.set_title('RT.vs.PT')
    ax.hlines([1, 2], colors='black', *ax.get_xlim())
    # ax.set(xticklabels=[])
    # ax.set(xlabel=None)
    # ax.tick_params(bottom=False)
    plt.savefig(f'{outdir}/RTvsPT_heatmap.png')
    plt.close()

    # calculate correlations for all sample based on tads_array
    total_size = 0
    tadarray_size_dict = {}
    for tadarray in tads_sub_array_dict:
        current_tadarray = tads_sub_array_dict[tadarray]
        max_size = max((sum(np.abs(v))) for k, v in current_tadarray.items())
        total_size += max_size
        tadarray_size_dict[tadarray] = max_size
    sample_tad_corr_dict = {}
    for tadarray in tads_sub_array_dict:
        current_tadarray = tads_sub_array_dict[tadarray]
        # pearson correlation should have same length and at least 2 length
        max_length = max((len(v)) for k, v in current_tadarray.items())
        max_size = tadarray_size_dict[tadarray]
        tadarray_weight = max_size / total_size
        if max_length == 1:
            max_length = 2
        new_tadarray = {}
        for sample in current_tadarray:
            sample_array = current_tadarray[sample].tolist() + (max_length - len(current_tadarray[sample])) * [0]
            # avoid NAN generated based no variance within sample, e.g., [1,1,1]
            if len(np.unique(sample_array)) == 1:
                sample_array[0] -= 0.1
                sample_array[1] += 0.1
            new_tadarray[sample] = sample_array
        samples_tad_corr_mat = pd.DataFrame(0, columns=samples, index=samples)
        for subset in itertools.combinations(samples, 2):
            sample1 = subset[0]
            sample2 = subset[1]
            pearsonr_result = pearsonr(new_tadarray[sample1],
                                       new_tadarray[sample2])
            samples_tad_corr_mat.loc[sample1, sample2] = pearsonr_result[0]
            samples_tad_corr_mat.loc[sample2, sample1] = pearsonr_result[0]
        # > 0.7 to 1; >0.3 <0.7 to 0.5; <0.3 to 0
        samples_tad_corr_mat[samples_tad_corr_mat >= 0.3] = 1
        samples_tad_corr_mat[(samples_tad_corr_mat > -0.3) & (samples_tad_corr_mat < 0.3)] = 0.5
        samples_tad_corr_mat[samples_tad_corr_mat <= -0.3] = 0
        sample_tad_corr_dict[tadarray] = tadarray_weight * samples_tad_corr_mat

    samples_tad_corr_mat_final = sum(sample_tad_corr_dict.values())
    samples_tad_corr_mat_final.values[np.diag_indices_from(samples_tad_corr_mat_final)] = 1
    get_tri_heatmap(samples_tad_corr_mat_final, loc='lower', show_text=False, vmin=0.6, vmax=0.9,
                    output='{}/samples_pearsonr_heatmap.png'.format(outdir))
    tissueType = ['NT', 'PT', 'RT']
    sample_tad_merge_corr_mat = pd.DataFrame(0, columns=tissueType, index=tissueType)
    sample_tad_merge_corr_mat.loc['NT', 'NT'] = samples_tad_corr_mat_final.loc['NT1', 'NT2']
    sample_tad_merge_corr_mat.loc['NT', 'PT'] = samples_tad_corr_mat_final.loc[
        ['NT1', 'NT2'], ['PT{}'.format(i) for i in range(1, 6)]].mean().mean()
    sample_tad_merge_corr_mat.loc['NT', 'RT'] = samples_tad_corr_mat_final.loc[
        ['NT1', 'NT2'], ['RT{}'.format(i) for i in range(1, 6)]].mean().mean()
    sample_tad_merge_corr_mat.loc['PT', 'PT'] = samples_tad_corr_mat_final.loc[['PT{}'.format(i) for i in range(1, 6)],
    ['PT{}'.format(i) for i in range(1, 6)]].mean().mean()
    sample_tad_merge_corr_mat.loc['PT', 'RT'] = samples_tad_corr_mat_final.loc[['PT{}'.format(i) for i in range(1, 6)],
    ['RT{}'.format(i) for i in range(1, 6)]].mean().mean()
    sample_tad_merge_corr_mat.loc['RT', 'RT'] = samples_tad_corr_mat_final.loc[['RT{}'.format(i) for i in range(1, 6)],
    ['RT{}'.format(i) for i in range(1, 6)]].mean().mean()
    get_tri_heatmap(sample_tad_merge_corr_mat, loc='upper', fontsize=30, show_text=True, vmin=0.65, vmax=0.85,
                    output='{}/samples_mean_pearsonr_heatmap.png'.format(outdir))
