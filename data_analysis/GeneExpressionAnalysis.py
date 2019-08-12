#!/usr/bin/env python3

import datetime
import itertools
import math
import os
import re
import subprocess
import urllib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from sklearn import decomposition


#############
## CLASSES ##
#############


class Preprocessor():
    
    def __init__(self):
        
        return
    
    
    def _extract_GEMS_tags(self, value):
        
        value = value.split('/')[-1]
        pattern = '(?:(R_UZG.*?)|(?:(?:R_)?(.*?)(_?[CM]?)?))_(S[0-9]+)_(L00[0-9])_R1_001[.]sam'
        tag0, tag1, sample_type, sample, lane = re.findall(string=value, pattern=pattern)[0]
        tag = ''.join([tag0, tag1])
    
        tag = tag.replace('_', '').upper()
        sample_type = sample_type.replace('_', '')
        
        if 'DIN' in tag:
            tag = tag.replace('DIN', 'DINANT')
        if 'ELI' in tag:
            tag = tag.replace('ELI', 'STELI')
        if 'RUZG-2' in tag:
            tag = tag.replace('-', '')
        if 'BEN-006' in tag:
            tag = 'BEN-006C'
            sample_type = ''
    
        return tag, sample_type, sample, lane
    
    
    def clean(self, df, rows_to_drop=['_unmapped', '_ambiguous', '_no_feature'], data_origin=None):
        
        # remove rows with reads we don't care about (ie not assigned to a specific gene)
        df.drop(rows_to_drop, inplace=True)
        
        # fill in all missing read counts with 0 (missing read values originate from genes observed only in a single file)
        df.fillna(0, inplace=True)
    
        # sort both columns and rows by mean read count (from high to low)
        # transpose dataframe to have genes as columns and samples as rows
        df = df.reindex(df.mean(axis=1).sort_values(ascending=False).index)
        df = df.T.reindex(df.mean(axis=0).sort_values(ascending=False).index)
    
        # drop duplicate rows (= samples), based on 5000 first gene expression values (using higher numbers can cause recursion errors)
        df = df.drop_duplicates(subset=df.columns[:5000])
        
        # extract meta information from sample paths
        run_pattern = '/([0-9]{6}_NB[0-9]{6}_[0-9]{4}_.{10})/'
        df.loc[:, 'run'] = pd.Series(df.index).apply(lambda x: re.findall(run_pattern, x)[0]).values
        
        lane_pattern = '_(L00[1-4])_R'
        df.loc[:, 'lane'] = pd.Series(df.index).apply(lambda x: re.findall(lane_pattern, x)[0]).values
        
        name_pattern = '.*/(.*?)_L00[1-4]_R'
        df.loc[:, 'name'] = pd.Series(df.index).apply(lambda x: re.findall(name_pattern, x)[0]).values
        
        if data_origin == 'GEMS':
            
            df.loc[:, 'tags'] = pd.Series(df.index).apply(lambda x: self._extract_GEMS_tags(x)[0]).values
            df.loc[:, 'sample_type'] = pd.Series(df.index).apply(lambda x: self._extract_GEMS_tags(x)[1]).values
            df.loc[:, 'sample'] = pd.Series(df.index).apply(lambda x: self._extract_GEMS_tags(x)[2]).values
            
            df.loc[:, 'name'] = df[['tags', 'sample_type', 'sample']].apply(lambda x: '_'.join(x), axis=1)
        
        elif data_origin == 'MMR':
            
            repeat_pattern = 'EXP[037]_([12])_'
            df.loc[:, 'replicates'] = pd.Series(df.index).apply(lambda x: re.findall(repeat_pattern, x)[0]).values
            
            sample_pattern = '(M[0-9]+)'
            df.loc[:, 'sample'] = pd.Series(df.index).apply(lambda x: re.findall(sample_pattern, x)[0]).values
            
            day_pattern = 'EXP([037])'
            df.loc[:, 'day'] = pd.Series(df.index).apply(lambda x: re.findall(day_pattern, x)[0]).values
            
        return df
    
    
    def combine_lane_counts(self, df, sample_id='name'):
        
        print('Collecting metadata...')
        meta_cols = [col for col in df.columns if 'gene' not in col]
        print(', '.join(meta_cols))
        print('Grouping rows by \'{}\' and summing counts'.format(sample_id))
        meta_df = df[meta_cols].drop(['lane'], axis=1).drop_duplicates().set_index(sample_id, drop=False)
        summed_lanes = df.groupby(sample_id).sum()    
        combined_df = summed_lanes.join(meta_df)
        print('Grouping done!\nFinal shape of df: {}\n'.format(combined_df.shape))
        
        return combined_df
    


class LabelProcessor():
    
    def __init__(self, df, data_origin='GEMS'):
        
        if data_origin == 'GEMS':
            self._make_label_df()    
            self.full_df = self._assign_groups(df)
        else:
            self.full_df = df
            
        self.internal_df = None
        
        
    def _make_label_df(self):
        
        self.meningitis_labels_file = '../data/Code meningitis stalen.xlsx'
        self.reuma_labels_file = '../data/reuma_studies_coded.xlsx'
        
        self.mening_df = pd.read_excel('../data/Code meningitis stalen.xlsx')
        self.mening_df = self.mening_df.iloc[:, :2].fillna('drop')
        self.mening_df = self.mening_df[self.mening_df['Data code'] != 'drop']        
        self.mening_df.columns = ['tags', 'group']
        
        self.reu_df = pd.read_excel('../data/reuma_studies_coded.xlsx')
        self.reu_df = self.reu_df[['groep', 'Nummer tube']]
        self.reu_df['Nummer tube'] = self.reu_df['Nummer tube'].apply(lambda x: x.split(' ')[0])
        self.reu_df['Nummer tube'] = self.reu_df['Nummer tube'].apply(lambda x: 'R'+x if x.startswith('UZG') else x)
        self.reu_df.columns = ['group', 'tags']
        
        self.label_df = pd.concat([self.reu_df, self.mening_df])
        self.label_df['tags'] = self.label_df['tags'].apply(lambda x: x.replace(' ', '').upper())
        
        
    def _assign_groups(self, df, based_on='tags', group_name='group'):
            
        df.loc[:, group_name] = df.loc[:, based_on].apply(self._find_group).values
        return df
    
    
    def _find_group(self, value):
        
        ls = self.label_df[self.label_df['tags']==value]['group'].values
        if len(ls) == 1:
            return ls[0]
        else:
            return 'Other'
    
    
    def add_samples(self, column, value):
        
        if self.internal_df is None:
            self.internal_df = self.full_df[self.full_df[column] == value]
            
        else:
            self.tmp_df = self.full_df[self.full_df[column] == value]
            self.internal_df = pd.concat([self.internal_df, self.tmp_df])
    
        self.internal_df = self.internal_df.fillna('')
        
    
    def remove_samples(self, column, value):
        
        self.internal_df = self.internal_df[self.internal_df[column] != value]
    
    
    def regroup_meta_tags(self, tag, new_tag, group_dict):
        
        self.internal_df.loc[:, new_tag] = ''
        for new_group, old_values in group_dict.items():
            for value in old_values:
                self.internal_df.loc[self.internal_df[tag] == value , new_tag] = new_group
        
        
    def combine_columns(self, upper_col, lower_col, new_col_name):
        
        self.internal_df.loc[:, new_col_name] = ''
        for upper_value in self.internal_df[upper_col].unique():
            c = 0
            for lower_value in self.internal_df[self.internal_df[upper_col]==upper_value][lower_col].unique():
                mask = (self.internal_df[upper_col] == upper_value) & (self.internal_df[lower_col] == lower_value)
                self.internal_df.loc[mask, new_col_name] = c
                c += 1
                
        
    def generate_DESeq_data(self, read_path='../data/readcounts_test.txt', col_data_path='../data/col_data_test.txt'):
        
        gene_cols = [col for col in self.internal_df.columns if 'gene' in col]
        meta_cols = [col for col in self.internal_df.columns if 'gene' not in col]
        
        if len(gene_cols) == 0:
            raise ValueError('''No gene columns found.\n
                             Gene columns are identified by looking for columns with the string "gene" in them.''')
                
        # gene read count
        read_df = self.internal_df[gene_cols].T
        read_df.columns = self.internal_df['name'].values
        read_df.rename_axis('genename', inplace=True)
        read_df.to_csv(read_path, sep='\t', index=True, header=True)
        
        # meta data
        col_data = self.internal_df[meta_cols]
        col_data.to_csv(col_data_path, sep='\t', index=False)
    
        return read_path, col_data_path
    
    
    
class DESeq2run():
    
    def __init__(self, read_path, col_path, sample_column, test='Wald'):
        
        self.read_path = read_path
        self.col_path = col_path
        self.sample_column = sample_column
        self.test = test
        
        self.data_dir = '/'.join(os.path.abspath(self.read_path).split('/')[:-1])
        self.meta_data = pd.read_csv(self.col_path, sep='\t')
        self.metatags = self.meta_data.columns.values
    
        self.normalized_counts = None
    
    
    def calc_DEGs(self, design, to_test, R_script_dir):
        
        results = []
        for contrast, subgroups in to_test.items():
            for contrast_subgroup1, contrast_subgroup2 in itertools.combinations(subgroups, 2):
                script_path, deg_path = self._generate_R_script(design, contrast, contrast_subgroup1, contrast_subgroup2, R_script_dir)
                results.append(deg_path)
                #continue
                print('Starting R process - {}: {} vs {}'.format(contrast, contrast_subgroup1, contrast_subgroup2))
                proc = subprocess.Popen('Rscript {}'.format(script_path).split(), stdout=subprocess.PIPE)
                (out, err) = proc.communicate()
        
        self.normalized_counts = pd.read_csv('{}/DESeq2_normalized_readcounts.txt'.format(self.data_dir))
        self.normalized_counts.set_index('Unnamed: 0', inplace=True)
        
        return results
                
                
    def _generate_R_script(self, design, contrast, contrast_subgroup1, contrast_subgroup2, R_script_dir):
        
        script_path = '{}/generated_deseq2_analysis_{}-{}vs{}.R'.format(R_script_dir,
                                                                       contrast,
                                                                       contrast_subgroup1,
                                                                       contrast_subgroup2)
        with open(script_path, 'w') as o:
            
            o.write('''# This script was automatically generated from {}.
                    # Date: {}
                    
                    library(DESeq2)
                    library("pheatmap")
                    
                    setwd("{}")
                    
                    datafilepath <- "{}"
                    data <- read.table(datafilepath, header = TRUE, sep = "\t", row.names=1)
                    
                    coldatafilepath <- "{}"
                    colData <- read.table(coldatafilepath, header = TRUE, sep = "\t", check.names = FALSE)

                    '''.format(os.path.realpath(__file__),
                               datetime.datetime.now(),
                               self.data_dir,
                               os.path.abspath(self.read_path),
                               os.path.abspath(self.col_path)
                               ).replace('    ', ''),
                    )
    
            for metatag in self.metatags:
                o.write('''colData${} <- as.factor(colData${})
                        '''.format(metatag, metatag).replace('    ', ''))
        
            deg_path = '{}/DESeq2_results_{}-{}vs{}.txt'.format(self.data_dir,
                                                               contrast,
                                                               contrast_subgroup1,
                                                               contrast_subgroup2)
            o.write('''
                    dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ {})
                    
                    num_pats <- length(unique(colData${}))
                    dds <- dds[ rowSums(counts(dds)) > 5*num_pats, ]
                    
                    dds <- DESeq(dds, test="{}")
                    
                    counts <- counts(dds, normalized=TRUE)
                    write.csv(as.data.frame(counts), file = "{}/DESeq2_normalized_readcounts.txt")
                    
                    res <- results(dds, contrast = c("{}", "{}", "{}"))
                    
                    res <- res[order(res$padj),]
                    
                    write.csv(as.data.frame(res), file = "{}")
                    '''.format(design, self.sample_column, self.test, self.data_dir,
                               contrast, contrast_subgroup1, contrast_subgroup2, deg_path).replace('    ', ''))
            
            return script_path, deg_path
    
   
    def PCA(self, color_by, show=False):
        
        pca = decomposition.PCA(n_components=4, random_state=0)
        
        X = pd.DataFrame(pca.fit_transform(self.normalized_counts))
        if color_by != None:
            Y = self.meta_data[color_by]
            X = pd.concat([X, Y], ignore_index=True, axis=1)
            X.columns = [str(x) for x in X.columns[:-1]]+['hue']
        else:
            X.columns = [str(x) for x in X.columns]
        
        for i in range(3):
            
            plt.subplots(figsize=[8, 8])
            for factor in X['hue'].unique():
                factor_mask = X['hue'] == factor
                plt.scatter(X.loc[factor_mask, str(i)], X.loc[factor_mask, str(i+1)], label=factor)
            plt.xlabel('PC{} - {}% variance'.format(i+1, pca.explained_variance_ratio_[i]*100))
            plt.ylabel('PC{} - {}% variance'.format(i+2, pca.explained_variance_ratio_[i+1]*100))
            plt.legend()
            plt.savefig('{}/PCA_normalized_all_genes_{}_PC{}_PC{}.png'.format(self.data_dir, color_by, i+1, i+2), dpi=600)
            
            if show == True:
                plt.show()        
            
            plt.clf()
        


class DESeq2results():
    
    def __init__(self, results_file, gene_mapping_file='../original_data/genes.txt'):
        
        self.data_dir = '/'.join(results_file.split('/')[:-1])
        
        # DESeq2 results
        self.full_deseq_result = pd.read_csv(results_file)
        self.full_deseq_result.set_index('Unnamed: 0', inplace=True)
        
        # gene mapper
        gene_df = pd.read_csv(gene_mapping_file, sep='\t', header=None, prefix='column_')
        self.gene_map = gene_df[['column_0', 'column_3']]

        
    def get_sig_DEGs(self, write_to, padj_threshold=0.05, fold_threshold=1):
        
        # select DEGs
        self.deseq_result = self.full_deseq_result.loc[self.full_deseq_result['padj'] <= padj_threshold, :]
        greater_fold_mask = self.deseq_result['log2FoldChange'] >= fold_threshold
        smaller_fold_mask = self.deseq_result['log2FoldChange'] <= -fold_threshold
        threshold_mask = np.array(greater_fold_mask) | np.array(smaller_fold_mask)
        self.deseq_result = self.deseq_result.loc[threshold_mask, :]
        
        # convert gene names
        alt_genes = pd.DataFrame(pd.Series(self.deseq_result.index, name='genes').apply(lambda x: self._lookup(self.gene_map, 'column_0', 'column_3', x)))
        alt_genes.index = self.deseq_result.index
        self.deseq_result = pd.concat([self.deseq_result, alt_genes], axis=1).sort_values(by='padj', ascending=True)
        
        self.deseq_result.to_csv(write_to, sep='\t')
        
        return self.deseq_result
          
        
    def _lookup(self, df, key, target, value):
        
        return df.loc[df[key] == value, target].values[0]
    
    
    def PCA(self, color_by=None, readcountfile=None, coldatafile=None, show=False, filename=None):
        
        pca = decomposition.PCA(n_components=4, random_state=0)

        if readcountfile == None:
            readcountfile = '{}/DESeq2_normalized_readcounts.txt'.format(self.data_dir)
            
        if coldatafile == None:
            coldatafile = '{}/col_data.txt'.format(self.data_dir)
            
        meta_data = pd.read_csv(coldatafile, sep='\t')
        
        norm_read_counts = pd.read_csv(readcountfile)
        norm_read_counts.set_index('Unnamed: 0', inplace=True)
        
        sig_genes = self.deseq_result
        sig_genes = sig_genes.index
        
        X = pd.DataFrame(pca.fit_transform(norm_read_counts.loc[sig_genes, :]))
        if color_by != None:
            Y = meta_data[color_by]
        else:
            Y = pd.Series(['' for x in range(X.shape[0])])
        X = pd.concat([X, Y], ignore_index=True, axis=1)
        X.columns = [str(x) for x in X.columns[:-1]]+['hue']
        
        for i in range(3):
            
            plt.subplots(figsize=[16, 16])
            for factor in X['hue'].unique():
                factor_mask = X['hue'] == factor
                plt.scatter(X.loc[factor_mask, str(i)], X.loc[factor_mask, str(i+1)], label=factor)
            plt.xlabel('PC{} - {}% variance'.format(i+1, pca.explained_variance_ratio_[i]*100))
            plt.ylabel('PC{} - {}% variance'.format(i+2, pca.explained_variance_ratio_[i+1]*100))
            plt.legend()
            if filename == None:
                plt.savefig('{}/PCA_normalized_sig_genes_{}_PC{}_PC{}.png'.format(self.data_dir, color_by, i+1, i+2), dpi=600)
            else:
                plt.savefig('{}/{}_PC{}_PC{}'.format(self.data_dir, filename, i+1, i+2), dpi=600)
            
            if show == True:
                plt.show()
            
            plt.clf()
        
        
    def heatmap(self, filename, readcountfile=None, coldatafile=None, label='group'):
    
        if readcountfile == None:
            readcountfile = '{}/DESeq2_normalized_readcounts.txt'.format(self.data_dir)
            
        if coldatafile == None:
            coldatafile = '{}/col_data.txt'.format(self.data_dir)
        
        meta_data = pd.read_csv(coldatafile, sep='\t')
        
        norm_read_counts = pd.read_csv(readcountfile)
        norm_read_counts.set_index('Unnamed: 0', inplace=True)
        
        sig_genes = self.deseq_result
        sig_genes = sig_genes.index
        
        X = norm_read_counts.loc[sig_genes, :].T
        Y = meta_data[label]
        
        X = X.assign(index=pd.Series(Y).values).set_index('index')
        X = X.sort_index()
        
        plt.subplots(figsize=(20, 20))
        plt.tight_layout()
        sns.heatmap(X.applymap(log))
        plt.savefig('{}/{}'.format(self.data_dir, filename), dpi=600)
        plt.clf()
        
        plt.subplots(figsize=(20, 20))
        
        colors = sns.color_palette()[: len(Y.unique())]
        color_dict = dict(zip(Y.unique(), colors))
        mapped_colors = [color_dict[x] for x in Y]
        sns.clustermap(X.applymap(log), row_colors=mapped_colors)
        plt.savefig('{}/clustered_{}'.format(self.data_dir, filename), dpi=600)
        plt.clf()
        
        
        
###############
## FUNCTIONS ##
###############


def log(value):
    
    try:
        return math.log(value)
    except ValueError:
        return 0
    
    
def collect_data(directory, file_pattern, index='genename'):
    
    '''
    Search through all files in directory and its subdirectories.
    Read in all files that match the regular expression file_pattern as pandas dataframes and concatenate them together.
    
    Arguments:
    ----------
    - directory: string - path to top directory to start search from
    - file_pattern: string - regular expression pattern to match file paths with
    '''
    
    # find all files in directory and subdirectories that match file_pattern
    raw_readcounts = []
    
    for dr, subdr, files in os.walk(directory):
        print('Searching files in {}...'.format(dr))
        for file in files:
            path = '{}/{}'.format(dr, file)
            if len(re.findall(file_pattern, path)) > 0:
                print('File found: {}'.format(file))
                raw_readcounts.append(path)
    print('Files found:\n{}'.format('\n'.join(raw_readcounts)))

    # read in each file as a pandas dataframe and concatenate them together
    to_concat = []
    
    print('\nReading in files:')
    for file in raw_readcounts:
        df = pd.read_csv(file, sep='\t')
        print('{} with shape: {}'.format(file, df.shape))
        df.set_index(index, inplace=True)
        df.fillna(0, inplace=True)
        unnamed_mask = [col for col in df.columns if 'Unnamed' not in col]
        to_concat.append(df.loc[:, unnamed_mask])
    
    print('Concatenating dataframes...')
    df = pd.concat(to_concat, axis=1, sort=True)
    print('Concatenation done!\nFinal shape of df: {}\n'.format(df.shape))
    df = df.reindex(sorted(df.columns), axis=1)
    return df
    

    
if __name__ == '__main__':
    
    path = '/Users/nicolasdeneuter/Dropbox/PhD/Projects/GOA/MMR/RNAseq/data/unsupervised/DESeq2_results_day-0vs3.txt'
    deg_result = DESeq2results(results_file=path)
    deg_result.get_sig_DEGs(write_to=path.replace('.txt', '_fold1_sig.txt'))