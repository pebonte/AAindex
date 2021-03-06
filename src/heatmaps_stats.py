#! /usr/bin/env python3

import os
import sys
import glob
import re
import math
import seaborn as sns; sns.set(style="whitegrid", font_scale = 1.2)
import scikit_posthocs as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from aaindex_processing import (retrieve_values_aaindex,
                                retrieve_loop_data,
                                get_min_and_max,
                                make_dataframe_from_aaindex_data,
                                make_aaindex_boxplots_by_family)
'''
Script that does statistical tests from alignments and aaindex values.

REQUIREMENT : Internet

'''

AMINO_ACID_LIST = []
aa_list = 'ARNDCQEGHILKMFPSTWYV'
for aa in aa_list:
    AMINO_ACID_LIST.append(aa)

if __name__ == '__main__':
    # Retrieving alignments data
    data = retrieve_loop_data('../data/alignment_plain/', '*.plain')
    alignment_data_dict, remaining = data
    # Files that won't be run by the script
    files_remaining_string = 'List of files with more than 6 cysteines : '
    for file_remaining in remaining:
        files_remaining_string += file_remaining + ', '
    print(files_remaining_string)
    # AAindex that will be retrieved
    list_aaindex_ids = ['BHAR880101', 'CASG920101', 'CHAM830107',
                        'CHOP780201', 'CHOP780202', 'CHOP780203',
                        'CIDH920105', 'DAYM780201', 'EISD860102',
                        'FASG760101', 'FAUJ880111', 'FAUJ880111',
                        'FAUJ880112', 'GOLD730101', 'GRAR740102',
                        'JANJ780101', 'JANJ780102', 'JANJ780103',
                        'JOND750101', 'JOND920102', 'KLEP840101',
                        'KRIW790101', 'KRIW790102', 'KYTJ820101',
                        'MITS020101', 'PONP930101', 'RACS820114',
                        'RADA880108', 'TAKK010101', 'TAKK010101',
                        'VINM940101', 'WARP780101', 'WOLR790101',
                        'ZIMJ680101']
    #list_aaindex_ids = ['GRAR740102']
    # Retrieving AAindex data
    aaindex_data = retrieve_values_aaindex(list_aaindex_ids)
    web_aaindex_data_dict, aaindex_names_data = aaindex_data
    # Group both data into a dataframe
    df_data = make_dataframe_from_aaindex_data(web_aaindex_data_dict,
                                               aaindex_names_data,
                                               alignment_data_dict,
                                               list_aaindex_ids)
    data, dataframe, coding_family_name_dict, list_aaindex_ids = df_data
    # Figure Parameters
    HEIGHT = 2100 * 1.5
    WIDTH = 1270 * 1.5
    # Format: diagonal, non-significant, p<0.001, p<0.01, p<0.05
    cmap = ['1', '#D6D6D6',  '#08306b',  '#4292c6', '#c6dbef']
    for index in list_aaindex_ids:
        aaindex_df = dataframe.loc[dataframe['AAindex'] == index]
        loop_list = list(set(data['Loop']))
        family_list = list(set(data['Family']))
        dirName = '../results/AAindex_stats/{}'.format(index)
        if not os.path.exists(dirName):
            os.mkdir(dirName)
        for loop in loop_list:
            try:
                print('Loop {} for AAindex {} : {}'.format(loop, index, aaindex_names_data[index]))
            except:
                print('Loop {} for AAindex {}'.format(loop, index))
            try:
                fig = plt.figure(figsize=(12, 7))
                loop_aaindex_df = aaindex_df.loc[aaindex_df['Loop'] == loop]
                # Mann-Whitney test
                pc = sp.posthoc_mannwhitney(loop_aaindex_df,
                                            val_col='Value',
                                            group_col='Family',
                                            p_adjust = 'fdr_bh')
                # Demi-matrix mask
                mask = np.zeros_like(pc)
                mask[np.triu_indices_from(mask)] = True
                ax = plt.axes(title = 'P-values loop {} for {}'.format(loop, index))
                ax.set_title('P-values loop {} for {}'.format(loop, index), fontsize=35)
                heatmap_args = {'cmap': cmap,
                                'mask': mask,
                                'ax': ax,
                                'xticklabels': True,
                                'yticklabels': True,
                                'clip_on': False,
                                'square': True,
                                'cbar_ax_bbox': [0.80, 0.35, 0.04, 0.3]}
                sp.sign_plot(pc, **heatmap_args)
                DPI = (fig.get_dpi()) * 2
                fig.set_size_inches(float(HEIGHT)/float(DPI),float(WIDTH)/float(DPI))
                # Saving figure
                fig.savefig('{}/AAindex_stats_aaindex_{}_loop_{}'.format(dirName, index, loop))
                plt.close(fig)
                print('OK\n')
                #plt.show()
            except:
                print('didnt work\n')
                pass
    print('\nFINISH')
