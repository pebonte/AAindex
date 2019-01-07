#! /usr/bin/env python3

import os
import sys
import glob
import re
import time
import math
import seaborn as sns; sns.set(style="whitegrid")
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

AMINO_ACID_LIST = []
aa_list = 'ARNDCQEGHILKMFPSTWYV'
for aa in aa_list:
    AMINO_ACID_LIST.append(aa)

data = retrieve_loop_data('../data/alignment_plain/', '*.plain')
alignment_data_dict, remaining = data

files_remaining_string = 'List of files with more than 6 cysteines : '
for file_remaining in remaining:
    files_remaining_string += file_remaining + ', '
print(files_remaining_string)

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
list_aaindex_ids = ['CHOP780207']
aaindex_data = retrieve_values_aaindex(list_aaindex_ids)
web_aaindex_data_dict, aaindex_names_data = aaindex_data
df_data = make_dataframe_from_aaindex_data(web_aaindex_data_dict,
                                        aaindex_names_data,
                                        alignment_data_dict,
                                        list_aaindex_ids)
data, dataframe, coding_family_name_dict, list_aaindex_ids = df_data

for family in list(sorted(list(alignment_data_dict['number_of_seq'].keys()))):
    make_aaindex_boxplots_by_family(family,
                                    list_aaindex_ids,
                                    dataframe,
                                    coding_family_name_dict,
                                    alignment_data_dict,
                                    aaindex_names_data)
    
