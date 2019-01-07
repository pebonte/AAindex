#! /usr/bin/env python3

import os
import sys
import glob
import re
import time
import math
from io import StringIO
import numpy as np
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

AMINO_ACID_LIST = []
aa_list = 'ARNDCQEGHILKMFPSTWYV'
for aa in aa_list:
    AMINO_ACID_LIST.append(aa)

chrome_options = Options()
chrome_options.add_argument("--headless")

def retrieve_values_aaindex(aaindex_ids):
    '''
    Function to retrieve aaindex values from the online database genome.jp

    ARGUMENTS:
    - aaindex_ids : List of aaaindex you want to retrieve data.
    '''
    start = time.time()
    # Create web driver object that will simulate a Chrome navigator
    br = webdriver.Chrome('../chromedriver', chrome_options=chrome_options)
    # Max time in order to find element in webpage
    br.implicitly_wait(120)
    # Set dictionnaries to store data
    aaindex_dico = {}
    aaindex_name = {}
    for aaindex in aaindex_ids:
        aaindex_values = {}
        aaindex_dico[aaindex] = aaindex_values
        # Open webpage in navigator
        br.get('https://www.genome.jp/dbget-bin/www_bget?aaindex:{}'.\
               format(aaindex))
        # Retrieve web data from html xpath in a web object data
        web_data = br.find_element_by_xpath('//*[@id="wrapper"]/pre')
        # Retrieve text data from web object data
        values_data = web_data.text
        # Regex to find values
        regex_values = re.compile('^( +-?[0-9]+(\.)?([0-9]+)?){10}')
        # Regex to find aaindex name
        regex_aaindex_name = re.compile('^D ([^\(]+)')
        # Transform data into a StringIO object (needed to parse by line)
        values_data = StringIO(values_data)
        values = ''
        for line in values_data:
            if regex_values.search(line):
                # Concatenate the two lines with aaindex values
                values += regex_values.search(line).group(0)
            if regex_aaindex_name.search(line):
                name = regex_aaindex_name.search(line).group(1)
                aaindex_name[aaindex] = name
        # Format and transform both values lines into a list
        values = values.replace('    ', ' ').replace('   ', ' ').\
                        replace('  ', ' ').replace(' ', '\t')
        values = values.split('\t')
        values.remove('')
        # Fill dictionnary for each amino acid with the corresponding value
        for idx, value in enumerate(values):
            aaindex_values[AMINO_ACID_LIST[idx]] = float(value)
        print(aaindex + '\t' + name)

    # Close navigator
    br.quit()
    end = time.time()
    print('Running time = {} s'.format(end - start))
    return aaindex_dico, aaindex_name

def retrieve_loop_data(path, ext):
    '''
    Function to retrieve data from alignements files with ext extension.

    ARGUMENTS:
    - path : Directory path of files.
    - ext  : Name of the alignment extension (.plain, .txt).
    '''
    # Dictionnaries and lists to store data
    files_remaining = []
    alignment_dict = {}
    dict_number_of_seq = {}
    alignment_dict['number_of_seq'] = dict_number_of_seq
    # Loop to find files with ext extension
    for file in glob.glob('{}{}'.format(path, ext)):
        alignment_list = []
        name_list = []
        loops_idx = {}
        # Read alignment as a dataframe with two columns : name, seq
        df = pd.read_csv(file,
                         sep = '\t',
                         header = None,
                         encoding = 'utf-8',
                         names = ['name', 'seq'])
        # Transform dataframe into a dictionnary with key = name, value = seq
        df_dict = pd.Series(df.seq.values, index=df.name).to_dict()
        for key in df_dict:
            # Storing alignment sequences in a list
            alignment_list.append(df_dict[key])
            # Storing alignment name in a list
            name_list.append(key)
        # Setting a sequence test as a reference
        test_seq = alignment_list[0]
        # Setting number of cysteines found in all alignment sequences
        number_of_cys = 0
        # Setting loop number
        loop_number = 1
        for idx, aa in enumerate(test_seq):
            # Every time we find a cysteine, checking if all other sequences
            # in the alignment have this cysteine
            if aa == 'C':
                not_a_cysteine = False
                for alig in alignment_list:
                    if (alig[idx] != 'C') and (not_a_cysteine == False):
                        not_a_cysteine = True
                # if all sequences have the same cysteine
                # => Setting start and end of loop
                if not_a_cysteine == False:
                    if loop_number == 1:
                        loops_idx['{}_start'.format(loop_number)] = idx
                    elif (loop_number > 1) :
                        loops_idx['{}_end'.format(int(loop_number - 1))] = int(idx - 1)
                        loops_idx['{}_start'.format(loop_number)] = idx
                        storing_end = idx
                    # if new cysteine found, increase loop number
                    loop_number += 1
                    number_of_cys += 1
        loops_idx['{}_end'.format(int(loop_number - 2))] = storing_end
        loops_idx.pop('{}_start'.format(loop_number - 1), None)
        # if number of cyteines is greater than 6, can't process alignments
        if number_of_cys > 6:
            files_remaining.append(file.replace(path , ''))
        elif number_of_cys == 6:
            dict_alig = {}
            family_name = file.replace(ext.replace('*', ''), '').replace(path , '')
            alignment_dict[family_name] = dict_alig
            dict_number_of_seq[family_name] = len(alignment_list)
            for i in range(len(alignment_list)):
                dict_alig_loop = {}
                dict_alig[name_list[i]] = dict_alig_loop
                test = False
                for j in range(1, 6):
                    loop_start = int(loops_idx['{}_start'.format(j)])
                    loop_end = int(loops_idx['{}_end'.format(j)] + 1)
                    # Storing loop sequence (from a cysteine to aa before cysteine) without gaps
                    dict_alig_loop['loop{}'.format(j)] =\
                        alignment_list[i][loop_start:loop_end].replace('-', '')
    # Returning data
    return alignment_dict, files_remaining

def get_min_and_max(df_aaindex):
    '''
    Function to retrieve dictionnary index min and max value to set y axis limits in plots.

    ARGUMENTS:
    - df_aaindex : Pandas dataframe with all values corresponding to an aaindex.
    '''
    # Getting min value
    min_value = min(list(df_aaindex['Value']))
    # Keeping 80% of the min value
    min_value = math.floor((math.floor(min_value)) * 0.80)
    # Getting max value
    max_value = max(list(df_aaindex['Value']))
    # Keeping 120% of the max value
    max_value = math.ceil((math.ceil(max_value)) * 1.20)
    return min_value, max_value

def make_dataframe_from_aaindex_data(web_aaindex_dict, aaindex_names, alignment_data, list_aaindex):
    '''
    Function to make a dataframe from aaindex data.

    ARGUMENTS:
    - web_aaindex_dict : AAindex dictionnary obtained via genome.jp.
    - aaindex_names    : Dictionnary with key = aaindex, value = name/description of the aaindex
    - alignment_data   : Dictionnary with alignment data.
    - list_aaindex     : List of aaaindex.
    '''
    # Lists with variables for pandas dataframe
    aaindex_value_list = []
    loop_value_list = []
    family_value_list = []
    sequence_name_value_list = []
    value_value_list = []
    # Dictionnary coding table for family name
    coding_family_name = {}
    # First index for family name coding table
    idx_family = 1
    for family in list(sorted(list(alignment_data.keys()))):
        if not family.startswith('number_of_seq'):
            print('{:-^60s}'.format(family))
            index_dict = {}
            for id_loop in range(1, 6):
                for sequence_name in alignment_data[family]:
                    seq_loop = (alignment_data[family][sequence_name]['loop{}'.format(id_loop)])
                    ################
                    # SEQUENCE SIZE #
                    ################
                    # Fill family name coding table (setting each family to an index)
                    if family not in coding_family_name:
                        coding_family_name[family] = idx_family
                        idx_family += 1
                    family_value_list.append(coding_family_name[family])
                    aaindex_value_list.append('Sequence_size')
                    loop_value_list.append(id_loop)
                    sequence_name_value_list.append(sequence_name)
                    value_value_list.append(len(seq_loop))
                    ############
                    # AAINDEXs #
                    ############
                    # Loop to retrieve data from dictionnaries
                    for index in list_aaindex:
                        if not index in index_dict:
                            temp_list_to_store_values = []
                            loop_dict = {}
                            index_dict[index] = loop_dict
                            if not 'loop{}'.format(id_loop) in loop_dict:
                                loop_dict['loop{}'.format(id_loop)] = temp_list_to_store_values
                        if not 'loop{}'.format(id_loop) in index_dict[index]:
                            temp_list_to_store_values = []
                            index_dict[index]['loop{}'.format(id_loop)] = temp_list_to_store_values
                        link_to_temp_list_to_store_values = index_dict[index]['loop{}'.format(id_loop)]
                        list_with_index_values = []
                        # For each amino acids, get aaindex value
                        for aa in seq_loop:
                            if aa in AMINO_ACID_LIST:
                                list_with_index_values.append(web_aaindex_dict[index][aa])
                        median_loop = np.median(list_with_index_values)
                        aaindex_value_list.append(index)
                        loop_value_list.append(id_loop)
                        family_value_list.append(coding_family_name[family])
                        sequence_name_value_list.append(sequence_name)
                        value_value_list.append(float(format(median_loop, '.2f')))
            # Delete dictionnary to free memory space
            del index_dict
    # Add 'Sequence_size' at the beginning of the list of aaindex
    list_aaindex.insert(0, 'Sequence_size')
    # Dictionnary with data list
    data = {'AAindex' : aaindex_value_list,
            'Loop' : loop_value_list,
            'Family' : family_value_list,
            'Sequence_name' : sequence_name_value_list,
            'Value' : value_value_list}
    # Pandas dataframe from dictionary
    dataframe = pd.DataFrame(data)
    return data, dataframe, coding_family_name, list_aaindex

def make_aaindex_boxplots_by_family(family_name, list_aaindex, dataframe, coding_family_name, alignment_data, aaindex_names):
    '''
    Function to make a dataframe from aaindex data.

    ARGUMENTS:
    - family_name        : Name of the family you want to plot data.
    - list_aaindex       : List of aaaindex.
    - dataframe          : Pandas dataframe with all value data (see make_dataframe_from_aaindex_data() function).
    - coding_family_name : Dictionnary coding table for family name with key = family name, value = index.
    - alignment_data     : Dictionnary with alignment data.
    - aaindex_names      : Dictionnary with key = aaindex, value = name/description of the aaindex.
    '''
    # Figure Parameters
    HEIGHT = 2100 * 1.5
    WIDTH = 1270 * 1.5
    family_df = dataframe.loc[dataframe['Family'] == coding_family_name[family_name]]
    n_seq_family = alignment_data['number_of_seq'][family_name]
    print('{} sequences in {} family'.format(n_seq_family, family_name))
    # If more than 20 plots, will create multiple figures
    plots_left_to_do = len(list_aaindex) + 1
    idx_fig = 1
    max_number_of_plots_on_fig = 12
    id_list = 0
    last_id_list = 0
    # when you see +1 or -1, it is because there is 1 additional plot for "number of sequences"
    # that isn't in the list_aaindex that need to be taken in account
    while plots_left_to_do > 0:
        if plots_left_to_do > max_number_of_plots_on_fig:
            nb_plots = max_number_of_plots_on_fig
            plots_left_to_do = plots_left_to_do - max_number_of_plots_on_fig + 1
            last_id_list += max_number_of_plots_on_fig - 1
        else:
            nb_plots = plots_left_to_do
            plots_left_to_do = 0
            last_id_list += nb_plots - 1
        # Setting number of rows and columns for matplotib figure
        ncols_grid = math.ceil(math.sqrt(nb_plots))
        if nb_plots <= (ncols_grid * (ncols_grid - 1)):
            nrows_grid = ncols_grid - 1
        else:
            nrows_grid = ncols_grid
        # Creating grid in figure where each cell will be a different plot
        gs = gridspec.GridSpec(nrows_grid, ncols_grid)
        # Making empty figure
        fig = plt.figure(figsize=(12, 7))
        # Setting figure title
        fig.suptitle(family_name + ' Knottins family', fontsize=15)
        ### First plot will be an histogram with number of sequences data
        # Plot will be in grid cell [0 ; 0]
        ax = plt.subplot(gs[0, 0])
        # Making histogram
        ax.hist(list(alignment_data['number_of_seq'].values()), bins=200)
        # Making vertical line corresponding to the actual number of sequences in the family concerned
        ax.axvline(n_seq_family, color='red', linestyle='dashed', linewidth=1)
        # Setting title
        ax.set_title('Histogram : Number of sequences \nby family', fontsize=13)
        # Adding plot to figure
        fig.add_subplot(ax)
        # Setting last row cell index and column cell index used
        id_row = 0
        id_col = 0

        while id_list < last_id_list:
            aaindex = list_aaindex[id_list]
            boxplot_aaindex = []
            # Selecting only data from family dataframe based on aaindex concerned
            df_aaindex_family = family_df.loc[family_df['AAindex'] == aaindex]
            # Selecting only data from full dataframe based on aaindex concerned
            df_aaindex_all_family = dataframe.loc[dataframe['AAindex'] == aaindex]
            for id_loop in range(1, 6):
                # Selecting only data from family + aaindex dataframe based loop concerned
                df_aaindex_by_loop = df_aaindex_family.loc[df_aaindex_family['Loop'] == id_loop]
                # Creating list with data from family + aaindex + loop dataframe
                # Each list concerned values from 1 family, 1 aaindex and 1 loop
                boxplot_aaindex.append(list(df_aaindex_by_loop['Value']))
            # Setting row cell index and column cell index depending on last available cells in grid
            # How ?:
            # Starts with :
            # => row 0 col 1 => row 0 col 2 => ... =>row 0 col (ncols_grid - 1)
            # Then:
            # => row 1 col 0 => row 0 col 1 => ... => row 1 col (ncols_grid - 1)
            # Then:
            # => row (nrows_grid - 1) col 0 => row (nrows_grid - 1) col 1 => ... => row (nrows_grid - 1) col (ncols_grid - 1)
            if id_col < (ncols_grid - 1):
                id_col += 1
                ax2 = plt.subplot(gs[id_row, id_col])
            elif (id_col == (ncols_grid - 1)):
                id_row += 1
                id_col = 0
                ax2 = plt.subplot(gs[id_row, id_col])
            # all aaindex except 'Sequence_size'
            if aaindex in aaindex_names:
                name = aaindex_names[aaindex]
                # Formatting the name/description of the AAINDEX
                # if length of name/description is greater than 30, will cut the name/description in 2 lines (with \n)
                if len(name) > 30:
                    stop = False
                    new_name = ''
                    for letter in range(25, len(name)):
                        if (name[letter] == ' ') and (stop == False) :
                            new_name = name[:(letter)] + '\n' + name[(letter + 1):]
                            stop = True
                    name = new_name
                # Getting max and min values
                min_aaindex_value, max_aaindex_value = get_min_and_max(df_aaindex_all_family)
                # Setting max and min values as y axis limits
                ax2.set_ylim([min_aaindex_value, max_aaindex_value])
            else:
                name = ''
                ax2.set_ylim([0, 40])
            # Setting plot title
            ax2.set_title(aaindex + '\n' + name, fontsize=13)
            # Making boxplot
            ax2.boxplot(boxplot_aaindex, showmeans=True)
            # Setting x tick labels
            ax2.set_xticklabels(['Loop1', 'Loop2', 'Loop3', 'Loop4', 'Loop5'], fontsize=7)
            # Adding plot to figure
            fig.add_subplot(ax2)
            id_list += 1
        # Setting up figures parameters
        plt.subplots_adjust(hspace=0.75)
        DPI = (fig.get_dpi()) * 2
        fig.set_size_inches(float(HEIGHT)/float(DPI),float(WIDTH)/float(DPI))
        # Saving figure
        fig.savefig('../results/AAindex_boxplots/AAindex_Data_for_{}_part{}'.format(family_name, idx_fig))
        plt.close(fig)
        idx_fig += 1
