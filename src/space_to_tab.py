#! /usr/bin/env python3

import os
import sys
import glob
import re
import time

''' Script that converts spaces to tab in plain alignment'''
if __name__ == '__main__':
    for file in glob.glob('../data/alignment_txt/*.txt'):
        with open(file, 'r') as txt_file:
            alignment = txt_file.readlines()
            with open(file.replace('txt', 'plain'), 'w') as out:
                for line in alignment:
                    line = re.sub('\ +', '\t', line)
                    out.write(line)
            os.remove(file)
