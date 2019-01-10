#! /usr/bin/env python3

import os
import sys
import glob
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait

'''
Script that convert .msf files to txt files.

REQUIREMENT : Internet
'''

def convert_msf_to_plain(msf_file):
    '''
    Function to make a dataframe from aaindex data.

    ARGUMENTS:
    - msf_file : Path to msf file
    '''
    alignment = ''
    # reading file
    print(msf_file)
    with open(msf_file, 'r') as al:
        alignment = al.readlines()
    new_file = msf_file.replace('KNOTTIN2017Aug30/alignments2017Aug30/', 'alignment_txt/').replace('.msf', '.txt')
    try:
        br = webdriver.Chrome('../chromedriver')
    except:
        print('Bad chromedriver version, please download one specific to your OS')
        sys.exit(1)
    br.implicitly_wait(120)
    # open web page
    br.get('https://www.ebi.ac.uk/Tools/msa/mview/')
    # scroll web page
    br.execute_script("window.scrollTo(0, 30)")
    # find text area
    text_area = br.find_element_by_xpath('//*[@id="sequence"]')
    # clear text area if text written in it
    text_area.clear()
    # send alignment from msf file to text area
    for line in alignment:
        text_area.send_keys(line)
    # scroll towards the bottom of the page
    br.execute_script("window.scrollTo(0, 700)")
    # Find 'more options' button
    button1 = br.find_element_by_xpath('//*[@id="jd_toolSubmissionForm"]/div[3]/fieldset/p/a[@class="jd_button shadow"]')
    # click on button
    button1.click()
    # Find output format windows
    select_button = Select(br.find_element_by_xpath('//*[@id="outputformat"]'))
    # select "plain" format
    select_button.select_by_value("plain")
    # scroll to bottom of the page
    br.execute_script("window.scrollTo(0, 1080)")
    # Find submit button
    submit = br.find_element_by_xpath('//*[@id="jd_submitButtonPanel"]/input[@value="Submit"]')
    # click on submit button
    submit.click()
    # Find download button
    ddl_button = br.find_element_by_id('alnFile')
    # click on 'download alignment' button
    ddl_button.click()
    # find contents in web page
    alignment_text = br.find_element_by_xpath('/html/body/pre')
    # extract text from webpage
    plain_alignment = alignment_text.text
    # save alignment
    with open(new_file, 'w') as out:
        out.write(plain_alignment)
    br.quit()


if __name__ == '__main__':
    for file in sorted(glob.glob('../data/KNOTTIN2017Aug30/alignments2017Aug30/*.msf')):
        convert_msf_to_plain(file)
