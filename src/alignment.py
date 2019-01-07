#! /usr/bin/env python3

import os
import sys
import glob
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait


for file in sorted(glob.glob('../data/KNOTTIN2017Aug30/alignments2017Aug30/*.msf')):
    alignment = ''
    with open(file, 'r') as al:
        alignment = al.readlines()

    new_file = file.replace('KNOTTIN2017Aug30/alignments2017Aug30/', 'alignment_txt/').replace('.msf', '.txt')
    print(new_file)
    br = webdriver.Chrome('../chromedriver')
    br.implicitly_wait(120)

    br.get('https://www.ebi.ac.uk/Tools/msa/mview/')

    br.execute_script("window.scrollTo(0, 30)")

    file = br.find_element_by_xpath('//*[@id="sequence"]')
    file.clear()
    for line in alignment:
        file.send_keys(line)

    br.execute_script("window.scrollTo(0, 700)")
    button1 = br.find_element_by_xpath('//*[@id="jd_toolSubmissionForm"]/div[3]/fieldset/p/a[@class="jd_button shadow"]')
    button1.click()


    select_button = Select(br.find_element_by_xpath('//*[@id="outputformat"]'))
    select_button.select_by_value("plain")

    br.execute_script("window.scrollTo(0, 1080)")

    submit = br.find_element_by_xpath('//*[@id="jd_submitButtonPanel"]/input[@value="Submit"]')
    submit.click()

    ddl_button = br.find_element_by_id('alnFile')
    ddl_button.click()

    alignment_text = br.find_element_by_xpath('/html/body/pre')
    plain_alignment = alignment_text.text

    with open(new_file, 'w') as out:
        out.write(plain_alignment)

    br.quit()
