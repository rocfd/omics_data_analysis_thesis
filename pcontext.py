#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:16:13 2020

Description: Script to extract the flanking regions of an epitope given
the epitope sequence and the parental protein sequence.

Output: 30AA pre and post epitope sequences.

@author: rocfarriolduran
"""
# imports
import sys
import csv
import os.path
from pathlib import Path

h = '''
    To right usage of this script:
        $ python3 pcontext.py
    The files in use have to be provided by stating "location/file_name.csv"
    to the input questions that appear in the console.

    <file_name> should be a .csv separated by ";"

    OPTIONAL you can provide coma separated goals
        For example: 0,16,32,64
        If you provide the goals....????????????

    The script returns a <file_name>_out.xlsx as output.

    You can also do the following (if using Mac/Linux OS):
        $ chmod +x pcontext.py
        $ ./pcontext.py <location/file_name.csv> <comaSeparateGoals>

    You need python 3 installed in your computer!!!

    '''



def fileExist(file):
    if file!="":
        if Path(file).is_file():
            return True
        elif file=='-h':
            print(h)
        else:
            print("\n>>>>>>> "+file + " File not exist or is not accessible\n")
            return False
    else:
        return False



data_file=""

num_args=len(sys.argv)


if num_args>=2:
    data_file = sys.argv[1]
if num_args>=3:
    print("Arguments not used.")


while not fileExist(data_file):
    data_file = input(''' Provide location and name of the file dataset.
        For example: data/file_test.csv
            (-h for help)
        Here: ''')
print("Input file: " + data_file)

#data_file= 'test_protein_parser.csv'

extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=data_file[:-lenextension]+"_Out"+extension

out_file=open(nameOutFile, "w")
out_file.writelines("Id;Mut;WT;Anchor;nonAnchor")

f=open(data_file, "r", encoding = 'utf-8-sig')
inputFile = csv.reader(f, delimiter=';')
y=0
for reg in inputFile :
    if y!=0:
        id=reg[0]
        mut=reg[1]
        wt=reg[2]
        anchor=reg[3]
        nonAnchor=reg[4]
        for aa in reg[1]:
            if aa ==
        post=""
        llpep=len(pep)
        start=seq.find(pep)
        if start!=-1:
            if start<=30 :
                if start==0:
                    pre=""
                else:
                    pre=seq[:start-1]
            else:
                pre=seq[start-31:start-1]
            post=seq[start+llpep+1:start+llpep+31]
            out_file.writelines(id+";"+seq+";"+pep+";"+str(start)+";"+pre+";"+str(len(pre))+";"+post+";"+str(len(post))+";"+str(len(pep))+";"+str(len(seq))+"\n")
        else:
            pre=""
            post=""
            out_file.writelines(id+";"+seq+";"+pep+";"+str(start)+";"+pre+";"+str(len(pre))+";"+post+";"+str(len(post))+";"+str(len(pep))+";"+str(len(seq))+"\n")
    y=y+1

out_file.close()
print("Output file: "+nameOutFile)

conda install difflib
import difflib

case_a = 'afrykbnerskojęzyczny'
case_b = 'afrykanerskojęzycznym'

output_list = [li for li in difflib.ndiff(case_a, case_b) if li[0] != ' ']
