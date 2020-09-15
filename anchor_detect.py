#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:16:13 2020

Description: characterization of mutations within epitopes (neoepitopes) in anchor and non-anchor positions regarding the binding HLA molecule. 

Functionality: given an epitope, its mutated aminoacid and the binding motif of the HLA restricted allotype (where the epitope), the script classifies
the mutation on wether it is located on a anchor position or on a non-anchor position. Importantly, anchor positions differ between HLA allotypes.

@author: rocfarriolduran
"""
# imports
import sys
import csv
import os.path
from pathlib import Path

h = '''

    Heeeeeelp please!!!!

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



data_file="anchorsDat.csv"
anchors_file="anchorsParams.csv"

num_args=len(sys.argv)


if num_args>=2:
    data_file = sys.argv[1]
if num_args==3:
    inputAnchors = sys.argv[2]


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

while not fileExist(anchors_file):
    anchors_file = input(''' Provide location and name of the file params.
        For example: data/file_test.csv
            (-h for help)
        Here: ''')
print("Anchors file: " + anchors_file)

out_file=open(nameOutFile, "w")
out_file.writelines("Id;PeptideA;PeptideB;Equals;NotEquals;Len;Anchor;NotAnchor\n")

f=open(data_file, "r", encoding = 'utf-8-sig')
inputFile = csv.reader(f, delimiter=';')
ff=open(anchors_file, "r", encoding = 'utf-8-sig')
inputAnchors = csv.reader(ff, delimiter=';')
y=0
par=[]
unSortedinputAnchors=[]
for row in inputAnchors:
    if y!=0:
        par.append(int(row[0]))
        par.append(int(row[1]))
        unSortedinputAnchors.append(par)
        par = []
    y=y+1

sortedinputAnchors = sorted(unSortedinputAnchors, key=lambda row:(row[0],row[1]), reverse=False)

lAnt = 0
lList = []
parAnchors = []
listAnchors = []
for row in sortedinputAnchors:
    if lAnt!=row[0]:
        parAnchors = []
        parAnchors.append(lAnt)
        parAnchors.append(lList)
        listAnchors.append(parAnchors)
        lAnt = row[0]
        lList = []
    else:
        lList.append(row[1])
y=0
for reg in inputFile :
    if y!=0:
        id=reg[0]
        pepA=reg[1]
        pepB=reg[2]
        llpepA=len(pepA)
        llpepB=len(pepB)
        Anchor=[]
        NotAnchor=[]
        Equals=""
        NotEquals=""
        PAnchor = []
        if llpepA!=llpepB:
            Len="Not equals"
        else:
            Len=str(llpepA)
            for row in listAnchors:
                if row[0]==llpepA:
                    PAnchor=row[1]
        for i in range(0,llpepA):
            print(str(i)+";"+pepA[i:i+1]+";"+pepB[i:i+1])
            if pepA[i:i+1]==pepB[i:i+1]:
                Equals=Equals+pepA[i:i+1]
                NotEquals=NotEquals+"*"
                Anchor.append(i)
            else:
                NotEquals=NotEquals+pepA[i:i+1]
                Equals=Equals+"*"
                NotAnchor.append(i)
        out_file.writelines(id+";"+pepA+";"+pepB+";"+Equals+";"+NotEquals+";"+Len+";"+str(Anchor)+";"+str(NotAnchor)+"\n")

    y=y+1

out_file.close()
print("Output file: "+nameOutFile)
