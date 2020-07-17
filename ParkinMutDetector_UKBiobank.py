#!/usr/bin/python
import pandas as pd
import time
import os
import numpy as np

columns = ['162206852', '162206917', '162394349', '162394435', '162864388', 'mut', 'sample']
results = pd.DataFrame()
results = pd.DataFrame(columns = columns)


#This cell is run as a batch job 
samplelist = list()
caseList = list()

#Intialize results table and exon copy number values

cnv = 0
insert = 0
files = 0
while files <= 487:
    folderpath = "/path/to/PRKN/output/" + str(files)
    caselist = os.listdir(folderpath)
    for item in caselist:
        patientpath =  ("/path/to/PRKN/output/" + str(files) +"/"+item)
        try: #Reads one microarray data file. This will be changed into a loop to do this iteratively
            print("Reading data for " + item)
            park2 = pd.read_csv(patientpath, sep='\t', names = ['Chr', 'RsID', 'Family',"Position","Ref","Alt","Log_R_Ratio","B_Allele_Freq"])  
            park2.columns = park2.columns.str.replace(' ', '_')
            print(patientpath)
        except FileNotFoundError:
            text = "no sample"  
    #The neurochip data is subsetted out for each exon and a bit of the surrounding intronic regions -+ 20000BP

     
        park2_pos = park2.set_index('Position')
        BAF = park2_pos.B_Allele_Freq
        SOI = BAF[[162206852, 162206917, 162394349, 162394435, 162864388]]
        SOI2 = SOI.tolist()

        mut = False
        for x in SOI2:
            if x < .75:
                mut = True

        SOI2.append(mut)

        SOI2.append(item)

        results.loc[len(results)] = SOI2

    files+=1
    print(files)

results2 = results.set_index('sample')

results2.to_csv('/path/to/PRKN/testukboutput3_0_487.csv')
