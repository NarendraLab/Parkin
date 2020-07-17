#!/usr/bin/python
import pandas as pd
import time
import os
import numpy as np

columns = ['sample','E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12','max']
results = pd.DataFrame()
results = pd.DataFrame(columns = columns)


#This cell is run as a batch job 
samplelist = list()
caseList = list()

#Intialize results table and exon copy number values

cnv = 0
insert = 0
files = 0
while files <= 1:
    folderpath = "path/to/PRKN/output/" + str(files)
    caselist = os.listdir(folderpath)
    for item in caselist:
        patientpath =  ("/path/to/PRKN/" + str(files) +"/"+item)
        try: #Reads one microarray data file. This will be changed into a loop to do this iteratively
            print("Reading data for " + item)
            park2 = pd.read_csv(patientpath, sep='\t', names = ['Chr', 'RsID', 'Family',"Position","Ref","Alt","Log_R_Ratio","B_Allele_Freq"])  
            park2.columns = park2.columns.str.replace(' ', '_')
            print(patientpath)
        except FileNotFoundError:
            text = "no sample"  
    #The neurochip data is subsetted out for each exon and all of the surrounding intronic regions

        exomel = []
        e1 = park2[(park2.Position < 173148803) & (park2.Position > 163006601)]
        e2 = park2[(park2.Position < 163006600) & (park2.Position > 162774071)]
        e3 = park2[(park2.Position < 162774070) & (park2.Position > 162652922)]
        e4 = park2[(park2.Position < 162652921) & (park2.Position > 162548686)]
        e5 = park2[(park2.Position < 162548685) & (park2.Position > 162434787)]
        e6 = park2[(park2.Position < 162434786) & (park2.Position > 162300638)]
        e7 = park2[(park2.Position < 162300637) & (park2.Position > 162098627)]
        e8 = park2[(park2.Position < 162098626) & (park2.Position > 161980212)]
        e9 = park2[(park2.Position < 161980211) & (park2.Position > 161888899)]
        e10 = park2[(park2.Position < 161888898) & (park2.Position > 161794533)]
        e11 = park2[(park2.Position < 161794532) & (park2.Position > 161776183)]
        e12 = park2[(park2.Position < 161776182) & (park2.Position > 151768452)]

#the subsetted data is entered into a dictionary for easy access and manipulation

        exomel = {'exx1': e1,'exx2': e2,'exx3': e3,'exx4': e4,'exx5': e5,'exx6': e6,'exx7': e7,'exx8': e8,'exx9': e9,'exx10': e10,'exx11': e11,'exx12': e12
                }  
  
        exons = 12
        count = 0
        dupexons2 = []
        dupexons3 = []
        dupexons3.append(item)

        while(count < 12):

            count +=1
            dup = False
            dupexons = []
            for items in exomel['exx'+str(count)].B_Allele_Freq.iteritems():
                if (items[1] > .125 and items[1] < .375) or (items[1] > .625 and items[1] < .875):
                    dup = True
                    if dup == True:
                        dupexons.append(count)
            sumdupexons = sum(dupexons)/count
            dupexons2.append(sumdupexons)

        maxdupexons2 = max(dupexons2)
        dupexons2.append(maxdupexons2)
        dupexons3 = dupexons3 + dupexons2
        if maxdupexons2 > 1:
            results.loc[len(results)] = dupexons3

    files+=1
    print(files)

results.to_csv('/path/to/PRKN/ParkinDup_ukboutput_0_1.csv')
