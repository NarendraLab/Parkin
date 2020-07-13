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
files = 401
while files <= 487:
    folderpath = "/data/LNG/CORNELIS_TEMP/WILL_PRKN/output/" + str(files)
    caselist = os.listdir(folderpath)
    for item in caselist:
        patientpath =  ("/data/LNG/CORNELIS_TEMP/WILL_PRKN/output/" + str(files) +"/"+item)
        try: #Reads one microarray data file. This will be changed into a loop to do this iteratively
            print("Reading data for " + item)
            park2 = pd.read_csv(patientpath, sep='\t', names = ['Chr', 'RsID', 'Family',"Position","Ref","Alt","Log_R_Ratio","B_Allele_Freq"])  
            park2.columns = park2.columns.str.replace(' ', '_')
            print(patientpath)
        except FileNotFoundError:
            text = "no sample"  
    #The neurochip data is subsetted out for each exon and a bit of the surrounding intronic regions -+ 20000BP

        exome = []
        exon1 = park2[(park2.Position < 163148803) & (park2.Position > 163148694)]
        exon2 = park2[(park2.Position < 162864505) & (park2.Position > 162864342)]
        exon3 = park2[(park2.Position < 162683797) & (park2.Position > 162683557)]
        exon4 = park2[(park2.Position < 162622284) & (park2.Position > 162622163)]
        exon5 = park2[(park2.Position < 162475206) & (park2.Position > 162475123)]
        exon6 = park2[(park2.Position < 162394449) & (park2.Position > 162394334)]
        exon7 = park2[(park2.Position < 162206940) & (park2.Position > 162206804)]
        exon8 = park2[(park2.Position < 161990448) & (park2.Position > 161990387)]
        exon9 = park2[(park2.Position < 161970035) & (park2.Position > 161969886)]
        exon10 = park2[(park2.Position < 161807909) & (park2.Position > 161807826)]
        exon11 = park2[(park2.Position < 161781237) & (park2.Position > 161781120)]
        exon12 = park2[(park2.Position < 161771243) & (park2.Position > 161768452)] 
        exomel = []

        exome = {'ex1': exon1,'ex2': exon2,'ex3': exon3,'ex4': exon4,'ex5': exon5,'ex6': exon6,'ex7': exon7,'ex8': exon8,'ex9': exon9,'ex10': exon10,'ex11': exon11,'ex12': exon12
            } #the subsetted data is entered into a dictionary for easy access and manipulation
 

        baseline = park2.Log_R_Ratio.median()
        del1 = baseline - park2.Log_R_Ratio.std()*2 
        exons = 12
        count = 0
        dupexons2 = []
        dupexons3 = []
        dupexons3.append(item)

        while(count < 12):

            count +=1
            dup = False
            dupexons = []
            for items in exome['ex'+ str(count)].Log_R_Ratio.iteritems():
                if (items[1] < del1):
                    dup = True
                    if dup == True:
                        dupexons.append(count)
            sumdupexons = sum(dupexons)/count
            dupexons2.append(sumdupexons)

    	maxdupexons2 = max(dupexons2)
        dupexons2.append(maxdupexons2)
        dupexons3 = dupexons3 + dupexons2
        if maxdupexons2 > 2:
            results.loc[len(results)] = dupexons3

    files+=1
    print(files)

results.to_csv('/data/LNG/CORNELIS_TEMP/WILL_PRKN/DelV1ukboutput.csv')
