#!/usr/bin/python
import pandas as pd
import time
import os
import numpy as np
park2 = pd.DataFrame()
results = pd.DataFrame()
#This cell is run as a batch job 
samplelist = list()
caseList = list()
#Intialize results table and exon copy number values
columns = ['Sample ID', 'Exon 1','Exon 2','Exon 3','Exon 4','Exon 5','Exon 6','Exon 7','Exon 8','Exon 9','Exon 10','Exon 11','Exon 12']
rows = {'e1': 0, 'e2': 0,'e3': 0,'e4': 0,'e5': 0,'e6': 0,'e7': 0,'e8': 0,'e9': 0,'e10': 0,'e11': 0,'e12': 0}
ncolumns = ['ID','Ex1qc','Ex2qc','Ex3qc','Ex4qc','Ex5qc','Ex6qc','Ex7qc','Ex8qc','Ex9qc','Ex10qc','Ex11qc','Ex12qc'] 
qrows = {'e1': 'NI', 'e2': 'NI','e3': 'NI','e4': 'NI','e5': 'NI','e6': 'NI','e7': 'NI','e8': 'NI','e9': 'NI','e10': 'NI','e11': 'NI','e12': 'NI'}
results = pd.DataFrame()
results = pd.DataFrame(columns = columns)
quc = pd.DataFrame()
quc = pd.DataFrame(columns = ncolumns)
cnv = 0
insert = 0
files = 0
while files <= 488:
    folderpath = "/data/LNG/CORNELIS_TEMP/WILL_PRKN/output/" + str(files)
    caselist = os.listdir(folderpath)
    for item in caselist:
        patientpath =  ("/data/LNG/CORNELIS_TEMP/WILL_PRKN/output/" + str(files) +"/"+item)
        try: #Reads one microarray data file. This will be changed into a loop to do this iteratively
            print("Reading data for " + item)
	    park2 = pd.read_csv(patientpath, sep='\t', header = 0, names = ['Chr', 'RsID', 'Family',"Position","Ref","Alt","LoGR","BAF"]) 
            park2.columns = park2.columns.str.replace(' ', '_')
            print(patientpath)
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
            e1 = park2[(park2.Position < 163168803) & (park2.Position > 163128694)]
            e2 = park2[(park2.Position < 162884505) & (park2.Position > 162844342)]
            e3 = park2[(park2.Position < 162703797) & (park2.Position > 162663557)]
            e4 = park2[(park2.Position < 162642284) & (park2.Position > 162602163)]
            e5 = park2[(park2.Position < 162495206) & (park2.Position > 162455123)]
            e6 = park2[(park2.Position < 162414449) & (park2.Position > 162374334)]
            e7 = park2[(park2.Position < 162226940) & (park2.Position > 162186804)]
            e8 = park2[(park2.Position < 162010448) & (park2.Position > 161970387)]
            e9 = park2[(park2.Position < 161990035) & (park2.Position > 161949886)]
            e10 = park2[(park2.Position < 161827909) & (park2.Position > 161787826)]
            e11 = park2[(park2.Position < 161801237) & (park2.Position > 161761120)]
            e12 = park2[(park2.Position < 161791243) & (park2.Position > 161748452)]
            exome = {'ex1': exon1,'ex2': exon2,'ex3': exon3,'ex4': exon4,'ex5': exon5,'ex6': exon6,'ex7': exon7,'ex8': exon8,'ex9': exon9,'ex10': exon10,'ex11': exon11,'ex12': exon12
                } #the subsetted data is entered into a dictionary for easy access and manipulation
            exomel = {'exx1': e1,'exx2': e2,'exx3': e3,'exx4': e4,'exx5': e5,'exx6': e6,'exx7': e7,'exx8': e8,'exx9': e9,'exx10': e10,'exx11': e11,'exx12': e12
                 }  
            exon2.drop([756,758],axis=0, inplace=True)
            exon7.drop([489,497],axis=0, inplace=True)
            baseline = park2.LoGR.mean()
            del1 = baseline - park2.LoGR.std()
            del2 = baseline - 1
	    in1 = baseline + park2.LoGR.std() *2
            exons = 12
            count = 1
            sv = False 
            while(count < 12):
                count +=1
                intensity = exome['ex'+ str(count)].LoGR.median()
                hetero = False
                for items in exomel['exx'+str(count)].BAF.iteritems():
                    if items[1] > .2 and items[1] < .8:
                        hetero = True
                        break
                if (intensity > del1) & (count != 9) & (intensity < in1):
                    rows['e'+str(count)] = 2
                    qrows['e'+str(count)] = 'P'
                else: 
                    if (intensity <= del1) & (intensity > del2) & (count  != 5) & (count != 9) & (count != 10) & (count != 12) & (hetero == False):
                        rows['e'+str(count)] = 1
                        qrows['e'+str(count)] = 'P'
                        cnv+=1
                        sv = True
                    if (count == 3 or count ==8 or count == 0):
                        qrows['e'+str(count)] = 'F - low certainty'
                    if (intensity <= del2) & (count != 9) & (count != 10):
                        qrows['e'+str(count)] = 'P'
                        rows['e'+str(count)] = 0
                        cnv+=1
                        sv = True
                    if intensity > .8:
                        qrows['e'+str(count)] = 'F'
                        rows['e'+str(count)] = '-'
                    if (count == 1 or count  == 5 or count == 9 or count == 10 or count == 12):
                        rows['e'+str(count)] = '-'
                        qrows['e'+str(count)] = 'NI'
            if sv == True:
            	print("Deletion at " + item)
                results = results.append({'Sample ID': item, 'Exon 1': '-', 'Exon 2': rows['e2'],'Exon 3': rows['e3'],'Exon 4': rows['e4'],'Exon 5': rows['e5'],'Exon 6': rows['e6'],'Exon 7': rows['e7'],'Exon 8': rows['e8'],'Exon 9': '-','Exon 10': rows['e10'],'Exon 11': rows['e11'], 'Exon 12': rows['e12']}, ignore_index=True)
                quc = quc.append({'ID': item, 'Ex1qc': qrows['e1'], 'Ex2qc': qrows['e2'],'Ex3qc': qrows['e3'],'Ex4qc': qrows['e4'],'Ex5qc': qrows['e5'],'Ex6qc': qrows['e6'],'Ex7qc': qrows['e7'],'Ex8qc': qrows['e8'],'Ex9qc': qrows['e9'],'Ex10qc': qrows['e10'],'Ex11qc': qrows['e11'], 'Ex12qc': qrows['e12']}, ignore_index=True)
        except IOError:
            print("no sample")
    files+=1
    print(files)
print('There are a total of ' + str(cnv) + ' separate deletions')
results = results.join(quc)
results.to_csv('/data/LNG/CORNELIS_TEMP/WILL_PRKN/DelV1ukboutput.csv')
 
