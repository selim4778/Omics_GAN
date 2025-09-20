import os
import re
import numpy as np
import pandas as pd
from pandas import DataFrame, Series

os.chdir("C:/Users/RezaMdSelim/OneDrive - Tulane University/Desktop/Other_project/Ashad_alam/4th_project/colon/docking")


directory = os.getcwd() + '/pre_protein'
output = os.getcwd() 
protein = os.listdir(directory)


Data = DataFrame()
for receptor in protein:
    curdir = directory + '/'+ receptor + '/log'
    subs = os.listdir(curdir)
    mylines = [] 
    drug = DataFrame()
    value = DataFrame()
    for file in subs: 
        a_file = open (curdir + '/' + file)    
        lines_to_read = [27]
        for position, line in enumerate(a_file):
            if position in lines_to_read:
                #print line
                mylines.append(line)
                x = line.split(" ")
                d = file.split('.')
                y = pd.DataFrame(x)
                #y = y.transpose()
                #y['name']=d[0]
                #d = file.split('.')
                z = y.iloc[11] + y.iloc[12]
                z['drug']=d[0]
                df = pd.DataFrame(z)
                df = df.transpose()            
        value = pd.concat([value, df], axis=0)
        w = value.set_index('drug')
        w.columns = [receptor]
    Data = pd.concat([Data, w], axis=1)  
    
    
Data.to_csv(output +'/'+ 'original_score_colon.csv')



###########  ranking

Data.dtypes
s = Data.astype(float)
s.dtypes

s['drug_score']= s.mean(axis=1, numeric_only=True)   
short_drug = s.sort_values(['drug_score'], ascending=True)
short_drug = short_drug.drop('drug_score', axis=1).T

short_drug['protein_score']= short_drug.mean(axis=1, numeric_only=True)   
short_drug_protein = short_drug.sort_values(['protein_score'], ascending=True)
short_drug_protein = short_drug_protein.drop('protein_score', axis=1)

short_drug_protein.to_csv(output +'/'+ 'original_score_colon_ranked.csv')


columns = short_drug_protein.columns

print(columns)

with open("columns_name.txt", "w") as file:
    file.write(columns)