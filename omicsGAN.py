import os
import sys
import pandas as pd 
import numpy as np 
import argparse
from omics1_ros import omics1
from omics2_ros import omics2
from omics3_ros import omics3
os.chdir('C:/Users/RezaMdSelim/Desktop/Other_project/Ashad_alam/4th_project/colon/new_run/run2')
curdir=os.getcwd()


total_update = int(sys.argv[1]) # 1st argument example 5
omics1_name = sys.argv[2]  # 2nd argument example mRNA.csv
omics2_name = sys.argv[3]  # 3rd argument example miRNA.csv
omics3_name = sys.argv[4]  # 3rd argument example miRNA.csv
adj_file = sys.argv[5]     # 4th argument example bipartite_targetscan_gene.csv
adj_file2 = sys.argv[6]     # 4th argument example bipartite_targetscan_gene.csv
adj_file3 = sys.argv[7]     # 4th argument example bipartite_targetscan_gene.csv
label = sys.argv[8]    # 4th argument example label.csv

omics1_result = []
omics2_result = []
omics3_result = []
for i in range(1, total_update+1):
    omics1_result.append(omics1(i,omics1_name,omics2_name,omics3_name,adj_file,adj_file2,label))
    omics2_result.append(omics2(i,omics1_name,omics2_name,omics3_name,adj_file,adj_file3,label))
    omics3_result.append(omics3(i,omics1_name,omics2_name,omics3_name,adj_file2,adj_file3,label))

omics1_result = np.array(omics1_result)
omics2_result = np.array(omics2_result)
omics3_result = np.array(omics3_result)

x1=pd.DataFrame(np.array(omics1_result)).transpose()
x2=pd.DataFrame(np.array(omics2_result)).transpose()
x3=pd.DataFrame(np.array(omics3_result)).transpose()
x1.to_csv(curdir+'best_omics1_result.inc', sep='\t')
x2.to_csv(curdir+'best_omics2_result.inc', sep='\t')
x3.to_csv(curdir+'best_omics3_result.inc', sep='\t')


keep_mRNA = (np.argsort(np.mean(omics1_result,axis=0))[::-1][0])
keep_miRNA = (np.argsort(np.mean(omics2_result,axis=0))[::-1][0])
keep_meth = (np.argsort(np.mean(omics3_result,axis=0))[::-1][0])

print('Best prediction for omics1: ',omics1_result[keep_mRNA])
print('Best update for omics1: ',keep_mRNA+1)
print('Best prediction for omics2: ',omics2_result[keep_miRNA])
print('Best update for omics2: ',keep_miRNA+1)
print('Best prediction for omics3: ',omics3_result[keep_meth])
print('Best update for omics3: ',keep_meth+1)

for i in range(1, total_update+1):
	if i!=keep_mRNA+1:
		os.remove(curdir+"omics1_"+str(i)+".csv") 

	if i!=keep_miRNA+1:
		os.remove(curdir+"omics2_"+str(i)+".csv")
        
	if i!=keep_meth+1:
		os.remove(curdir+"omics3_"+str(i)+".csv")

