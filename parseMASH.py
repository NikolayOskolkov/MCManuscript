#!/home/nikolay/miniconda3/bin/python3

import os
import pandas as pd

os.chdir('/home/nikolay/WABI/T_van_der_Valk/IGV_PhyloNorway')

#df = pd.read_csv('mash_dist_hippuris.txt', header = None, sep = '\t')
#df = pd.read_csv('mash_dist_plants_plus_bacteria.txt', header = None, sep = '\t')
df = pd.read_csv('mash_dist_plants_plus_bacteria_plus_hippuris.txt', header = None, sep = '\t')
df.columns = ['RefID', 'QueryID', 'MashDistance', 'Pvalue', 'MatchHash']
print(df.iloc[0:5, 0:5])

unique_RefIDs = list(dict.fromkeys(list(df.QueryID)))

dist_matrix = pd.DataFrame(columns = unique_RefIDs, index = unique_RefIDs)
for i in range(len(unique_RefIDs)):
	print('Processing {} references'.format(i))
	for j in range(len(unique_RefIDs)):
		dist_matrix.iloc[i, j] = float(df[(df.QueryID == unique_RefIDs[i]) & (df.RefID == unique_RefIDs[j])].MashDistance) 

print(dist_matrix.iloc[0:5, 0:5])
#dist_matrix.to_csv('dist_hippuris.txt', sep = '\t')
#dist_matrix.to_csv('dist_plants_plus_bacteria.txt', sep = '\t')
dist_matrix.to_csv('dist_plants_plus_bacteria_plus_hippuris.txt', sep = '\t')

