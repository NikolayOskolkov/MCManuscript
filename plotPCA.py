# Run this script as: python plotPCA.py

import os; import re
import matplotlib.pyplot as plt
import pandas as pd; import numpy as np
from sklearn.decomposition import PCA; import matplotlib.patches as mpatches

import seaborn as sns
sns.set(font_scale = 1.5)

# Read the matrix of pairwise Mash distances
os.chdir('/home/nikolay/WABI/T_van_der_Valk/IGV_PhyloNorway')
df = pd.read_csv('dist_phylonorway_endo_bact.txt.gz', compression='gzip', sep = '\t')
labels = ['red' if re.search("bact", s) else 'blue' for s in df.columns]
print(df.iloc[0:5, 0:5])
print(labels[0:5])

# Compute PCA and plot PC1 vs. PC2
X = np.log10(df + 1)
X_reduced = PCA(n_components = 2).fit_transform(X)

figure = plt.figure(figsize = (20, 15))
plt.scatter(X_reduced[:, 0], X_reduced[:, 1], c = labels, s = 10)
plt.title('PCA: PHYLONORWAY', fontsize = 25)
plt.xlabel('PC1', fontsize = 22); plt.ylabel('PC2', fontsize = 22)
blue_patch = mpatches.Patch(color='blue', label='Endogenous regions')
red_patch = mpatches.Patch(color='red', label='Microbial-like regions')
plt.legend(handles=[blue_patch, red_patch], fontsize = 20)
plt.show()
