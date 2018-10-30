# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 10:27:13 2018

@author: jihon
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from DeepMASS import compare_structure

from sklearn.metrics import roc_curve, auc

kegg_pairs = pd.read_csv('data/kegg_pairs_refine.csv')
random_pairs = pd.read_csv('data/random_pairs_simulated.csv')

sim_kegg = []
sim_random = []

for i in kegg_pairs.index:
    smi1 = kegg_pairs['smiles1'][i]
    smi2 = kegg_pairs['smiles2'][i]
    similarity = compare_structure(smi1, smi2)
    if similarity > 0:
        sim_kegg.append(similarity)
    
for i in random_pairs.index:
    smi1 = random_pairs['smiles1'][i]
    smi2 = random_pairs['smiles2'][i]
    similarity = compare_structure(smi1, smi2)
    if similarity > 0:
        sim_random.append(similarity)

box_data = pd.DataFrame([sim_kegg, sim_random])    
box_data = box_data.transpose()    
        
sims = np.append(np.array(sim_kegg), np.array(sim_random))
rela = np.append(np.ones(len(sim_kegg)), np.zeros(len(sim_random)))
tpr, fdr, _ = roc_curve(rela, sims)
auroc = auc(tpr, fdr)

# box plot
plt.figure(figsize=(6, 4))
plt.boxplot(box_data.dropna().values)
plt.xticks(range(1,3), ['Kegg Pairs', 'Random Pairs'])
plt.ylabel('Structure Similarity')

# roc plot
plt.figure(figsize=(6, 4))
plt.plot(tpr, fdr, color='black',
         lw=2, label='ROC curve (area = %0.2f)' % auroc)
plt.plot([0, 1], [0, 1], color='red', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.show()
