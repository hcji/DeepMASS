# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 15:57:00 2018

@author: jihon
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit.Chem import MolFromSmiles, Draw


def plot_ms(spectrum):
    plt.figure(figsize=(6, 4))
    plt.vlines(spectrum['mz'], np.zeros(spectrum.shape[0]), np.array(spectrum['intensity']), 'black') 
    plt.axhline(0, color='black')
    plt.show()


def plot_compare_ms(spectrum1, spectrum2, tol=0.05):
    c_mz = []
    c_int = []
    for i in spectrum1.index:
        diffs = abs(spectrum2['mz'] - spectrum1['mz'][i])
        if min(diffs) < tol:
            c_mz.append(spectrum1['mz'][i])
            c_mz.append(spectrum2['mz'][np.argmin(diffs)])
            c_int.append(spectrum1['intensity'][i])
            c_int.append(-spectrum2['intensity'][np.argmin(diffs)])
    c_spec = pd.DataFrame({'mz':c_mz, 'intensity':c_int})   
    
    plt.figure(figsize=(6, 6))
    plt.vlines(spectrum1['mz'], np.zeros(spectrum1.shape[0]), np.array(spectrum1['intensity']), 'gray')
    plt.axhline(0, color='black')
    plt.vlines(spectrum2['mz'], np.zeros(spectrum2.shape[0]), -np.array(spectrum2['intensity']), 'gray')
    plt.vlines(c_spec['mz'], np.zeros(c_spec.shape[0]), c_spec['intensity'], 'red')
    plt.xlabel('m/z')
    plt.ylabel('Relative Intensity')
    plt.show()
    
def visualize_molecule(smiles):
    figs = []
    for smi in smiles:
        figs.append(Draw.MolToImage(MolFromSmiles(smi), size = (160, 160)))
    return figs
    
if __name__ == '__main__':
    from DeepMASS import read_ms
    spectrum1 = read_ms('data/spectra/measured_spectra/40V/C00002.csv')
    spectrum2 = read_ms('data/spectra/measured_spectra/40V/C00008.csv')
    plot_compare_ms(spectrum1, spectrum2)