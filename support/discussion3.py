# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 14:36:41 2018

@author: jihon
"""

import numpy as np
from DeepMASS import run_one_example, compare_structure

mass = 146.19
formula = 'C6H14N2O2'
spectrum_file = 'experiment/spectra/Lysine.csv'

result, neighbors = run_one_example(mass, formula, spectrum_file, energy='40V', thres = 0.5, database='structureDB')

result['image'][0]
i = 0
neighbors['Image'][np.argsort(-neighbors['similarity'])[i]]
compare_structure(neighbors['smiles'][np.argsort(-neighbors['similarity'])[i]], result['SMILES'][0])