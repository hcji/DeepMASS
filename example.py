# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:53:39 2018

@author: jihon
"""

# search formula
from DeepMASS import search_formula
formulas = search_formula(146.105, ppm = 50)

# predict structure
from DeepMASS import run_one_example
mass = 146.105
formula = 'C6H14N2O2'
spectrum_file = 'experiment/spectra/Lysine.csv'
result, neighbors = run_one_example(mass, formula, spectrum_file, energy='40V', thres = 0.5, database='structureDB')