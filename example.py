# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:53:39 2018

@author: jihon
"""

# search formula
from DeepMASS import search_formula
formulas = search_formula(146.105, ppm = 50)

# compare isotope pattern
import pandas as pd
from DeepMASS import isotope_pattern, compare_isotope
iso = pd.DataFrame({'mass': [146.105, 107.108], 'intensity': [0.924, 0.05]})
iso1 = isotope_pattern('C6H14N2O2')
score = compare_isotope(iso, iso1)

# predict structure
import time
from DeepMASS import run_one_example
start_time = time.time()
mass = 146.105
formula = 'C6H14N2O2'
spectrum_file = 'experiment/spectra/Lysine.csv'
result, neighbors = run_one_example(mass, formula, spectrum_file, energy='40V', thres = 0.5, database='structureDB')
end_time = time.time()
print('The running time is ' + str(round(end_time-start_time, 2)) + 's')