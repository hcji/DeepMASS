# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:53:39 2018

@author: jihon
"""

from DeepMASS import run_one_example

mass = 174.112
formula = 'C6H14N4O2'
spectrum_file = 'experiment/spectra/Arginine.csv'

result, neighbors = run_one_example(mass, formula, spectrum_file, energy='40V', thres = 0.5, database='structureDB')
