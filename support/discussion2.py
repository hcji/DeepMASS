# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:47:57 2018

@author: jihon
"""


import numpy as np
import matplotlib.pyplot as plt
from DeepMASS import read_ms, ms2vec, fft_ccor
from MassPlot import plot_compare_ms

spectrum1 = read_ms('data/spectra/measured_spectra/40V/C00006.csv')
spectrum2 = read_ms('data/spectra/measured_spectra/40V/C00003.csv')
# spectrum3 = spectrum2.copy()
# spectrum3['mz'] = spectrum3['mz'] - 12

plot_compare_ms(spectrum1, spectrum2, shift=1000)
# plot_compare_ms(spectrum1, spectrum3)

vec1 = ms2vec(spectrum1)
vec2 = ms2vec(spectrum2)
ccor = fft_ccor(vec1, vec2, norm=False).A.squeeze()
shifts = np.array(range(len(ccor)))/100 - 1000
keep = (shifts > -10) & (shifts < 40)

plt.figure(figsize=(6, 6))
plt.plot(shifts[keep], ccor[keep])
plt.xlabel('m/z shift')
plt.ylabel('Cross correlation coefficient')
plt.show()

