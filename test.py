# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:09:55 2018

@author: jihon
"""

from DeepMASS import build_model, fine_tune, leave_one_out_test
build_model(Test=True, Save=False)
fine_tune(Test=True, Save=False)
leave_one_out_test(energy='40V', database='structureDB')