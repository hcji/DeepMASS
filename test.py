# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 14:09:55 2018

@author: jihon
"""

from DeepMASS import build_model, fine_tune, leave_one_out_test
build_model(Test=False, Save=True)
fine_tune(Test=False, Save=True)
leave_one_out_test(energy='40V', database='structureDB')