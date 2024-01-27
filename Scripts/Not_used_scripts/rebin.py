#!/usr/bin/env python
import numpy as np
import argparse
import os
import time

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)

os.chdir('/hpcwork/jara0052/sichen/analysis_7.0/B1130_pass7_nov18_ext_7.7_all/147-1000GeV_v2/transferdata')
path = os.getcwd()

positive_259_330 = np.load('positive/259_330GeV/pattern0/positive_259_330_pattern_0.npy')
positive_330_525 = np.load('positive/330_525GeV/pattern0/positive_330_525_pattern_0.npy')
negative_259_330 = np.load('negative/259_330GeV/pattern0/negative_259_330_pattern_0.npy')
negative_330_525 = np.load('negative/330_525GeV/pattern0/negative_330_525_pattern_0.npy')

positive_259_450 = np.row_stack((positive_259_330,positive_330_525[np.where(positive_330_525[:,-1]<450)[0],:]))
negative_259_450 = np.row_stack((negative_259_330,negative_330_525[np.where(negative_330_525[:,-1]>-450)[0],:]))

if os.path.isdir('positive/259_450GeV/pattern0/'):
    pass
else:
    mkdir('positive/259_450GeV/pattern0/')
if os.path.isdir('negative/259_450GeV/pattern0/'):
    pass
else:
    mkdir('negative/259_450GeV/pattern0/')


np.save(path+'/positive/259_450GeV/pattern0/positive_259_450_pattern_0.npy',positive_259_450)
np.save(path+'/negative/259_450GeV/pattern0/negative_259_450_pattern_0.npy',negative_259_450)



