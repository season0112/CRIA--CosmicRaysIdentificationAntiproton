## The purpose of this scipt is to test ACQT7.9 and ACQT7.8 data. (07.2021)
# Major differences is :
# 1. MC TOFBeta
# 2. Hit Z coordinate
# 3. GBL Tracker Fit

import numpy as np
from root_numpy import root2array, tree2array
from ROOT import TFile
import matplotlib.pyplot as plt

'''
pass78_Negative = root2array("/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/pass7.8/negative/results/merged.root", "AntiprotonLowEnergyTree")
pass78_Positive = root2array("/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/pass7.8/positive/results/merged.root", "AntiprotonLowEnergyTree")
pass79_Negative = root2array("/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/pass7.8/negative/results/merged.root", "AntiprotonLowEnergyTree")
pass79_Positive = root2array("/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/pass7.8/positive/results/merged.root", "AntiprotonLowEnergyTree")
'''

pass78_Negative = np.load("/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/pass7.8/negative/results/merged.npy", allow_pickle=True)
pass79_Negative = np.load("/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/pass7.9/negative/results/merged.npy", allow_pickle=True)

pass78negative_L1Charge = []
pass78negative_L2Charge = []
pass78negative_L9Charge = []
pass78negative_TOFClusterEnergy1 = []
pass78negative_TOFClusterEnergy2 = []
pass78negative_TOFClusterEnergy3 = []
pass78negative_TOFClusterEnergy4 = []
for i in range(pass78_Negative.shape[0]):
    pass78negative_L1Charge.append(pass78_Negative["TrackerCharges"][i][0])
    pass78negative_L2Charge.append(pass78_Negative["TrackerCharges"][i][1])
    pass78negative_L9Charge.append(pass78_Negative["TrackerCharges"][i][8])
    if (pass78_Negative["TOFClusterEnergy"][i].shape[0]==4):
      pass78negative_TOFClusterEnergy1.append(pass78_Negative["TOFClusterEnergy"][i][0])
      pass78negative_TOFClusterEnergy2.append(pass78_Negative["TOFClusterEnergy"][i][1])
      pass78negative_TOFClusterEnergy3.append(pass78_Negative["TOFClusterEnergy"][i][2])
      pass78negative_TOFClusterEnergy4.append(pass78_Negative["TOFClusterEnergy"][i][3])
pass78negative_L1Charge = np.array(pass78negative_L1Charge)
pass78negative_L2Charge = np.array(pass78negative_L2Charge)
pass78negative_L9Charge = np.array(pass78negative_L9Charge)
pass78negative_TOFClusterEnergy1 = np.array(pass78negative_TOFClusterEnergy1)
pass78negative_TOFClusterEnergy2 = np.array(pass78negative_TOFClusterEnergy2)
pass78negative_TOFClusterEnergy3 = np.array(pass78negative_TOFClusterEnergy3)
pass78negative_TOFClusterEnergy4 = np.array(pass78negative_TOFClusterEnergy4)


pass79negative_L1Charge = []
pass79negative_L2Charge = []
pass79negative_L9Charge = []
pass79negative_TOFClusterEnergy1 = []
pass79negative_TOFClusterEnergy2 = []
pass79negative_TOFClusterEnergy3 = []
pass79negative_TOFClusterEnergy4 = []
for i in range(pass79_Negative.shape[0]):
    pass79negative_L1Charge.append(pass79_Negative["TrackerCharges"][i][0])
    pass79negative_L2Charge.append(pass79_Negative["TrackerCharges"][i][1])
    pass79negative_L9Charge.append(pass79_Negative["TrackerCharges"][i][8])
    if (pass79_Negative["TOFClusterEnergy"][i].shape[0]==4):
      pass79negative_TOFClusterEnergy1.append(pass79_Negative["TOFClusterEnergy"][i][0])
      pass79negative_TOFClusterEnergy2.append(pass79_Negative["TOFClusterEnergy"][i][1])
      pass79negative_TOFClusterEnergy3.append(pass79_Negative["TOFClusterEnergy"][i][2])
      pass79negative_TOFClusterEnergy4.append(pass79_Negative["TOFClusterEnergy"][i][3])
pass79negative_L1Charge = np.array(pass79negative_L1Charge)
pass79negative_L2Charge = np.array(pass79negative_L2Charge)
pass79negative_L9Charge = np.array(pass79negative_L9Charge)
pass79negative_TOFClusterEnergy1 = np.array(pass79negative_TOFClusterEnergy1)
pass79negative_TOFClusterEnergy2 = np.array(pass79negative_TOFClusterEnergy2)
pass79negative_TOFClusterEnergy3 = np.array(pass79negative_TOFClusterEnergy3)
pass79negative_TOFClusterEnergy4 = np.array(pass79negative_TOFClusterEnergy4)





varname = "UpperTofCharge"

plt.figure(figsize=(9,9))
plt.hist(pass78_Negative[varname],bins=200, range=(0.8, 1.5), color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79_Negative[varname],bins=200, range=(0.8, 1.5), color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel(varname,fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/'+ varname +'.png')

## Tracker Layer Charges
plt.figure(figsize=(9,9))
plt.hist(pass78negative_L1Charge, bins=200, range=(0.6, 1.8),  color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79negative_L1Charge, bins=200, range=(0.6, 1.8),  color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel('L1Charge',fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/L1Charge.png')

plt.figure(figsize=(9,9))
plt.hist(pass78negative_L2Charge, bins=200, range=(0.6, 1.8),  color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79negative_L2Charge, bins=200, range=(0.6, 1.8),  color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel('L2Charge',fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/L2Charge.png')

plt.figure(figsize=(9,9))
plt.hist(pass78negative_L9Charge, bins=200, range=(0.6, 1.8),  color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79negative_L9Charge, bins=200, range=(0.6, 1.8),  color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel('L9Charge',fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/L9Charge.png')

## TOFClusterEnergy
plt.figure(figsize=(9,9))
plt.hist(pass78negative_TOFClusterEnergy1, bins=200, range=(1, 5),  color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79negative_TOFClusterEnergy1, bins=200, range=(1, 5),  color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel('TOFClusterEnergy1',fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/TOFClusterEnergy1.png')

plt.figure(figsize=(9,9))
plt.hist(pass78negative_TOFClusterEnergy2, bins=200, range=(1, 5),  color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79negative_TOFClusterEnergy2, bins=200, range=(1, 5),  color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel('TOFClusterEnergy2',fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/TOFClusterEnergy2.png')

plt.figure(figsize=(9,9))
plt.hist(pass78negative_TOFClusterEnergy3, bins=200, range=(1, 5),  color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79negative_TOFClusterEnergy3, bins=200, range=(1, 5),  color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel('TOFClusterEnergy3',fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/TOFClusterEnergy3.png')

plt.figure(figsize=(9,9))
plt.hist(pass78negative_TOFClusterEnergy4, bins=200, range=(1, 5),  color="b", facecolor='b',edgecolor='b', alpha=1.0, label="Pass78_Negative")
plt.hist(pass79negative_TOFClusterEnergy4, bins=200, range=(1, 5),  color="r", facecolor='r',edgecolor='r', alpha=0.5, label="Pass79_Negative")
plt.legend(loc='best',fontsize=15)
plt.xlabel('TOFClusterEnergy4',fontsize=22)
plt.savefig('/hpcwork/jara0052/sichen/pass7.8_pass7.9_CompareTest/TOFClusterEnergy4.png')








