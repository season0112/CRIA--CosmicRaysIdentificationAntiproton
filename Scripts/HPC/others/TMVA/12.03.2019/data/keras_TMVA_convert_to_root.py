import ROOT
from ROOT import TMVA, TFile, TTree, TCut, TString
import matplotlib.pyplot as plt
import numpy as np
from root_numpy import array2tree, array2root

inputlist = np.arange(0,4)
inputlist = np.append(inputlist,np.arange(7,19))


B1042_1and9_n = np.load('B1042_pr.pl1.flux.l1a9.2016000_7.6_all_negative_147_175_pattern_0.npy')
B1042_1and9_n2_80 = B1042_1and9_n[:,inputlist].copy()
B1042_1and9_n2_80 = B1042_1and9_n2_80[0:int(B1042_1and9_n2_80.shape[0]*0.8),:]
sign_B1042_1and9_n2_80 = B1042_1and9_n2_80[:,-1]/np.abs(B1042_1and9_n2_80[:,-1])
B1042_1and9_n2_80[:,4] = B1042_1and9_n2_80[:,4] * sign_B1042_1and9_n2_80
B1042_1and9_n2_80[:,5] = B1042_1and9_n2_80[:,5] * sign_B1042_1and9_n2_80
B1042_1and9_n2_80[:,6] = B1042_1and9_n2_80[:,6] * sign_B1042_1and9_n2_80
B1042_1and9_n2_80.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(B1042_1and9_n2_80, 'B1042_1and9_n2_80.root', 'tree', 'recreate')

B1042_1and9_n2_20 = B1042_1and9_n[:,inputlist].copy()
B1042_1and9_n2_20 = B1042_1and9_n2_20[-int(B1042_1and9_n2_20.shape[0]*0.2):,:]
sign_B1042_1and9_n2_20 = B1042_1and9_n2_20[:,-1]/np.abs(B1042_1and9_n2_20[:,-1])
B1042_1and9_n2_20[:,4] = B1042_1and9_n2_20[:,4] * sign_B1042_1and9_n2_20
B1042_1and9_n2_20[:,5] = B1042_1and9_n2_20[:,5] * sign_B1042_1and9_n2_20
B1042_1and9_n2_20[:,6] = B1042_1and9_n2_20[:,6] * sign_B1042_1and9_n2_20
B1042_1and9_n2_20.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(B1042_1and9_n2_20, 'B1042_1and9_n2_20.root', 'tree', 'recreate')


#


B1042_1or9_n = np.load('B1042_pr.pl1.flux.l1o9.2016000_7.6_all_negative_147_175_pattern_0.npy')
B1042_1or9_n2 = B1042_1or9_n[:,inputlist].copy()
sign_B1042_1or9_n2 = B1042_1or9_n2[:,-1]/np.abs(B1042_1or9_n2[:,-1])
B1042_1or9_n2[:,4] = B1042_1or9_n2[:,4] * sign_B1042_1or9_n2
B1042_1or9_n2[:,5] = B1042_1or9_n2[:,5] * sign_B1042_1or9_n2
B1042_1or9_n2[:,6] = B1042_1or9_n2[:,6] * sign_B1042_1or9_n2
B1042_1or9_n2.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(B1042_1or9_n2, 'B1042_1or9_n2.root', 'tree', 'recreate')

#


B1042_anti_p = np.load('B1042_antipr.pl1.1800_7.6_all_positive_147_175_pattern_0.npy')
B1042_anti_p2 = B1042_anti_p[:,inputlist].copy()
sign_B1042_anti_p2 = B1042_anti_p2[:,-1]/np.abs(B1042_anti_p2[:,-1])
B1042_anti_p2[:,4] = B1042_anti_p2[:,4] * sign_B1042_anti_p2
B1042_anti_p2[:,5] = B1042_anti_p2[:,5] * sign_B1042_anti_p2
B1042_anti_p2[:,6] = B1042_anti_p2[:,6] * sign_B1042_anti_p2
B1042_anti_p2.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(B1042_anti_p2, 'B1042_anti_p2.root', 'tree', 'recreate')


#


B1042_1or9_p = np.load('B1042_pr.pl1.flux.l1o9.2016000_7.6_all_positive_147_175_pattern_0.npy')
B1042_1or9_p2 = B1042_1or9_p[:,inputlist].copy()
sign_B1042_1or9_p2 = B1042_1or9_p2[:,-1]/np.abs(B1042_1or9_p2[:,-1])
B1042_1or9_p2[:,4] = B1042_1or9_p2[:,4] * sign_B1042_1or9_p2
B1042_1or9_p2[:,5] = B1042_1or9_p2[:,5] * sign_B1042_1or9_p2
B1042_1or9_p2[:,6] = B1042_1or9_p2[:,6] * sign_B1042_1or9_p2
B1042_1or9_p2.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(B1042_1or9_p2[0:B1042_1and9_n2_80.shape[0],:], 'B1042_1or9_p2_train.root', 'tree', 'recreate')
testnumber = B1042_1and9_n2_20.shape[0] + B1042_1or9_n2.shape[0] + B1042_anti_p2.shape[0]
array2root(B1042_1or9_p2[-testnumber:,:], 'B1042_1or9_p2_test.root', 'tree', 'recreate')








