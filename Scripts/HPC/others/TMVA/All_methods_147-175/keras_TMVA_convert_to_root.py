import ROOT
from ROOT import TMVA, TFile, TTree, TCut, TString
import matplotlib.pyplot as plt
import numpy as np
from root_numpy import array2tree, array2root


negative = np.load('negative_147_175_pattern_0.npy')
negative2 = negative[:,0:16].copy()
negative2[:,4] = np.fabs(negative2[:,4])
negative2[:,5] = np.fabs(negative2[:,5])
negative2[:,6] = np.fabs(negative2[:,6])
negative2.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(negative2, 'negative.root', 'tree', 'recreate')

positive = np.load('positive_147_175_pattern_0.npy')
positive2 = positive[:,0:16].copy()
positive2[:,4] = np.fabs(positive2[:,4])
positive2[:,5] = np.fabs(positive2[:,5])
positive2[:,6] = np.fabs(positive2[:,6])
positive2.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(positive2, 'positive.root', 'tree', 'recreate')


########################################################################################################################

negative_ele = np.load('negative_147_175_pattern_0_ele.npy')
negative2_ele = negative_ele[:,0:16].copy()
negative2_ele[:,4] = np.fabs(negative2_ele[:,4])
negative2_ele[:,5] = np.fabs(negative2_ele[:,5])
negative2_ele[:,6] = np.fabs(negative2_ele[:,6])
negative2_ele.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(negative2_ele, 'negative_ele.root', 'tree', 'recreate')

positive_ele = np.load('positive_147_175_pattern_0_ele.npy')
positive2_ele = positive_ele[:,0:16].copy()
positive2_ele[:,4] = np.fabs(positive2_ele[:,4])
positive2_ele[:,5] = np.fabs(positive2_ele[:,5])
positive2_ele[:,6] = np.fabs(positive2_ele[:,6])
positive2_ele.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(positive2_ele, 'positive_ele.root', 'tree', 'recreate')



#########################################################################################################################

negative_test = np.load('negative_147_175_pattern_0_test.npy')
negative2_test = negative_test[:,0:16].copy()
positive_test = np.load('positive_147_175_pattern_0_test.npy')
positive2_test = positive_test[:,0:16].copy()

negative2_test[:,4] = np.fabs(negative2_test[:,4])
negative2_test[:,5] = np.fabs(negative2_test[:,5])
negative2_test[:,6] = np.fabs(negative2_test[:,6])
positive2_test[:,4] = np.fabs(positive2_test[:,4])
positive2_test[:,5] = np.fabs(positive2_test[:,5])
positive2_test[:,6] = np.fabs(positive2_test[:,6])

negative2_test.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(negative2_test, 'negative_test.root', 'tree', 'recreate')
positive2_test.dtype = [('a1','float32'), ('a2','float32'), ('a3','float32'), ('a4','float32'), ('a5','float32'),('a6','float32'),('a7','float32'),('a8','float32'),('a9','float32'),('a10','float32'),('a11','float32'),('a12','float32'),('a13','float32'),('a14','float32'),('a15','float32'),('a16','float32')]
array2root(positive2_test, 'positive_test.root', 'tree', 'recreate')






