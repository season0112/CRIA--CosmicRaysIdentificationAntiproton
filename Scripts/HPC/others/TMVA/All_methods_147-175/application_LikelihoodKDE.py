import ROOT
from ROOT import TMVA, TFile, TTree, TCut, TString
import matplotlib.pyplot as plt
import numpy as np
import array 
from root_numpy import array2tree, array2root, tree2array

reader = ROOT.TMVA.Reader()

v1 = array.array('f',[0])
v2 = array.array('f',[0])
v3 = array.array('f',[0])
v4 = array.array('f',[0])
v5 = array.array('f',[0])
v6 = array.array('f',[0])
v7 = array.array('f',[0])
v8 = array.array('f',[0])
v9 = array.array('f',[0])
v10 = array.array('f',[0])
v11 = array.array('f',[0])
v12 = array.array('f',[0])
v13 = array.array('f',[0])
v14 = array.array('f',[0])
v15 = array.array('f',[0])
v16 = array.array('f',[0])

reader.AddVariable("a1",v1)
reader.AddVariable("a2",v2)
reader.AddVariable("a3",v3)
reader.AddVariable("a4",v4)
reader.AddVariable("a5",v5)
reader.AddVariable("a6",v6)
reader.AddVariable("a7",v7)
reader.AddVariable("a8",v8)
reader.AddVariable("a9",v9)
reader.AddVariable("a10",v10)
reader.AddVariable("a11",v11)
reader.AddVariable("a12",v12)
reader.AddVariable("a13",v13)
reader.AddVariable("a14",v14)
reader.AddVariable("a15",v15)
reader.AddVariable("a16",v16)


#reader.BookMVA("BDT","dataset_pymva/weights/TMVAClassification_BDT.weights.xml")
reader.BookMVA("LikelihoodKDE","dataset_pymva/weights/TMVAClassification_LikelihoodKDE.weights.xml")
########################################################################################
positivefile_test = TFile.Open("positive_test.root")
negativefile_test = TFile.Open("negative_test.root")
negativefile = TFile.Open("negative.root")
# Get signal and background trees from file
signal = positivefile_test.Get('tree')
background1 = negativefile_test.Get('tree')
background2 = negativefile.Get('tree')

signal.SetEntries(background1.GetEntries()+background2.GetEntries())

signal_array = tree2array(signal)
background1_array = tree2array(background1)
background2_array = tree2array(background2)
background_array = np.append(background1_array, background2_array)
print(background1_array.shape)
print(background2_array.shape)
print(background_array.shape)

mva_LikelihoodKDE_signal_array = np.array([])
mva_LikelihoodKDE_background_array = np.array([])

for ievt in range(signal_array.shape[0]):
	v1[0] = signal_array['a1'][ievt]
	v2[0] = signal_array['a2'][ievt]
	v3[0] = signal_array['a3'][ievt]
	v4[0] = signal_array['a4'][ievt]
	v5[0] = signal_array['a5'][ievt]
	v6[0] = signal_array['a6'][ievt]
	v7[0] = signal_array['a7'][ievt]
	v8[0] = signal_array['a8'][ievt]
	v9[0] = signal_array['a9'][ievt]
	v10[0] = signal_array['a10'][ievt]
	v11[0] = signal_array['a11'][ievt]
	v12[0] = signal_array['a12'][ievt]
	v13[0] = signal_array['a13'][ievt]
	v14[0] = signal_array['a14'][ievt]
	v15[0] = signal_array['a15'][ievt]
	v16[0] = signal_array['a16'][ievt]
	mva_LikelihoodKDE = reader.EvaluateMVA("LikelihoodKDE")
	mva_LikelihoodKDE_signal_array = np.append(mva_LikelihoodKDE_signal_array,mva_LikelihoodKDE)      
	print("ievt="+str(ievt))
	print("mva_LikelihoodKDE="+str(mva_LikelihoodKDE))
np.save("mva_LikelihoodKDE_signal_array.npy", mva_LikelihoodKDE_signal_array)

for ievt in range(background_array.shape[0]):
        v1[0] = background_array['a1'][ievt]
        v2[0] = background_array['a2'][ievt]
        v3[0] = background_array['a3'][ievt]
        v4[0] = background_array['a4'][ievt]
        v5[0] = background_array['a5'][ievt]
        v6[0] = background_array['a6'][ievt]
        v7[0] = background_array['a7'][ievt]
        v8[0] = background_array['a8'][ievt]
        v9[0] = background_array['a9'][ievt]
        v10[0] = background_array['a10'][ievt]
        v11[0] = background_array['a11'][ievt]
        v12[0] = background_array['a12'][ievt]
        v13[0] = background_array['a13'][ievt]
        v14[0] = background_array['a14'][ievt]
        v15[0] = background_array['a15'][ievt]
        v16[0] = background_array['a16'][ievt]
        mva_LikelihoodKDE = reader.EvaluateMVA("LikelihoodKDE")
        mva_LikelihoodKDE_background_array = np.append(mva_LikelihoodKDE_background_array,mva_LikelihoodKDE)
        print("ievt="+str(ievt))
        print("mva_LikelihoodKDE="+str(mva_LikelihoodKDE))
np.save("mva_LikelihoodKDE_background_array.npy", mva_LikelihoodKDE_background_array)



