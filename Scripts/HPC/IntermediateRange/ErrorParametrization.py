import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F
import uproot
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import binning
import argparse
import math
from iminuit import Minuit
from iminuit.util import describe, make_func_code
from iminuit.cost import LeastSquares
import uncertainties
import ParameterizationTemplates_function
import inspect, re
import os
import ErrorParametrization_CovNumber

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def Plot(x, y, covnumber, parametername, parameterindex):
    plt.figure(figsize=(18,9))

    plt.plot(x, y, "*", markersize=10)
    if covnumber>0:
        plt.plot(x, smooth(y, covnumber),  "o", markersize=10)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel(parametername + str(parameterindex) + ' error',fontsize=30)

    plt.savefig('ErrorParametrization/' + str(parametername) + str(parameterindex) + str(".pdf"))
    plt.close()



def main():

    #filenamelist = ['ParameterError_TOFAntiproton', 'ParameterError_TOFElectron', 'ParameterError_TOFPion', 'ParameterError_TRDAntiproton', 'ParameterError_TRDElectron', 'ParameterError_TRDPion']
    filenamelist = ['ParameterError_TOFAntiproton', 'ParameterError_TOFElectron', 'ParameterError_TOFPion', 'ParameterError_TRDAntiproton', 'ParameterError_TRDElectron', 'ParameterError_TRDPion']

    for filename in filenamelist:
        npyfile = np.load(filename+'.npy')
        
        ParameterError_new = []

        ## Plot
        for i in range(npyfile.shape[1]):
            print("filename:"+ str(filename) + ", " + "parameterindex" + str(i))
            CovNumber = ErrorParametrization_CovNumber.CovNumber(filename, i)
            print("CovNumber:"+ str(CovNumber))
            Plot(BinningCenterAll, npyfile[:,i], CovNumber, filename, i)
            if CovNumber>0:
                ParameterError_new.append(smooth(npyfile[:,i], CovNumber)) 
            else:
                ParameterError_new.append(npyfile[:,i])
           
        np.save(filename+'_new.npy', np.array(ParameterError_new).T)    


if __name__ == '__main__':

    #### Parser Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--Range'         , help='Rigidity Range, low, intermediate or high')
    parser.add_argument('--issversion'    , help='ISS data version, pass7.8 or 2016paper')
    arguments = parser.parse_args()

    if (arguments.issversion and arguments.issversion):
        ISSversion = arguments.issversion
        RigidityRange  = arguments.Range
    else:
        print("You need to provide all parameters!")
        os._exit(0)

    IntermediatePath = os.getenv('HPCINTERMEDIATEDIR')
    LowPath          = os.getenv('HPCLOWENERGYDIR')

    if RigidityRange == 'intermediate':
        BinningAll       = binning.binslow[10:]
        BinningAll[19]   = '16.6_18'
        BinningCenterAll = binning.published2016binnings_center[9:29]
    elif RigidityRange == 'low':
        BinningAll       = binning.binslow[1:17]
        BinningAll[0]    = '1_1.16'
        BinningCenterAll = binning.published2016binnings_center[0:16]

    main()









