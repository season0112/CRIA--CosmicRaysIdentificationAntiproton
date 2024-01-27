import numpy as np
import binning
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import PythonPlotDefaultParameters

def Load(lowworkpath):
    TRD = np.loadtxt(lowworkpath + '/totalall/TRDLogLikelihood_CutValue_eff_0.94_1.txt')
    TOF = np.loadtxt(lowworkpath + '/totalall/TOFBeta_CutValue_eff_0.95_1.txt')
    return TRD, TOF


def Plot_TRD(lowworkpath, TRD):
    plt.figure(figsize=(35,18))
    ax=plt.gca()

    plt.plot(binning.published2016binnings_center[0:TRD.shape[0]]     , TRD, color='black', linewidth=10.0)
    plt.errorbar( binning.published2016binnings_center[0:TRD.shape[0]], TRD, yerr=0, fmt='o', markersize=35, color='r', label='')

    plt.xlabel("Rigidity (GV)"             , fontsize=70, horizontalalignment='right', x=1.0)
    plt.ylabel( r'$\Lambda_{\rm{LowEdge}}$', fontsize=80, horizontalalignment='right', y=1.0)
    plt.xticks(fontsize=70)
    plt.yticks(fontsize=60)
    ax.tick_params(axis='both', which='both', direction='in', length=30, width=10)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))

    plt.savefig(lowworkpath + "/totalall/LowerEdge_TRD.pdf")


def Plot_TOF(lowworkpath, TOF):
    plt.figure(figsize=(35,18))
    ax=plt.gca()

    plt.plot( binning.published2016binnings_center[0:TOF.shape[0]]    , TOF, color='black', linewidth=10.0)
    plt.errorbar( binning.published2016binnings_center[0:TOF.shape[0]], TOF, yerr=0, fmt='o', markersize=35, color='r', label='')

    plt.xlabel("Rigidity (GV)"             , fontsize=70, horizontalalignment='right', x=1.0)
    plt.ylabel( r'1/$\beta_{\rm{LowEdge}}$', fontsize=80, horizontalalignment='right', y=1.0, labelpad=15)

    plt.xticks(fontsize=70)
    plt.yticks(fontsize=60)

    ax.tick_params(axis='both', which='both', direction='in', length=30, width=10)

    plt.savefig(lowworkpath + "/totalall/LowerEdge_TOF.pdf")


def main():
    lowworkpath = os.getenv('HPCLOWENERGYDIR')

    TRD, TOF = Load(lowworkpath)

    Plot_TRD(lowworkpath, TRD)
    Plot_TOF(lowworkpath, TOF)


if __name__ == '__main__':

    main()






