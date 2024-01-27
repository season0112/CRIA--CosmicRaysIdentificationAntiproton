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

def Plot_Trd_OneGaussOneNovo(x_Center_trd, x_More_trd, y_Value_trd, y_trdfit_value, y_trdfit_value_more, y_trdfit_value_uncer_more, y_trdfit_error_uncer_more, labelname, binname, reduced_chi_squared, RigidityRange, Para_TRD, Para_TRD_uncertainty):

    GaussPart_more             = ParameterizationTemplates_function.Gauss_function      (x_More_trd, Para_TRD[0], Para_TRD[2]    , Para_TRD[3])
    NovosibirskPart_more       = ParameterizationTemplates_function.Novosibirsk_function(x_More_trd, Para_TRD[1], Para_TRD[4]    , Para_TRD[5]    , Para_TRD[6])
    #GaussPart_uncer_more       = ParameterizationTemplates_function.Gauss_function      (x_More_trd, Para_TRD_uncertainty[0], Para_TRD_uncertainty[2], Para_TRD_uncertainty[3])
    #NovosibirskPart_uncer_more = ParameterizationTemplates_function.Novosibirsk_function(x_More_trd, Para_TRD_uncertainty[1], Para_TRD_uncertainty[4], Para_TRD_uncertainty[5], Para_TRD_uncertainty[6])    
    #print("GaussPart_more" + str(GaussPart_more))
    #print("NovosibirskPart_more" + str(NovosibirskPart_more))
    #print("y_trdfit_value_more" + str(y_trdfit_value_more))

    #### Plot
    fig, ax = plt.subplots(figsize=(18, 9))

    plt.plot(x_More_trd, y_trdfit_value_more      , color = 'blue',  label = 'Fit')
    #plt.plot(x_More_trd, GaussPart_more          , color = 'red',   linestyle='dashed', label = 'Gauss')
    #plt.plot(x_More_trd, NovosibirskPart_more    , color = 'green', linestyle='dashed', label = 'Novosibirsk')

    #plt.plot(x_More_trd, y_trdfit_value_more     ,  color = 'm'   , label = labelname + str("(new)"))
    #plt.plot(x_More_trd, GaussPart_new_more      ,  color = 'lime', linestyle='dashed', label = 'Gauss (new)')
    #plt.plot(x_More_trd, NovosibirskPart_new_more,  color = 'cyan', linestyle='dashed', label = 'Novosibirsk (new)')

    plt.errorbar(x_Center_trd, y_Value_trd, xerr=0, yerr=np.sqrt(y_Value_trd), label = 'Data', fmt='o', markersize=10, color='black')

    #plt.text(0.2, 0.2, "chie/dof:" + "{:.2f}".format(reduced_chi_squared), fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('TRDLikelihood',fontsize=30)
    plt.legend( loc='best', fontsize=20 )

    plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_linear.pdf" )

    #plt.yscale('log')
    #plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_log.pdf" )
    plt.close()


def Plot_Trd_OneGauss(x, x_more, y, y_fit, y_fit_more, y_fit_value_new_more, y_fit_error_new_more, labelname, binname, reduced_chi_squared, RigidityRange):

    ####
    fig, ax = plt.subplots(figsize=(18, 9))

    plt.plot(x_more, y_fit_more          , color = 'blue', label = 'Fit')
    #plt.plot(x_more, y_fit_value_new_more, color = 'm'   , label = labelname + str("(new)"))

    plt.errorbar(x, y, xerr=0, yerr=np.sqrt(y), label = 'Data', fmt='o', markersize=10, color='black')

    #plt.text(0.2, 0.2, "chie/dof:" + "{:.2f}".format(reduced_chi_squared), fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('TRDLikelihood',fontsize=30)
    plt.legend( loc='best', fontsize=20 )

    plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_linear.pdf" )

    #plt.yscale('log')
    #plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_log.pdf" )
    plt.close()


def Plot_Trd_OneNovo(x, x_more, y, y_fit, y_fit_more, y_fit_value_new_more, y_fit_error_new_more, labelname, binname, reduced_chi_squared, RigidityRange):
    fig, ax = plt.subplots(figsize=(18, 9))

    plt.plot(x_more, y_fit_more,            color = 'blue',  label = 'Fit')
    #plt.plot(x_more, y_fit_value_new_more    ,  color = 'm'   , label = labelname + str("(new)"))

    plt.errorbar(x, y, xerr=0, yerr=np.sqrt(y), label = 'Data', fmt='o', markersize=10, color='black')

    #plt.text(0.2, 0.2, "chie/dof:" + "{:.2f}".format(reduced_chi_squared), fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('TRDLikelihood',fontsize=30)
    plt.legend( loc='best', fontsize=20 )

    plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_linear.pdf" )

    #plt.yscale('log')
    #plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_log.pdf" )
    plt.close()


def Plot_Trd_ExponentialFun(x       , x_more, y       , y_fit      , y_fit_more      , y_fit_value_new_more, y_fit_error_new_more, labelname, binname, reduced_chi_squared, RigidityRange):

    fig, ax = plt.subplots(figsize=(18, 9))

    plt.plot(x_more, y_fit_more          ,  color = 'blue',  label = 'Fit')
    #plt.plot(x_more, y_fit_value_new_more,  color = 'm'   , label = labelname + str("(new)"))

    plt.errorbar(x, y, xerr=0, yerr=np.sqrt(y), label = 'Data', fmt='o', markersize=10, color='black')

    #plt.text(0.2, 0.2, "chie/dof:" + "{:.2f}".format(reduced_chi_squared), fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('TRDLikelihood',fontsize=30)
    plt.legend( loc='best', fontsize=20 )

    plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_linear.pdf" )

    #plt.yscale('log')
    #plt.savefig( "Plot/ParametrizationPlot/TRD_" + labelname + "_" + binname + RigidityRange + "_log.pdf" )
    plt.close()


def Plot_Tof_OneGauss(x_Center_tof, x_More_tof, y_Value_tof, y_toffit_value, y_toffit_value_more, y_toffit_value_uncer_more, y_toffit_error_uncer_more, labelname, binname, reduced_chi_squared_TOF, RigidityRange):

    fig, ax = plt.subplots(figsize=(18, 9))

    plt.plot(x_More_tof, y_toffit_value_more      , color = 'blue', label = 'Fit')
    #plt.plot(x_More_tof, y_toffit_value_uncer_more, color = 'm'   , label = 'Fit' + str("(new)"))

    plt.errorbar(x_Center_tof, y_Value_tof, xerr=0, yerr=np.sqrt(y_Value_tof), label = 'Data', fmt='o', markersize=10, color='black')

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    #plt.xlabel(r'$\frac{1}{TOFBeta}$ - $\frac{1}{Beta(R,m_p)}$',fontsize=30)
    plt.xlabel( '1/TOFBeta - 1/Beta(R,m_p)' , fontsize=30)
    plt.legend( loc='best', fontsize=20 )

    plt.savefig( "Plot/ParametrizationPlot/TOF_" + labelname + "_" + binname + RigidityRange + "_linear.pdf" )

    #plt.yscale('log')
    #plt.savefig( "Plot/ParametrizationPlot/TOF_" + labelname + "_" + binname + RigidityRange + "_log.pdf" )
    plt.close()


def Plot_Tof_OneNovo(x_Center_tof, x_More_tof, y_Value_tof, y_toffit_value, y_toffit_value_more, y_toffit_value_uncer_more, y_toffit_error_uncer_more, labelname, binname, reduced_chi_squared_TOF, RigidityRange):

    fig, ax = plt.subplots(figsize=(18, 9))

    plt.plot(x_More_tof, y_toffit_value_more      , color = 'blue', label = 'Fit')
    #plt.plot(x_More_tof, y_toffit_value_uncer_more, color = 'm'   , label = labelname + str("(new)"))

    plt.errorbar(x_Center_tof, y_Value_tof, xerr=0, yerr=np.sqrt(y_Value_tof), label = 'Data', fmt='o', markersize=10, color='black')

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel( '1/TOFBeta - 1/Beta(R,m_p)' , fontsize=30)
    plt.legend( loc='best', fontsize=20 )

    plt.savefig( "Plot/ParametrizationPlot/TOF_" + labelname + "_" + binname + RigidityRange + "_linear.pdf" )

    #plt.yscale('log')
    #plt.savefig( "Plot/ParametrizationPlot/TOF_" + labelname + "_" + binname + RigidityRange + "_log.pdf" )
    plt.close()


def Plot_GeneratedParameters( parameterset, name, spicename, binname, RigidityRange):

    for i in range(len(parameterset)):
        fig, ax = plt.subplots(figsize=(18, 9))

        plt.hist(parameterset[i,:], 30, density=False, facecolor='g', alpha=0.75)

        plt.xlabel('p'+str(i+1), fontsize=30)
        plt.text(0.9, 0.9, binname + str("GV"), fontsize=20, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        plt.yscale('log')
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)

        plt.savefig('Plot/ParametersetPlot/' + name + '_p' + str(i+1) + "_" + binname + "_" + spicename + '_' + RigidityRange + ".pdf")


def Plot_ParavsRigidty(Rigidity, ParameterSet, spice, dimension, RigidityRange):
    ParameterSet = np.array(ParameterSet)

    for i in range(ParameterSet.shape[1]):

        plt.figure(figsize=(18,9))

        plt.plot(Rigidity, ParameterSet[:,i], "*",  markersize=10)

        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.xlabel('Rigidity (GV)', fontsize=30)
        plt.ylabel('Parameter' + str(i), fontsize=30)

        plt.savefig( 'Plot/ParavsRigidtyPlot/Parameter_' + str(i) + "_" + dimension + "_" + str(spice) + RigidityRange + str(".pdf") )
        plt.close()
    

def Plot_Chi2(x, y, name, RigidityRange):

    plt.figure(figsize=(18,9))

    plt.plot(x, y, "*", markersize=10)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Chi2/dof',fontsize=30)

    #plt.ylim(-1, 500)

    #plt.savefig('Plot/chi2Plot/' +name + RigidityRange + str(".pdf"))
    plt.savefig(name + RigidityRange + str(".pdf"))
    plt.close()


def Plot_GeneratedMultiTemplates(x_Center, x_More, y_Value, parameterset, value, binname, RigidityRange, UsedFunction, GenerateNumber, Dimensionname):

    #### Plot
    fig, ax = plt.subplots(figsize=(18, 9))

    plt.errorbar(x_Center, y_Value, xerr=0, yerr=np.sqrt(y_Value), label = 'Data', fmt='o', markersize=10, color='red')

    for i in range(GenerateNumber):
        plt.plot(x_More, UsedFunction(x_More, *parameterset[:,i]), color = 'blue', linestyle='dotted') 


    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel(Dimensionname, fontsize=30)
    plt.legend( loc='best', fontsize=20 )

    plt.savefig( "Plot/MultiTemplatePlot/MultiTemplate_" + Dimensionname + "_" + value + "_" + binname + RigidityRange + "_linear.pdf" )

    plt.yscale('log')
    plt.savefig( "Plot/MultiTemplatePlot/MultiTemplate_" + Dimensionname + "_" + value + "_" + binname + RigidityRange + "_log.pdf" )
    plt.close()
