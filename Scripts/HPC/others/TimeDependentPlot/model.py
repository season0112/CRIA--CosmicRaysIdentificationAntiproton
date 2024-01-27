import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import argparse
import binning
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import scipy.integrate as integrate
from iminuit import Minuit
from iminuit.cost import LeastSquares
import sys
import uncertainties
import PythonPlotDefaultParameters
from decimal import *
import decimal
import model
from uncertainties import ufloat
from uncertainties import unumpy
from iminuit.util import describe, make_func_code
import matplotlib.ticker as mticker



#### Force-Field Approximation model
# definition of simple power-law
def powerlaw(x, xref, C, gamma):
    #return C * np.power(x/xref, gamma)
    return np.array(C * (x/xref) ** gamma, dtype=np.float64)
# convert energy to rigidity for given mass (Z=1 only!)
def rigidity(E, m):
    return np.sqrt(E**2 - m**2)
# convert rigidity to energy for given mass
def energy(R, m):
    return np.sqrt(R**2 + m**2)
# power-law in rigidity with force-field approximation:
# J_TOA(E) = (E^2-m^2) / ((E+phi)^2-m^2) * J_LIS(E+phi)
# where E is total energy (!)
#
# transformation of variables: dN/dR = dN/dE * dE/dR
# here we have dE/dR = R/E
def powerlaw_solarmod(R, phi, C, gamma):
    #m = 0.938272 # change mass here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m = 0.000511
    E = energy(R, m)
    Elis = E + phi
    Rlis = rigidity(Elis, m)
    return R**2/Rlis**2 * powerlaw(Rlis, 1.0, C, gamma) * R/E * Elis/Rlis
# LaffertyWyatt
def LaffertyWyatt(E1,E2,gamma):
    denom = E2**(1-gamma)-E1**(1.0-gamma)
    num = (E2-E1)*(1.0-gamma)
    return (num / denom)**(1.0/gamma)


# exact fit, where model is integrated numerically over each bin
# before (integral) fluxes are compared
# unfortunately, this is CPU-intensive (when coded in python)
class ExactLeastSquares:
    def __init__(self, electron_R, Rmin, Rmax, flux, err_total, R1_fit, R2_fit):
        self.electron_R = electron_R
        self.Rmin = Rmin
        self.Rmax = Rmax
        self.flux = flux
        self.err_total = err_total
        self.R1_fit = R1_fit
        self.R2_fit = R2_fit
    def __call__(self, phi, C, gamma):
        # calculate chi2 over fit range
        indexlist = np.where( (np.array(self.electron_R)>self.R1_fit) & (np.array(self.electron_R)<self.R2_fit) )[0]
        Rmin_new=[]
        Rmax_new=[]
        for i in indexlist:
            Rmin_new.append(self.Rmin[np.where(self.Rmin < self.electron_R[i])[0][-1]])
            Rmax_new.append(self.Rmax[np.where(self.Rmax > self.electron_R[i])[0][0]])
        integrand = lambda R: powerlaw_solarmod(R, phi, C, gamma)
        #model_binned = np.array( [integrate.quad(integrand, R1, R2)[0] / (R2-R1) for R1, R2 in zip(self.Rmin[indexlist], self.Rmax[indexlist])] )
        model_binned = np.array( [integrate.quad(integrand, R1, R2)[0] / (R2-R1) for R1, R2 in zip(Rmin_new, Rmax_new)] )
        #electron_R = [electron_R[i] for i in indexlist]
        pulls = (np.array(self.flux)[indexlist] - np.array(model_binned)) / np.array(self.err_total)[indexlist]
        chi2 = np.sum(pulls**2)
        return chi2


####
def binnin_revised(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, index, binmerge, group):

    Emin = binning.published2016binnings[index]
    Emax = binning.published2016binnings[index+int(binmerge)]

    #Emin = 1.72
    #Emax = 2.00
    #Emin = 4.54
    #Emax = 5.00

    #### Parameters
    # start values for first fits (will be updated in each time bin)
    phi_ele=1.64
    C_ele=3104
    gamma_ele=-3.7

    phi_pos=0.87
    C_pos=46.8
    gamma_pos=-3.36

    # global fit range
    if rigidityrange == "intermediate":
        #R1_fit = binning.published2016binnings[index-9]
        #R2_fit = binning.published2016binnings[index+10]
        R1_fit = binning.published2016binnings[index-9]
        R2_fit = binning.published2016binnings[index+17]
    elif rigidityrange == "low":
        #R1_fit = binning.published2016binnings[index-3]
        #R2_fit = binning.published2016binnings[index+20]
        R1_fit = binning.published2016binnings[index-1]
        R2_fit = binning.published2016binnings[index+21]

    print("Emin and Emax: " + str(Emin) + "_to_" + str(Emax))
    print("fit range: " + str(R1_fit) + "_to_" + str(R2_fit))


    # loop over proton data and store results
    #(1) rough method
    time_array = binning.Bartals1Unixtime_edge
    electron_flux_t_interpolated = []
    electron_err_t_interpolated = []
    positron_flux_t_interpolated = []
    positron_err_t_interpolated = []
    '''
    #(2) exact method
    electron_flux_t_interpolated_exact = []
    positron_flux_t_interpolated_exact = []
    electron_err_t_interpolated_exact = []
    positron_err_t_interpolated_exact = []
    '''

    if group == "Aachen":
        timerange = [x for x in range(0,88) if (x != 46 and x != 47)] # 46,47 are empty due to TTCSoff.
    elif group == "Perugia":
        timerange = timerange=range(100)  # 46,47 are not included in the result, therefore no need to pick them out.

    for timeindex in timerange: # 46,47 are empty due to TTCSoff.
        electron_R_all=[]
        electron_flux = []
        electron_fluxerror = []
        positron_R=[]
        positron_flux = []
        positron_fluxerror = []

        if group == "Aachen":
            for i in range(lepton_result.Get("SxFluxTools_Electron_1305417600_1307750400").Get("grElectron_1305417600_1307750400_FinalFlux").GetN()): #grElectron_1324080000_1326412800_FinalFlux.GetN()=52
                electron_R_all.append(lepton_result.Get("SxFluxTools_Electron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))).Get("grElectron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))+"_FinalFlux").GetX()[i])
                electron_flux.append(lepton_result.Get("SxFluxTools_Electron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))).Get("grElectron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))+"_FinalFlux").GetY()[i])
                electron_fluxerror.append(lepton_result.Get("SxFluxTools_Electron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))).Get("grElectron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))+"_FinalFlux").GetErrorY(i)) # double TGraphAsymmErrors::GetErrorY(int bin)

                positron_R.append(lepton_result.Get("SxFluxTools_Positron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))).Get("grPositron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))+"_FinalFlux").GetX()[i])
                positron_flux.append(lepton_result.Get("SxFluxTools_Positron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))).Get("grPositron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))+"_FinalFlux").GetY()[i])
                positron_fluxerror.append(lepton_result.Get("SxFluxTools_Positron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))).Get("grPositron_"+str(int(time_array[timeindex]))+"_"+str(int(time_array[timeindex+1]))+"_FinalFlux").GetErrorY(i))

        elif group == "Perugia":
            for i in range(3,46): #In total:45 points, 1.01-40.00GeV, but (1).last points 37.31-40.0 not 37.31-39.39,so remove 1 point; (2). binning.LeptonLowEnergyBin[i+1] in loop, so stop at 46. (binning.LeptonLowEnergyBin[3]=1.01, binning.LeptonLowEnergyBin[46]=35.36, binning.LeptonLowEnergyBin[47]=37.31) 
                g_electronflux = lepton_result_Maura.Get("ele_maura_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[i])) + "_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[i+1])) + ";1")
                g_positronflux = lepton_result_Maura.Get("pos_maura_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[i])) + "_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[i+1])) + ";1")

                #electron_R_all.append(g_electronflux.GetX()[timeindex])
                electron_R_all.append( (binning.LeptonLowEnergyBin[i] + binning.LeptonLowEnergyBin[i+1])/2 )
                electron_flux.append(g_electronflux.GetY()[timeindex])
                electron_fluxerror.append(g_electronflux.GetErrorY(timeindex))

                #positron_R.append(g_positronflux.GetX()[timeindex])
                positron_R.append( (binning.LeptonLowEnergyBin[i] + binning.LeptonLowEnergyBin[i+1])/2 )
                positron_flux.append(g_positronflux.GetY()[timeindex])
                positron_fluxerror.append(g_positronflux.GetErrorY(timeindex))

        '''
        # Exact method
        binningedge = binning.LeptonLowEnergyBin
        Rmin = binningedge[0:-1]
        Rmax = binningedge[1:]
        if timeindex<90:
            print(f'{timeindex}  (Exact method)')

            lsq_exact_ele = ExactLeastSquares(electron_R_all, Rmin, Rmax, electron_flux, electron_fluxerror, R1_fit, R2_fit)
            m_exact_ele = Minuit(lsq_exact_ele, phi=phi_ele, C=C_ele, gamma=gamma_ele)    
            m_exact_ele.errordef = 1.0
            m_exact_ele.migrad()
            m_exact_ele.hesse()

            lsq_exact_pos = ExactLeastSquares(positron_R, Rmin, Rmax, positron_flux, positron_fluxerror, R1_fit, R2_fit)
            m_exact_pos = Minuit(lsq_exact_pos, phi=phi_pos, C=C_pos, gamma=gamma_pos)
            m_exact_pos.errordef = 1.0
            m_exact_pos.migrad()
            m_exact_pos.hesse()
        '''

        # rough method: do a simple chi2 fit
        # fit range
        indexlist = np.where( (np.array(electron_R_all)>R1_fit) & (np.array(electron_R_all)<R2_fit) )[0]
        electron_R = [electron_R_all[i] for i in indexlist]
        electron_flux = [electron_flux[i] for i in indexlist]
        electron_fluxerror = [electron_fluxerror[i] for i in indexlist]        
        positron_R = [positron_R[i] for i in indexlist] # assume electron_R and positron_R are same.
        positron_flux = [positron_flux[i] for i in indexlist]
        positron_fluxerror = [positron_fluxerror[i] for i in indexlist]
        # electron fit
        leastsquare_rough_ele = LeastSquares(electron_R, electron_flux, electron_fluxerror, model.powerlaw_solarmod)
        m_rough_ele = Minuit(leastsquare_rough_ele, phi=phi_ele, C=C_ele, gamma=gamma_ele) # Minuit(fcn: Callable, *args, grad: Optional[Callable] = None, name: Optional[Sequence[str]] = None, **kwds)
        m_rough_ele.errordef = 1.0 # Access FCN increment above the minimum that corresponds to one standard deviation. Default value is 1.0. errordef should be 1.0 for a least-squares cost function and 0.5 for a negative log-likelihood function.
        m_rough_ele.migrad() # Run Migrad minimization.
        m_rough_ele.hesse() # Run HESSE algorithm to compute asymptotic errors.
        # positron fit
        leastsquare_rough_pos = LeastSquares(positron_R, positron_flux, positron_fluxerror, model.powerlaw_solarmod)
        m_rough_pos = Minuit(leastsquare_rough_pos, phi=phi_pos, C=C_pos, gamma=gamma_pos)
        m_rough_pos.errordef = 1.0
        m_rough_pos.migrad()
        m_rough_pos.hesse()

        # update start values for next time bin
        #(1). rough method
        phi_ele, C_ele, gamma_ele = m_rough_ele.values
        phi_pos, C_pos, gamma_pos = m_rough_pos.values

        # model curves with best-fit parameters
        xR = np.logspace(np.log10(R1_fit), np.log10(R2_fit), 200)
        #(1). rough method
        pl_approx_rough_ele = lambda R: model.powerlaw_solarmod(R, *m_rough_ele.values)
        pl_approx_rough_pos = lambda R: model.powerlaw_solarmod(R, *m_rough_pos.values)
        #(2). exact method
        #pl_approx_exact_ele = lambda R: model.powerlaw_solarmod(R, *m_exact_ele.values)
        #pl_approx_exact_pos = lambda R: model.powerlaw_solarmod(R, *m_exact_pos.values)

        # integrate local parameterization over rigidity range 
        #(1). rough method
        interpolated_flux_ele = integrate.quad(pl_approx_rough_ele, Emin, Emax)[0] / (Emax - Emin)
        interpolated_flux_pos = integrate.quad(pl_approx_rough_pos, Emin, Emax)[0] / (Emax - Emin)
        electron_flux_t_interpolated.append(interpolated_flux_ele)
        positron_flux_t_interpolated.append(interpolated_flux_pos)
        #(2). exact method
        #interpolated_flux_ele_exact = integrate.quad(pl_approx_exact_ele, Emin, Emax)[0] / (Emax - Emin)
        #interpolated_flux_pos_exact = integrate.quad(pl_approx_exact_pos, Emin, Emax)[0] / (Emax - Emin)
        #electron_flux_t_interpolated_exact.append(interpolated_flux_ele_exact)
        #positron_flux_t_interpolated_exact.append(interpolated_flux_pos_exact)
        
        # positron/electron data is stored as a function of energy E
        # proton data is stored as a function of rigidity R
        # for positrons and electrons, we have E=R
        E = model.LaffertyWyatt(Emin,Emax,gamma=3.0)
        R = E

        
        # Error propagation1 (same as exercise)
        # error propagation for integrated flux, calculate relative error at fiducial point
        # and assume that the relative error can be used for the integrated flux, too
        # (1). rough method
        params_ele = uncertainties.correlated_values(m_rough_ele.values, np.array(m_rough_ele.covariance))
        wrap_func_ele = uncertainties.wrap(model.powerlaw_solarmod)
        val_ele = wrap_func_ele(R, *params_ele)
        rel_error_ele = val_ele.s / val_ele.n #".s" is ".std_dev"; ".n" is ".nominal_value"
        electron_err_t_interpolated.append(interpolated_flux_ele * rel_error_ele)

        params_pos = uncertainties.correlated_values(m_rough_pos.values, np.array(m_rough_pos.covariance))
        wrap_func_pos = uncertainties.wrap(model.powerlaw_solarmod)
        val_pos = wrap_func_pos(R, *params_pos)
        rel_error_pos = val_pos.s / val_pos.n 
        positron_err_t_interpolated.append(interpolated_flux_pos * rel_error_pos)
        # (2). exact method
        #params_ele_exact = uncertainties.correlated_values(m_exact_ele.values, np.array(m_exact_ele.covariance))
        #wrap_func_ele_exact = uncertainties.wrap(model.powerlaw_solarmod)
        #val_ele_exact = wrap_func_ele_exact(R, *params_ele_exact)
        #rel_error_ele_exact = val_ele_exact.s / val_ele_exact.n 
        #electron_err_t_interpolated_exact.append(interpolated_flux_ele_exact * rel_error_ele_exact)

        #params_pos_exact = uncertainties.correlated_values(m_exact_pos.values, np.array(m_exact_pos.covariance))
        #wrap_func_pos_exact = uncertainties.wrap(model.powerlaw_solarmod)
        #val_pos_exact = wrap_func_pos_exact(R, *params_pos_exact)
        #rel_error_pos_exact = val_pos_exact.s / val_pos_exact.n 
        #positron_err_t_interpolated_exact.append(interpolated_flux_pos_exact * rel_error_pos_exact)
        

        ''' 
        # Error propagation2 (find closest bin and assume relative error is same)
        Error_Index_ele = electron_R.index(min(electron_R, key=lambda x: abs(x - R))) 
        #print("closest bin is" + str(electron_R[Error_Index_ele]))
        rel_error_ele = electron_fluxerror[Error_Index_ele] / electron_flux[Error_Index_ele]
        electron_err_t_interpolated.append(interpolated_flux_ele * rel_error_ele)

        Error_Index_pos = positron_R.index(min(positron_R, key=lambda x: abs(x - R)))
        rel_error_pos = positron_fluxerror[Error_Index_pos] / positron_flux[Error_Index_pos]
        positron_err_t_interpolated.append(interpolated_flux_pos * rel_error_pos)
        '''

        '''
        # plot electron and positron vs rigidity for one example time bin,
        # and compare the two methods, and the match between model fit(s) and the actual data,
        # to make sure that our method works
        if timeindex==20 or timeindex==70:
        #if timeindex<200:
            plt.figure()    
            fig_ele, ax_ele, = plt.subplots(figsize=PythonPlotDefaultParameters.figsize)
            #plt.plot(xR, pl_approx_exact_ele(xR), '--', linewidth=1., color='black', label='local power-law approximation (exact method)')
            plt.plot(xR, pl_approx_rough_ele(xR), '--', linewidth=1., color='black', label='local power-law approximation')
            plt.errorbar(electron_R, electron_flux, yerr=electron_fluxerror, fmt='o', color='green', label='data points')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Rigidity (GV)')
            plt.ylabel(r'electron flux ($\mathrm{m}^{-2}\,\,\mathrm{s}^{-1}\,\,\mathrm{sr}^{-1}\,\,\mathrm{GV}^{-1}$)')
            plt.title( str(datetime.utcfromtimestamp(int(time_array[timeindex])).strftime('%Y-%m-%d')) + "_" + str(datetime.utcfromtimestamp(int(time_array[timeindex+1])).strftime('%Y-%m-%d')))
            plt.legend(prop={'size': 14})
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            ax_ele.tick_params(direction='in', which='major', length=10, width=5)
            ax_ele.tick_params(direction='in', which='minor', length=5, width=3)
            ax_ele.xaxis.set_minor_formatter(mticker.ScalarFormatter())
            YaxisMin = plt.ylim()[0]
            YaxisMax = plt.ylim()[1]
            plt.vlines(Emin, YaxisMin, YaxisMax, colors = "r", linestyles = "dashed")
            plt.vlines(Emax, YaxisMin, YaxisMax, colors = "r", linestyles = "dashed")
            ax_ele.set_ylim(YaxisMin, YaxisMax)
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "onebintoshow_"+ str(timeindex) + "_" + str(Emin) + "_" + str(Emax) + "_electron_flux.pdf")
            plt.close()

            plt.figure()    
            fig_pos, ax_pos, = plt.subplots(figsize=PythonPlotDefaultParameters.figsize)
            #plt.plot(xR, pl_approx_exact_pos(xR), '--', linewidth=1., color='black', label='local power-law approximation (exact method)')
            plt.plot(xR, pl_approx_rough_pos(xR), '--', linewidth=1., color='black', label='local power-law approximation')
            plt.errorbar(positron_R, positron_flux, yerr=positron_fluxerror, fmt='o', color='green', label='data points')
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Rigidity (GV)')
            plt.ylabel(r'positron flux ($\mathrm{m}^{-2}\,\,\mathrm{s}^{-1}\,\,\mathrm{sr}^{-1}\,\,\mathrm{GV}^{-1}$)')
            plt.title( str(datetime.utcfromtimestamp(int(time_array[timeindex])).strftime('%Y-%m-%d')) + "_" + str(datetime.utcfromtimestamp(int(time_array[timeindex+1])).strftime('%Y-%m-%d')))
            plt.legend(prop={'size': 14})
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            ax_pos.tick_params(direction='in', which='major', length=10, width=5)
            ax_pos.tick_params(direction='in', which='minor', length=5, width=3)
            ax_pos.xaxis.set_minor_formatter(mticker.ScalarFormatter())
            YaxisMin = plt.ylim()[0]
            YaxisMax = plt.ylim()[1]
            plt.vlines(Emin, YaxisMin, YaxisMax, colors = "r", linestyles = "dashed")
            plt.vlines(Emax, YaxisMin, YaxisMax, colors = "r", linestyles = "dashed")
            ax_pos.set_ylim(YaxisMin, YaxisMax)
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "onebintoshow_"+ str(timeindex) + "_" + str(Emin) + "_" + str(Emax) + "_positron_flux.pdf")
            plt.close()
        '''
    electronflux = unumpy.uarray(electron_flux_t_interpolated, electron_err_t_interpolated)
    positronflux = unumpy.uarray(positron_flux_t_interpolated, positron_err_t_interpolated)

    ratio = electronflux/positronflux

    return unumpy.nominal_values(ratio), unumpy.std_devs(ratio)



