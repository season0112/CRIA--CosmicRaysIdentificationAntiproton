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


def Gauss_function(x, a, mean_gauss, sigma_gauss):
    return a * (1/(sigma_gauss*np.sqrt(2*math.pi))) * np.exp(-(x-mean_gauss)**2/(2*sigma_gauss**2))

def Novosibirsk_function(x, b, miu_Novo, sigma_Novo, tau_Novo):
    lamda = 1 + tau_Novo * (x-miu_Novo) * (np.sinh(tau_Novo*np.sqrt(np.log(4))))/(sigma_Novo*tau_Novo*np.sqrt(np.log(4))) # np.log(value, base). Defalt: np.log(4) = math.log(4, e).
    #print("lamda is "+ str(lamda))
    expindex = (-1/2) * ( ((np.log(lamda))**2)/(tau_Novo**2) + tau_Novo**2)
    #print( "expindex is " + str(expindex) )
    return b * np.exp( expindex )

def Novosibirsk_function_ACsoft(x, b, miu_Novo, sigma_Novo, tau_Novo):
    lamda = max(10**(-20), 1 + tau_Novo * (x-miu_Novo) * (np.sinh(tau_Novo*np.sqrt(np.log(4)))) / (sigma_Novo * tau_Novo * np.sqrt(np.log(4))) ) # np.log(value, base). Defalt: np.log(4) = math.log(4, e).
    #print("lamda is "+ str(lamda))
    return max(0.0, b * np.exp( (-1/2) * ( ((np.log(lamda))**2)/(tau_Novo**2) + tau_Novo**2) ))

def FitFunction(x, a, b, mean_gauss, sigma_gauss, miu_Novo, sigma_Novo, tau_Novo):
    #return Gauss_function(x, a, mean_gauss, sigma_gauss) + Novosibirsk_function(b, x, miu_Novo, sigma_Novo, tau_Novo)  ## this way has bug, only gaussion part works.
    return a * (1/(sigma_gauss*np.sqrt(2*math.pi))) * np.exp(-(x-mean_gauss)**2/(2*sigma_gauss**2)) + b * np.exp( (-1/2) * ( ((np.log(1 + tau_Novo * (x-miu_Novo) * (np.sinh(tau_Novo*np.sqrt(np.log(4))))/(sigma_Novo*tau_Novo*np.sqrt(np.log(4)))))**2)/(tau_Novo**2) + tau_Novo**2) ) 


def ExponentialFun(x, k1, k2, b, c):
    return k1 * np.exp(k2*x+b) + c

def func_None():
    print("Cannot find functions.")
