import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick 
from scipy.signal import savgol_filter
import matplotlib
import PythonPlotDefaultParameters


def Plot_Miu_of_Gaussion(x, y, yerror, value, RigidityBin, Dimension, yaxisfactor, rigidityrange):

    fig_mine, ax = plt.subplots()

    if value == "Antiproton":
        color = 'blue'
        plt.ylabel(r"Mean $\mu_{\overline{\rm{p}}," + Dimension + "}$", horizontalalignment='right', y=1.0)
        if rigidityrange == 'low': 
            ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y', useMathText=True)
            ax.yaxis.get_offset_text().set_fontsize(60)
    elif value == "Electron":
        color = 'red'
        plt.ylabel(r"Mean $\mu_{\rm{e}^{-}," + Dimension + "}$"       , horizontalalignment='right', y=1.0) 
    elif value == "Pion":
        color = 'green'
        plt.ylabel(r"Mean $\mu_{\rm{secondaries}," + Dimension + "}$" , horizontalalignment='right', y=1.0)

    plt.errorbar(x, -np.array(y), yerr=yerror, marker='o', linestyle="None", markerfacecolor=color, ecolor=color, markeredgecolor=color, markersize=30)

    if plt.ylim()[0]>0:
        ax.set_ylim(plt.ylim()[0]*(1-(yaxisfactor-1)), plt.ylim()[1]*yaxisfactor)
    elif plt.ylim()[0]<0: 
        ax.set_ylim(plt.ylim()[0]*yaxisfactor, plt.ylim()[1]*(1-(yaxisfactor-1)))

    plt.savefig( value + str("_") + Dimension + str("_") + str(RigidityBin) + "_miu.pdf")
    plt.close()


def Plot_Sigma_of_Gaussion(x, y, yerror, value, RigidityBin, Dimension, yaxisfactor, rigidityrange):

    fig_mine, ax = plt.subplots()

    if value == "Antiproton":
        color = 'blue'
        plt.ylabel(r"Sigma $\sigma_{\overline{\rm{p}}," + Dimension + "}$"   , horizontalalignment='right', y=1.0)
        if rigidityrange == 'low':
            yaxisfactor = 1.05
    elif value == "Electron":
        color = 'red'
        plt.ylabel(r"Sigma $\sigma_{\rm{e}^{-}," + Dimension + "}$"          , horizontalalignment='right', y=1.0)
    elif value == "Pion":
        color = 'green'
        plt.ylabel(r"Sigma $\sigma_{\rm{secondaries}," + Dimension + "}$", horizontalalignment='right', y=1.0)
        if rigidityrange == 'low':
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))

    plt.errorbar(x, y, yerr=yerror, marker='o', linestyle="None", markersize=30, markerfacecolor=color, ecolor=color, markeredgecolor=color)

    if plt.ylim()[0]>0:
        ax.set_ylim(plt.ylim()[0]*(1-(yaxisfactor-1)), plt.ylim()[1]*yaxisfactor)
    elif plt.ylim()[0]<0:
        ax.set_ylim(plt.ylim()[0]*yaxisfactor, plt.ylim()[1]*(1-(yaxisfactor-1)))

    plt.savefig( value + str("_") + Dimension + str(RigidityBin) + "_sigma.pdf")
    plt.close()


def Plot_Chi2(x, y, name, RigidityRange, RigidityBin):

    plt.figure(figsize=(18,9))

    plt.plot(x, y, "*", markersize=10)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('Rigidity (GV)', fontsize=30, fontname='Arial')
    plt.ylabel('Chi2/dof'     , fontsize=30)

    #plt.ylim(-1, 500)

    plt.savefig(name + RigidityRange + str(RigidityBin) + str(".pdf"))
    plt.close()


def Plot_Template(BinCenter, y_fit, reduced_chi_squared_TRD, BinContent, value, RigidityBin, time_index, Dimension, fit_start, fit_end): 
    fig_mine, ax = plt.subplots()

    plt.errorbar(BinCenter, y_fit      , yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="red" , ecolor="red" , linewidth=5, label="Fit"+"(Chi2/ndf="+str(round(reduced_chi_squared_TRD,2)) +")" )
    plt.errorbar(BinCenter, BinContent, yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="blue", ecolor="blue", linewidth=5, label="Data")

    plt.vlines(BinCenter[fit_start], 0, np.max(y_fit), color='black')
    plt.vlines(BinCenter[fit_end-1  ], 0, np.max(y_fit), color='black')

    plt.legend( loc='best', fontsize=40 )
    plt.savefig(str("Template_") + Dimension + str("_") + value + str("_") + str(RigidityBin) + str("_") + str(time_index) + ".pdf")
    plt.yscale('log')
    #ax.axes.set_ylim([0.1, plt.ylim()[1] ])
    plt.savefig(str("Template_") + Dimension + str("_") + value + str("_") + str(RigidityBin) + str("_") + str(time_index) + "_LogY.pdf")
    plt.close()


def Plot_TemplatetFitPlot(BinCenter, x_more, y_value_antiproton, y_value_electron, y_fit_antiproton, y_fit_electron, reduced_chi2_antiproton, reduced_chi2_electron, RigidityBin, time_index, Dimension):

    BinCenter = -BinCenter
    x_more    = -x_more

    '''
    print("shape:")
    print(BinCenter.shape)
    print(BinCenter)
    BinWidth = BinCenter[0]-BinCenter[1] 
    BinEdge_low   = BinCenter[-1]-BinWidth/2
    BinEdge_high  = BinCenter[0]+BinWidth/2
    BinEdge = np.arange(BinEdge_low, BinEdge_high+BinWidth, BinWidth)
    print(y_value_electron.shape)
    print(BinEdge.shape)
    print(BinEdge)
    h_electron   = np.histogram(y_value_electron  , bins=BinEdge)
    h_antiproton = np.histogram(y_value_antiproton, bins=BinEdge)
    '''

    # Plot
    fig, ax = plt.subplots()

    scaler_antiproton = 1./max(y_value_antiproton)
    scaler_electron   = 1./max(y_value_electron)
    yyyy= y_value_electron  *scaler_electron

    #plt.errorbar(BinCenter, y_value_electron  *scaler_electron  , yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="red"  , ecolor="red" , linewidth=5, label="Electron" )
    plt.step(BinCenter    , y_value_electron  *scaler_electron  , where='mid', linewidth=8, color='red'          , label="Electron"  )
    plt.plot(       x_more, y_fit_electron    *scaler_electron  , linewidth=7, color='red', linestyle=(0, (5, 1)), label="Electron Fit (Chi2/ndf=" + str(round(reduced_chi2_electron,2)) + ")" ) # linestyle: 'densely dashed'=(0, (5, 1))

    #plt.errorbar(BinCenter, y_value_antiproton*scaler_antiproton, yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="blue" , ecolor="blue", linewidth=5, label="Antiproton" )
    plt.step(BinCenter    , y_value_antiproton*scaler_antiproton, where='mid', linewidth='8', color='blue'         , label="Antiproton"  )
    plt.plot(       x_more, y_fit_antiproton  *scaler_antiproton, linewidth=7, color='blue' , linestyle=(0, (5, 1)), label="Antiproton Fit (Chi2/ndf=" + str(round(reduced_chi2_antiproton,2)) + ")" ) 

    handles, labels = plt.gca().get_legend_handles_labels()
    order = [0,2,1,3]
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

    plt.xlabel( r'$\Lambda_{\bf{TRD}}$', horizontalalignment='right', x=1.0)  
    plt.ylabel( 'Normalized events'    , horizontalalignment='right', y=1.0)  

    plt.savefig(str("TemplateFit_") + Dimension + str("_") + str(RigidityBin) + str("_") + str(time_index) + ".pdf")


def GetFileName(rigidityrange, value):
    if rigidityrange == 'intermediate':
        Dimension = 'TRD'
        if value == "Antiproton":
            FileName = 'data_pass7_positive_TRDeff_1.00'
        elif value == "Electron":
            FileName = 'template_electron_Data_TRDeff_1.00'
        elif value == "Pion":
            FileName = 'template_pion_Data_TRDeff_1.00'
    elif rigidityrange == 'low':
        Dimension = 'TOF'
        if value == "Antiproton":
            FileName = 'TofTRD_data_pass7_positive_template_TRDeff_0.98_TOFeff_0.97'
        elif value == "Electron":
            FileName = 'TofTRD_template_ElectronData_TRDeff_0.98_TOFeff_0.97'
        elif value == "Pion":
            FileName = 'TofTRD_template_PionData_TRDeff_0.98_TOFeff_0.97'
    return FileName, Dimension


def TunningForFitSettings(rigidityrange, value, UsedFunction, R_index, time_index, InitialValue, fit_start, fit_end):
    if rigidityrange == 'intermediate':
        print("progress")

    elif rigidityrange == 'low':
        if value == "Antiproton":
            print("progress")
        elif value == "Electron":
            if UsedFunction =='Gaussion':
                if R_index == 5:
                    if time_index == 96:
                        InitialValue = [100, -0.028, 0.02]
                        fit_start = 2
                        fit_end   = 17
                    if time_index == 102:
                        InitialValue = [100, -0.028, 0.02]
                        fit_start = 2
                        fit_end   = 17
                    if time_index == 120:
                        InitialValue = [80, -0.028, 0.02]
                        fit_start = 2
                        fit_end   = 13
                if R_index == 6:
                    if time_index == 12:
                        InitialValue = [220, -0.01, 0.02]
                        fit_start = 4
                        fit_end   = 15
                    if time_index == 42:
                        InitialValue = [140, -0.02, 0.02]
                        fit_start = 5
                        fit_end   = 15
                    if time_index == 90:
                        InitialValue = [175, -0.02, 0.02]
                        fit_start = 4
                        fit_end   = 16
                    if time_index == 96:
                        InitialValue = [110, -0.02, 0.02]
                        fit_start = 2
                        fit_end   = 17
                    if time_index == 102:
                        InitialValue = [110, -0.02, 0.02]
                        fit_start = 1
                        fit_end   = 17
                    if time_index == 114:
                        InitialValue = [80, -0.02, 0.02]
                        fit_start = 1
                        fit_end   = 17
        elif value == "Pion":
            if UsedFunction =='Gaussion':
                if R_index == 5:
                    if time_index == 6:
                        InitialValue = [500, -0.025, 0.02]
                        fit_start = 2
                        fit_end   = 19
                    if time_index == 30:
                        InitialValue = [350, -0.025, 0.02]
                        fit_start = 2
                        fit_end   = 19


    return InitialValue, fit_start, fit_end


def InitialValue(rigidityrange, value, UsedFunction, R_index):

    if rigidityrange == 'intermediate':
        Dimension = 'TRD'
        if value == "Antiproton":
            FileName = 'data_pass7_positive_TRDeff_1.00'
            if UsedFunction =='Gaussion':
                InitialValue = [1000000, -1.1, 0.2]  
            elif UsedFunction =='Novosibirsk':
                InitialValue = [1000000, -1.1, 0.2, 0.01] 
            fit_start = 18
            fit_end   = 30
        elif value == "Electron":
            FileName = 'template_electron_Data_TRDeff_1.00'
            if UsedFunction =='Gaussion':
                InitialValue = [3000, -0.3, 0.1] 
            elif UsedFunction =='Novosibirsk':
                InitialValue = [3000, -0.3, 0.1, 0.01]  
            fit_start = 22 #25
            fit_end   = -1
        elif value == "Pion":
            FileName = 'template_pion_Data_TRDeff_1.00'

    elif rigidityrange == 'low':
        Dimension = 'TOF'
        if value == "Antiproton":
            FileName = 'TofTRD_data_pass7_positive_template_TRDeff_0.98_TOFeff_0.97'
            if UsedFunction =='Gaussion':
                InitialValue = [1000000, 0, 0.2] 
            elif UsedFunction =='Novosibirsk':
                InitialValue = [200000, 0, 0.2, 0.01] 
            fit_start = 5
            fit_end   = 18
        elif value == "Electron":
            FileName = 'TofTRD_template_ElectronData_TRDeff_0.98_TOFeff_0.97'
            if UsedFunction =='Gaussion':
                if R_index < 7:
                    InitialValue = [200, -0.05, 0.02]  
                    fit_start = 0
                    fit_end   = 9
                elif R_index == 7:
                    InitialValue = [175, 0.00, 0.02]
                    fit_start = 0
                    fit_end   = 22

            elif UsedFunction =='Novosibirsk':
                if R_index == 0:
                    if time_index < 10:
                        InitialValue = [20, -0.15, 0.02, 0.01] 
                        fit_start = 0
                        fit_end   = 9
                    else:
                        InitialValue = [20, -0.13, 0.01, 0.01] 
                        fit_start = 0
                        fit_end   = 9
                elif 0 < R_index < 5:
                    InitialValue = [200, -0.04, 0.05, 0.01]  
                    fit_start = 0
                    fit_end   = 20
                else:
                    InitialValue = [200, -0.03, 0.02, 0.01]  
                    fit_start = 0
                    fit_end   = 20
        elif value == "Pion":
            FileName = 'TofTRD_template_PionData_TRDeff_0.98_TOFeff_0.97'
            if UsedFunction =='Gaussion':
                InitialValue = [400, -0.05, 0.05] 
                fit_start = 0
                fit_end   = 20
            elif UsedFunction =='Novosibirsk':
                InitialValue = [400, -0.05, 0.05, 0.01] 
                fit_start = 0
                fit_end   = 20

    return FileName, Dimension, InitialValue, fit_start, fit_end




