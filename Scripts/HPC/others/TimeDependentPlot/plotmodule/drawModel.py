import draw
from draw import *
from scipy.interpolate import make_interp_spline
import time
import PythonPlotDefaultParameters

def drawModelAslam(PbarOverProAslamModel_2426_2480_1BR, PbarOverProAslamModel_2480_2506_6BR, PbarOverProAslamModel_2481_2520_1BR, index, rigidityrange, lowworkpath, intermediateworkpath, binmerge, yaxisfactor):

    if rigidityrange == "low":
        path = lowworkpath +"/totalall/results/plot/"
    elif rigidityrange == "intermediate":
        path = intermediateworkpath +"/total/results/plot/"

    acceptable_index = [1, 3, 5, 9]

    if index in acceptable_index:

        pointnumber          = PbarOverProAslamModel_2426_2480_1BR[:,index].shape[0] + 1
        Date_Model_all       = draw.dates_1Bartel[1:pointnumber] + draw.dates_1Bartel[pointnumber:][0:50:7]
        Date_Model_all_unixtime = [i.timestamp() for i in Date_Model_all]
        Prediction_Model_all      = np.concatenate((PbarOverProAslamModel_2426_2480_1BR[:,index]  , PbarOverProAslamModel_2480_2506_6BR[:,index])  , axis=0)
        PredictionError_Model_all = np.concatenate((PbarOverProAslamModel_2426_2480_1BR[:,index+1], PbarOverProAslamModel_2480_2506_6BR[:,index+1]), axis=0)

        #### Plot
        fig, ax = plt.subplots(figsize=(30, 18))

        ## plot BR 2426 to 2480 in 1BR
        plt.errorbar( draw.dates_1Bartel[1:pointnumber][0:46]+draw.dates_1Bartel[1:pointnumber][48:], np.append(PbarOverProAslamModel_2426_2480_1BR[:,index][0:46], PbarOverProAslamModel_2426_2480_1BR[:,index][48:])*100000, yerr=np.append(PbarOverProAslamModel_2426_2480_1BR[:,index+1][0:46], PbarOverProAslamModel_2426_2480_1BR[:,index+1][48:])*100000, marker='s', linestyle="None", markersize=12, markerfacecolor="b",ecolor="b")  ## first bartel rotation removed (2426) 

        ## plot BR 2480 to 2506 in 6BR 
        plt.errorbar( draw.dates_1Bartel[pointnumber:][0:50:7], PbarOverProAslamModel_2480_2506_6BR[:,index]*100000, yerr=PbarOverProAslamModel_2480_2506_6BR[:,index+1]*100000, marker='s', linestyle="None", markersize=12, markerfacecolor="r",ecolor="r")

        ## plot BR 2481 to 2520 in 1BR
        scaler_2481To2520 = PbarOverProAslamModel_2426_2480_1BR[:,index][-1]/PbarOverProAslamModel_2481_2520_1BR[:,index][0]
        print("scaler_2481To2520:" + str(scaler_2481To2520))
        plt.errorbar( draw.dates_1Bartel_Aslam[pointnumber: pointnumber+PbarOverProAslamModel_2481_2520_1BR[:,index].shape[0]], PbarOverProAslamModel_2481_2520_1BR[:,index]*100000*scaler_2481To2520, yerr=PbarOverProAslamModel_2481_2520_1BR[:,index+1]*100000, marker='s', linestyle="None", markersize=12, markerfacecolor="g",ecolor="g")  

        ## Plot all
        #plt.errorbar( Date_Model_all[0:46]+Date_Model_all[48:], np.append(Prediction_Model_all[0:46], Prediction_Model_all[48:])*100000, yerr=np.append(PredictionError_Model_all[0:46], PredictionError_Model_all[48:])*100000, marker='o', linestyle="None", markersize=20, markerfacecolor="r",ecolor="r")

        ax.tick_params(axis='x', labeltop=False, labelbottom=True , labelleft=False, labelright=False, direction='in', length=10, width=5)
        ax.tick_params(axis='y', labeltop=False, labelbottom=False, labelleft=True , labelright=False, direction='in', length=10, width=5)
        ymin = plt.ylim()[0]
        ymax = plt.ylim()[1]
        ax.axes.set_ylim([ ymin-((ymax-ymin)/2)*yaxisfactor,  ymax+((ymax-ymin)/2)*yaxisfactor ])

        plt.xticks(fontsize=50)
        plt.yticks(fontsize=60)
        #plt.ylabel("antiproton to proton relative flux ratio",fontsize=50)

        #plt.legend(loc='best',fontsize=50)

        #plt.text(-0.13, 0.7, r'$\rm{\overline{p}/p}$ ($\times$ $10^{-5}$)', transform = ax.transAxes, rotation=90, fontsize=60)
        plt.text(-0.13, 0.7, r'$\rm{\overline{p}/p}$', transform = ax.transAxes, rotation=90, fontsize=60)

        plt.savefig( path + "PbarOverProAslamModel" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf" )



def drawCompareWithModelAslam(PbarOverProAslamModel_2426_2480_1BR, PbarOverProAslamModel_2480_2506_6BR, PbarOverProAslamModel_2481_2520_1BR, index, rigidityrange, lowworkpath, intermediateworkpath, binmerge, yaxisfactor):

    if rigidityrange == "low":
        path = lowworkpath +"/totalall/results/plot/"
    elif rigidityrange == "intermediate":
        path = intermediateworkpath +"/total/results/plot/"

    acceptable_index = [1, 3, 5, 9]

    if index in acceptable_index:

        ## Scaler by mean value
        scaler = np.mean(draw.B6result[0:6]*100000) / np.mean(PbarOverProAslamModel_2426_2480_1BR[0:36,index]*100000)
        scaler_2481To2520 = PbarOverProAslamModel_2426_2480_1BR[:,index][-1]/PbarOverProAslamModel_2481_2520_1BR[:,index][0]

        pointnumber                = PbarOverProAslamModel_2426_2480_1BR[:,index].shape[0] + 1 # 56
        pointnumber_1B             = (PbarOverProAslamModel_2426_2480_1BR[:,index].shape[0] + 1) + PbarOverProAslamModel_2481_2520_1BR[:,index].shape[0] # 96
        Date_Model_all             = draw.dates_1Bartel[1:pointnumber] + draw.dates_1Bartel[pointnumber:][0:50:7]
        Date_Model_all_1B          = draw.dates_1Bartel_Aslam[1:pointnumber_1B]
        Date_Model_all_unixtime    = [i.timestamp() for i in Date_Model_all]
        Date_Model_all_unixtime_1B = [i.timestamp() for i in Date_Model_all_1B]

        Prediction_Model_all         = np.concatenate((PbarOverProAslamModel_2426_2480_1BR[:,index]*100000*scaler  , PbarOverProAslamModel_2480_2506_6BR[:,index]*100000*scaler), axis=0)
        PredictionError_Model_all    = np.concatenate((PbarOverProAslamModel_2426_2480_1BR[:,index+1]*100000*scaler, PbarOverProAslamModel_2480_2506_6BR[:,index+1]*100000*scaler), axis=0)
        Prediction_Model_all_1B      = np.concatenate((PbarOverProAslamModel_2426_2480_1BR[:,index]*100000*scaler  , PbarOverProAslamModel_2481_2520_1BR[:,index]*100000*scaler*scaler_2481To2520), axis=0)
        PredictionError_Model_all_1B = np.concatenate((PbarOverProAslamModel_2426_2480_1BR[:,index+1]*100000*scaler, PbarOverProAslamModel_2481_2520_1BR[:,index+1]*100000*scaler*scaler_2481To2520), axis=0)

        #### Plot
        fig, ax = plt.subplots()

        ## Plot Model: plot BR 2426 to 2480 in 1BR
        #plt.errorbar( draw.dates_1Bartel[1:pointnumber]       , PbarOverProAslamModel_2426_2480_1BR[:,index]*100000*scaler, yerr=0, marker='s', linestyle="None", markersize=20, markerfacecolor="r",ecolor="r") ## first bartel rotation removed (2426)

        ## Plot Model: plot BR 2480 to 2506 in 6BR
        #plt.errorbar( draw.dates_1Bartel[pointnumber:][0:50:7], PbarOverProAslamModel_2480_2506_6BR[:,index]*100000*scaler, yerr=0, marker='s', linestyle="None", markersize=20, markerfacecolor="r",ecolor="r")

        ## Plot Model: plot BR 2481 to 2520 in 1BR 
        #plt.errorbar( draw.dates_1Bartel_Aslam[pointnumber: pointnumber+PbarOverProAslamModel_2481_2520_1BR[:,index].shape[0]], PbarOverProAslamModel_2481_2520_1BR[:,index]*100000*scaler_2481To2520*scaler, yerr=PbarOverProAslamModel_2481_2520_1BR[:,index+1]*100000, marker='s', linestyle="None", markersize=12, markerfacecolor="g",ecolor="g")

        ## Plot Model: All time Range
        #plt.errorbar( Date_Model_all[0:46]+Date_Model_all[48:], np.append(Prediction_Model_all[0:46], Prediction_Model_all[48:]), yerr=np.append(PredictionError_Model_all[0:46], PredictionError_Model_all[48:]), marker='o', linestyle="None", markersize=20, markerfacecolor="r",ecolor="r", label='Model Prediction')

        ## Plot Spline line for model
        x_smooth_unixtime = draw.dates_1Bartel_unixsecond_Aslam[1:97] 
        y_smooth = make_interp_spline(Date_Model_all_unixtime_1B[0:46]+Date_Model_all_unixtime_1B[48:], np.append(Prediction_Model_all_1B[0:46], Prediction_Model_all_1B[48:]) )(x_smooth_unixtime)
        x_smooth = [datetime.fromtimestamp(i) for i in x_smooth_unixtime]
        plt.plot(x_smooth, y_smooth, color='red', linewidth=10, label='Model Prediction')

        ## Plot this analysis
        plt.errorbar(draw.dates6B, draw.B6result*100000, yerr=draw.Error_Total6B*100000, marker='o', linestyle="None", markersize=30, markerfacecolor="blue", ecolor="blue", linewidth=5, label='This analysis' )

        ymin = plt.ylim()[0]
        ymax = plt.ylim()[1]
        ax.axes.set_ylim([ ymin-((ymax-ymin)/2)*(yaxisfactor-1.0),  ymax+((ymax-ymin)/2)*(yaxisfactor-1.0) ])

        plt.legend()

        plt.text(-0.13, 0.7, r'$\rm{\overline{p}/p}$ ($\times$ $10^{5}$)', transform = ax.transAxes, rotation=90, fontsize=60)        

        plt.savefig( path + "CompareWithPbarOverProAslamModel" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf" )




