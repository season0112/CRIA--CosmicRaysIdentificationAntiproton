#!/usr/bin/env python
from root_numpy import root2array, tree2array, array2root
from ROOT import TFile
import numpy as np
import os
import argparse


def main():

    if Sign == "positive":
        positive = root2array(str(rootfile), Treename)
        os.chdir(os.path.abspath('..'))

        #### Split data in Rigidity bins in all time range
        for i in range( np.int(StartIndex), np.int(EndIndex)):
            positive_tem = positive[np.where((binning[i]<positive['Rigidity']) & (positive['Rigidity']<binning[i+1]))[0]]
            if (binning[i] == 93 or binning[i] == 108 or binning[i] == 125 or binning[i] == 147 or binning[i] == 175 or binning[i] == 211 or binning[i] == 259 or binning[i] == 450):
                array2root(positive_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "{:.0f}".format(binning[i]) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
            elif (binning[i] == 80.5):
                array2root(positive_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ str(round(binning[i],2)) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
            else:
                array2root(positive_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ str(round(binning[i],2)) + "_" + str(round(binning[i+1],2)) + '.root', Treename, 'recreate')
        if arguments.rigidityrange == "highrange":
            positive_tem = positive[np.where((211<positive['Rigidity']) & (positive['Rigidity']<250))[0]]
            array2root(positive_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "211" + "_" + "250" + '.root', Treename, 'recreate')
            positive_tem = positive[np.where((250<positive['Rigidity']) & (positive['Rigidity']<330))[0]]
            array2root(positive_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "250" + "_" + "330" + '.root', Treename, 'recreate')
            positive_tem = positive[np.where((259<positive['Rigidity']) & (positive['Rigidity']<330))[0]]
            array2root(positive_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "259" + "_" + "330" + '.root', Treename, 'recreate')
            positive_tem = positive[np.where((330<positive['Rigidity']) & (positive['Rigidity']<525))[0]]
            array2root(positive_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "330" + "_" + "525" + '.root', Treename, 'recreate')
        
        #### Split ISS data in Rigidity bins in PRL time range
        if arguments.datamode == "MC":
            print("the dataset is MC. No Timestamp process will be done.")
        elif arguments.datamode == "ISSdata":
            positive_May2015 = positive[np.where(positive['TimeStamp']<PRL_splitdata)[0]]
            array2root(positive_May2015, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015.root', Treename, 'recreate')
            for i in range(np.int(StartIndex), np.int(EndIndex)):
                positive_tem2 = positive_May2015[np.where((binning[i]<positive_May2015['Rigidity']) & (positive_May2015['Rigidity']<binning[i+1]))[0]]
                if (binning[i] == 93 or binning[i] == 108 or binning[i] == 125 or binning[i] == 147 or binning[i] == 175 or binning[i] == 211 or binning[i] == 259 or binning[i] == 450):
                    array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + "{:.0f}".format(binning[i]) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                elif (binning[i] == 80.5):
                    array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + str(round(binning[i],2)) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                else:
                    array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + str(round(binning[i],2)) + "_" + str(round(binning[i+1],2)) + '.root', Treename, 'recreate')
            if arguments.rigidityrange == "highrange":
                positive_tem2 = positive_May2015[np.where((211<positive_May2015['Rigidity']) & (positive_May2015['Rigidity']<250))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + "211" + "_" + "250" + '.root', Treename, 'recreate')
                positive_tem2 = positive_May2015[np.where((250<positive_May2015['Rigidity']) & (positive_May2015['Rigidity']<330))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + "250" + "_" + "330" + '.root', Treename, 'recreate')
                positive_tem2 = positive_May2015[np.where((259<positive_May2015['Rigidity']) & (positive_May2015['Rigidity']<330))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + "259" + "_" + "330" + '.root', Treename, 'recreate')
                positive_tem2 = positive_May2015[np.where((330<positive_May2015['Rigidity']) & (positive_May2015['Rigidity']<525))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + "330" + "_" + "525" + '.root', Treename, 'recreate')
         
        ##### Split ISS data in Rigidity bins in PhysReport time range
        if arguments.datamode == "MC":
            print("the dataset is MC. No Timestamp process will be done.")
        elif arguments.datamode == "ISSdata":
            positive_Nov2017 = positive[np.where(positive['TimeStamp']<PhyRep_splitdata)[0]]
            array2root(positive_Nov2017, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017.root', Treename, 'recreate')
            for i in range(np.int(StartIndex), np.int(EndIndex)):
                positive_tem2 = positive_Nov2017[np.where((binning[i]<positive_Nov2017['Rigidity']) & (positive_Nov2017['Rigidity']<binning[i+1]))[0]]
                if (binning[i] == 93 or binning[i] == 108 or binning[i] == 125 or binning[i] == 147 or binning[i] == 175 or binning[i] == 211 or binning[i] == 259 or binning[i] == 450):
                    array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + "{:.0f}".format(binning[i]) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                elif (binning[i] == 80.5):
                    array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + str(round(binning[i],2)) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                else:
                    array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + str(round(binning[i],2)) + "_" + str(round(binning[i+1],2)) + '.root', Treename, 'recreate')
            if arguments.rigidityrange == "highrange":
                positive_tem2 = positive_Nov2017[np.where((211<positive_Nov2017['Rigidity']) & (positive_Nov2017['Rigidity']<250))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + "211" + "_" + "250" + '.root', Treename, 'recreate')
                positive_tem2 = positive_Nov2017[np.where((250<positive_Nov2017['Rigidity']) & (positive_Nov2017['Rigidity']<330))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + "250" + "_" + "330" + '.root', Treename, 'recreate')
                positive_tem2 = positive_Nov2017[np.where((259<positive_Nov2017['Rigidity']) & (positive_Nov2017['Rigidity']<330))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + "259" + "_" + "330" + '.root', Treename, 'recreate')
                positive_tem2 = positive_Nov2017[np.where((330<positive_Nov2017['Rigidity']) & (positive_Nov2017['Rigidity']<525))[0]]
                array2root(positive_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + "330" + "_" + "525" + '.root', Treename, 'recreate')


    elif Sign == "negative":
        negative = root2array(str(rootfile), Treename)
        os.chdir(os.path.abspath('..'))
        
        #### Split data in Rigidity bins in all time range
        for i in range( np.int(StartIndex), np.int(EndIndex)):
            negative_tem = negative[np.where((-binning[i]>negative['Rigidity']) & (negative['Rigidity']>-binning[i+1]))[0]]
            if (binning[i] == 93 or binning[i] == 108 or binning[i] == 125 or binning[i] == 147 or binning[i] == 175 or binning[i] == 211 or binning[i] == 259 or binning[i] == 450):
                array2root(negative_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "{:.0f}".format(binning[i]) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
            elif (binning[i] == 80.5):
                array2root(negative_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ str(round(binning[i],2)) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
            else:
                array2root(negative_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ str(round(binning[i],2)) + "_" + str(round(binning[i+1],2)) + '.root', Treename, 'recreate')
        if arguments.rigidityrange == "highrange":
            negative_tem = negative[np.where((-211>negative['Rigidity']) & (negative['Rigidity']>-250))[0]]
            array2root(negative_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "211" + "_" + "250" + '.root', Treename, 'recreate')
            negative_tem = negative[np.where((-250>negative['Rigidity']) & (negative['Rigidity']>-330))[0]]
            array2root(negative_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "250" + "_" + "330" + '.root', Treename, 'recreate')
            negative_tem = negative[np.where((-259>negative['Rigidity']) & (negative['Rigidity']>-330))[0]]
            array2root(negative_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "259" + "_" + "330" + '.root', Treename, 'recreate')
            negative_tem = negative[np.where((-330>negative['Rigidity']) & (negative['Rigidity']>-525))[0]]
            array2root(negative_tem, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_all_'+ "330" + "_" + "525" + '.root', Treename, 'recreate')

        #### Split ISS data in Rigidity bins in PRL time range
        if arguments.datamode == "MC":
            print("the dataset is MC. No Timestamp process will be done.")
        elif arguments.datamode == "ISSdata":
            negative_May2015 = negative[np.where(negative['TimeStamp']<PRL_splitdata)[0]]
            array2root(negative_May2015, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015.root', Treename, 'recreate')
            for i in range(np.int(StartIndex), np.int(EndIndex)):
                negative_tem2 = negative_May2015[np.where((-binning[i]>negative_May2015['Rigidity']) & (negative_May2015['Rigidity']>-binning[i+1]))[0]]
                if (binning[i] == 93 or binning[i] == 108 or binning[i] == 125 or binning[i] == 147 or binning[i] == 175 or binning[i] == 211 or binning[i] == 259 or binning[i] == 450):
                    array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + "{:.0f}".format(binning[i]) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                elif (binning[i] == 80.5):
                   array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + str(round(binning[i],2)) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                else:
                    array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015' + '_' + str(round(binning[i],2)) + "_" + str(round(binning[i+1],2)) + '.root', Treename, 'recreate')
            if arguments.rigidityrange == "highrange":  
                negative_tem2 = negative_May2015[np.where((-211>negative_May2015['Rigidity']) & (negative_May2015['Rigidity']>-250))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015_'+ "211" + "_" + "250" + '.root', Treename, 'recreate')
                negative_tem2 = negative_May2015[np.where((-250>negative_May2015['Rigidity']) & (negative_May2015['Rigidity']>-330))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015_'+ "250" + "_" + "330" + '.root', Treename, 'recreate')
                negative_tem2 = negative_May2015[np.where((-259>negative_May2015['Rigidity']) & (negative_May2015['Rigidity']>-330))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015_'+ "259" + "_" + "330" + '.root', Treename, 'recreate')
                negative_tem2 = negative_May2015[np.where((-330>negative_May2015['Rigidity']) & (negative_May2015['Rigidity']>-525))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_May2015_'+ "330" + "_" + "525" + '.root', Treename, 'recreate')
     

        #### Split ISS data in Rigidity bins in PhysReport time range
        if arguments.datamode == "MC":
            print("the dataset is MC. No Timestamp process will be done.")
        elif arguments.datamode == "ISSdata":
            negative_Nov2017 = negative[np.where(negative['TimeStamp']<PhyRep_splitdata)[0]]
            array2root(negative_Nov2017, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017.root', Treename, 'recreate')
            for i in range(np.int(StartIndex), np.int(EndIndex)):
                negative_tem2 = negative_Nov2017[np.where((-binning[i]>negative_Nov2017['Rigidity']) & (negative_Nov2017['Rigidity']>-binning[i+1]))[0]]
                if (binning[i] == 93 or binning[i] == 108 or binning[i] == 125 or binning[i] == 147 or binning[i] == 175 or binning[i] == 211 or binning[i] == 259 or binning[i] == 450):
                    array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + "{:.0f}".format(binning[i]) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                elif (binning[i] == 80.5):
                    array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + str(round(binning[i],2)) + "_" + "{:.0f}".format(binning[i+1]) + '.root', Treename, 'recreate')
                else:
                    array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017' + '_' + str(round(binning[i],2)) + "_" + str(round(binning[i+1],2)) + '.root', Treename, 'recreate')
            if arguments.rigidityrange == "highrange":
                negative_tem2 = negative_Nov2017[np.where((-211>negative_Nov2017['Rigidity']) & (negative_Nov2017['Rigidity']>-250))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017_'+ "211" + "_" + "250" + '.root', Treename, 'recreate')
                negative_tem2 = negative_Nov2017[np.where((-250>negative_Nov2017['Rigidity']) & (negative_Nov2017['Rigidity']>-330))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017_'+ "250" + "_" + "330" + '.root', Treename, 'recreate')
                negative_tem2 = negative_Nov2017[np.where((-259>negative_Nov2017['Rigidity']) & (negative_Nov2017['Rigidity']>-330))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017_'+ "259" + "_" + "330" + '.root', Treename, 'recreate')
                negative_tem2 = negative_Nov2017[np.where((-330>negative_Nov2017['Rigidity']) & (negative_Nov2017['Rigidity']>-525))[0]]
                array2root(negative_tem2, "finebinroot/"+rootfile.split('/')[-1][:-5]+'_Nov2017_'+ "330" + "_" + "525" + '.root', Treename, 'recreate')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--path', help='current path')
    parser.add_argument('--rigidityrange', help='which rigidityrange you choose, high? intermediate?')
    parser.add_argument('--rigiditystart', type= int, help='which rigidity')
    parser.add_argument('--rigidityend', type = int, help='which rigidity')
    parser.add_argument('--inputfile', help='input')
    parser.add_argument('--datamode', help='the dataset is MC or ISS data')
    parser.add_argument('--rigiditysign', help='which rigidity sign')
    arguments = parser.parse_args()

    os.chdir(arguments.path)

    StartIndex = arguments.rigiditystart
    EndIndex   = arguments.rigidityend
    Sign       = arguments.rigiditysign


    binning = np.array([ 0.8, 1.0, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0, 14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93, 108, 125, 147, 175, 211, 259, 450])
    rootfile = arguments.inputfile

    if arguments.rigidityrange == "intermediaterange":
        Treename = "AntiprotonIntermediateEnergyTree"
    elif arguments.rigidityrange == "highrange":
        Treename = "ExampleAnalysisTree"
    elif arguments.rigidityrange == "lowrange":
        Treename = "AntiprotonLowEnergyTree"

    PRL_splitdata    = 1432598400 # May 26,      2015 12:00:00 AM
    PhyRep_splitdata = 1510484400 # November 12, 2017 11:00:00 AM

    main()





