#!/usr/bin/env python 
import numpy as np
import os
from root_numpy import array2root, root2array
import argparse
import binning
import gc


def LoadBartelTimeStamp():
    with open ( os.getenv('MY_ANALYSIS') + "/Scripts/HPC/others/time.txt") as f:
      lines = []
      for line in f.readlines():
        line = line.strip("\n")
        lines.append(str(line))
    return lines

def LoadSixMonthsTimeStamp():
    with open (os.getenv('MY_ANALYSIS') + "/Scripts/HPC/others/time_6months.txt") as f:
      lines = []
      for line in f.readlines():
        line = line.strip("\n")
        lines.append(str(line))
    return lines


def main():
    ## index0 is bartel rotation 2426. pass7ext ends at 2530, which is index 104.
    if arguments.issversion == "pass7":
        Bartel_Rotation_Max = 95
        Six_Monthes_Max     = 15  # extension1: 0-14, 15 is 47%.

    elif arguments.issversion == "pass7ext": # to 29.10.2019
        Bartel_Rotation_Max = 117 # extension1:last bartel rotation is 102-104: 102,103 is okay, 104 is 66.7%.  extension2: to 116, 117 is empty.  115-116:10/18/2019-11/14/2019, 116-117:11/14/2019-12/11/2019.
        Six_Monthes_Max     = 19  # extension2: to 19, 18-19 is just a little bit.

    elif arguments.issversion == "pass7ext2": # to 06.06.2020 (from 29.10.2019 to 27.01.2020: ttcs off)
        Bartel_Rotation_Max = 126 # 123:05/21/2020, 124:06/17/2020, 126: 06/09/2020
        Six_Monthes_Max     = 20  # 18-19：10/04/2019-04/01/2020, 19-20:04/01/2020-09/28/2020 

    elif arguments.issversion == "pass7ext3": # to 03.05.2021 
        Bartel_Rotation_Max = 138 ## 135:07/05/2021, 136,137,138 should be empty.
        Six_Monthes_Max     = 21  ##  20:27/03/2021, 21:23/09/2021


    for rigidity_edge in range(arguments.Rigidity_start,arguments.Rigidity_end):

        if rigidity_edge == 28:
            rigidity = "{:g}".format(binning.published2016binnings[rigidity_edge]) + "_" + "{:g}".format(binning.published2016binnings[rigidity_edge+1]) # for intermediate range 
        else:
            rigidity = str(binning.published2016binnings[rigidity_edge]) + "_" + str(binning.published2016binnings[rigidity_edge+1])  # for low range
        

        print("begin rigidity:" + rigidity)

        negative = root2array(workpath + "/" + Dirname + "/B1130_pass7_7.8_all_Tree_negative_" + rigidity + ".root", treename)
        positive = root2array(workpath + "/" + Dirname + "/B1130_pass7_7.8_all_Tree_positive_" + rigidity + ".root", treename)

        #positive_sorted = positive[np.argsort(positive["TimeStamp"])]   ## np.argsort(positive["TimeStamp"]): sorted index of TimeStamp, 0,1,2,3,4,5,6.............
        #negative_sorted = negative[np.argsort(negative["TimeStamp"])]

        ####
        if arguments.mode == "Bartel":
            lines = LoadBartelTimeStamp()

            for i in range(0,Bartel_Rotation_Max):
                print("begin processing positive:" + str(i))
                index = np.where((int(lines[i]) < positive["TimeStamp"]) & (positive["TimeStamp"] <= int(lines[i+1])))[0]     
                np.save(workpath + "/" + Dirname + "/numpyfiles/positive/positive_" + rigidity + "_" +str(i)+".npy",positive[index])
                array2root(positive[index], workpath + "/" + Dirname + "/rootfiles/positive/positive_" + rigidity + "_" +str(i)+".root", 'tree','recreate')

            for i in range(0,Bartel_Rotation_Max):
                index = np.where((int(lines[i]) < negative["TimeStamp"]) & (negative["TimeStamp"] <= int(lines[i+1])))[0]
                np.save(workpath + "/" + Dirname + "/numpyfiles/negative/negative_" + rigidity + "_" +str(i)+".npy",negative[index])
                array2root(negative[index], workpath + "/" + Dirname + "/rootfiles/negative/negative_" + rigidity + "_" +str(i)+".root", 'tree','recreate')

        elif arguments.mode == "6months":
            lines = LoadSixMonthsTimeStamp()

            for i in range(0,Six_Monthes_Max):
                print("begin processing positive:" + str(i))
                index = np.where((int(lines[i]) < positive["TimeStamp"]) & (positive["TimeStamp"] <= int(lines[i+1])))[0]
                np.save(workpath + "/" + Dirname + "/numpyfiles_6months/positive/positive_" + rigidity + "_" +str(i)+".npy",positive[index])
                array2root(positive[index], workpath + "/" + Dirname + "/rootfiles_6months/positive/positive_" + rigidity + "_" +str(i)+".root", 'tree','recreate')

            for i in range(0,Six_Monthes_Max):
                index = np.where((int(lines[i]) < negative["TimeStamp"]) & (negative["TimeStamp"] <= int(lines[i+1])))[0]
                np.save(workpath + "/" + Dirname + "/numpyfiles_6months/negative/negative_" + rigidity + "_" +str(i)+".npy",negative[index])
                array2root(negative[index], workpath + "/" + Dirname + "/rootfiles_6months/negative/negative_" + rigidity + "_" +str(i)+".root", 'tree','recreate')
        else:
             print("mode choice wrong")
             sys.exit()

        #### Release Merory 
        del negative, positive 
        gc.collect()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-M','--mode',help='which split mode you need, Bartel Rotations or 6 months')
    parser.add_argument('-C','--cluster',help='which cluster')
    parser.add_argument('-V','--issversion',help='which issversion you choose')
    parser.add_argument('-S','--Rigidity_start', type=int, help='Rigidity_start')
    parser.add_argument('-E','--Rigidity_end', type=int, help='Rigidity_end')
    parser.add_argument('-R','--rigidityrange', help='Rigidity_end')
    arguments = parser.parse_args()

    juamsintermediateworkpath = os.getenv('JUAMSINTERMEDIATEENERGYDIR')
    hpcintermediateworkpath = os.getenv('HPCINTERMEDIATEDIR')
    hpclowworkpath = os.getenv('HPCLOWENERGYDIR')

    if arguments.cluster == "HPC":
        intermediateworkpath = hpcintermediateworkpath
        lowworkpath = hpclowworkpath
    elif arguments.cluster == "JUAMS":
        intermediateworkpath = juamsintermediateworkpath

    if arguments.rigidityrange == "IntermediateRange":
        workpath = intermediateworkpath
        Dirname = "total"
        treename = "AntiprotonIntermediateEnergyTree"
    elif arguments.rigidityrange == "LowRange":
        workpath = lowworkpath
        Dirname = "totalall"
        treename = "AntiprotonLowEnergyTree"

    main()








