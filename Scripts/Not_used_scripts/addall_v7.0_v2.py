#!/usr/bin/env python
import numpy as np
import argparse
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument('-D','--dataset', help='which_dataset_you_choose')
parser.add_argument('-E','--energy_range', choices=['14.1_38.9GeV', '38.9-147GeV', '147-1130GeV','1000-2000GeV','2.97-1130GeV','0.8_1130GeV'], help='which_energy_you_choose')
parser.add_argument('--cluster', help='which cluster you run your job.')
arguments = parser.parse_args()

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)

if arguments.cluster == "JUAMS":
    workpath = os.getenv('JUAMSHIGHENERGYDATADIR')
elif arguments.cluster == "HPC":
    workpath = "/hpcwork/jara0052/sichen/analysis_7.1_fixwithIGRFcut"
    workpath = os.getenv('HPCHIGHENERGYDATADIR')

begintime = time.time()

if arguments.energy_range == '147-1130GeV':
    binnings = np.array([147.0, 175.0, 211.0, 250.0, 330.0, 525.0, 643.0, 822.0, 1130.0])
elif arguments.energy_range == '38.9-147GeV':
    binnings = np.array([38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0])
elif arguments.energy_range == '14.1_38.9GeV':
    binnings = np.array([14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7,26.7,28.8,31.1,33.5,36.1,38.9])
elif arguments.energy_range == '1000-2000GeV':
    binnings = np.array([1000,1300,1600,2000])
elif arguments.energy_range == "2.97-1130GeV":
    binnings = np.array([2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0, 14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0, 175.0, 211.0, 250.0, 330.0, 525.0, 643.0, 822.0, 1130.0])
elif arguments.energy_range == "0.8_1130GeV":
    binnings = np.array([0.8, 1, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0, 14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0, 175.0, 211.0, 250.0, 330.0, 525.0, 643.0, 822.0, 1130.0])
else:
    print("no energy range is given")

trackerpattern = 0 

binleft = 0
path = workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' + '/results/rawdata/' + 'positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern)+'/'
print('filenum:',len([lists for lists in os.listdir(path) if os.path.isfile(os.path.join(path, lists))]))
number_of_files = len([lists for lists in os.listdir(path) if os.path.isfile(os.path.join(path, lists))])

#######################################################################################
for binleft in range(binnings.shape[0]-1):
    a0=np.array([])
    for i in range(0,number_of_files):
        a_tem = np.load(workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' + '/results/rawdata/' + 'positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern)+'/'+'positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_'+str(trackerpattern)+ "_" + "%05d" % i +'.npy')
        if a_tem.shape[0] != 0:
            if a0.shape[0] == 0:
                a0 = a_tem
            else:
                a0 = np.row_stack((a0,a_tem))
    if os.path.isdir(workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' +'/results/rawdata/' + 'transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern)):
        pass
    else:
        mkdir(workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' +'/results/rawdata/' + 'transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern))
    np.save(workpath+"/" + arguments.dataset + "/" + arguments.energy_range+ '_v2'  + '/results/rawdata/' + 'transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern)+'/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_'+str(trackerpattern)+'.npy',a0)

    a0 = np.array([])
    for i in range(0,number_of_files):
        a_tem = np.load(workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' + '/results/rawdata/' + 'negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern)+'/'+'negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_'+str(trackerpattern)+ "_" + "%05d" % i +'.npy')
        if a_tem.shape[0] != 0:
            if a0.shape[0] == 0:
                a0 = a_tem
            else:
                a0 = np.row_stack((a0,a_tem))
    if os.path.isdir(workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' +'/results/rawdata/' + 'transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern)):
        pass
    else:
        mkdir(workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' +'/results/rawdata/' + 'transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern))
    np.save(workpath+"/" + arguments.dataset + "/" + arguments.energy_range + '_v2' +'/results/rawdata/' + 'transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern'+str(trackerpattern)+'/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_'+str(trackerpattern)+'.npy',a0)

#### rebinning_259_JUAMS ###################################################################################

if (arguments.energy_range == '147-1130GeV' or arguments.energy_range == '0.8_1130GeV'):
    ##
    positive_211_250 = np.load(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/211.0_250.0GeV/pattern0/positive_211.0_250.0_pattern_0.npy')
    positive_250_330 = np.load(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/250.0_330.0GeV/pattern0/positive_250.0_330.0_pattern_0.npy')
    negative_211_250 = np.load(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/211.0_250.0GeV/pattern0/negative_211.0_250.0_pattern_0.npy')
    negative_250_330 = np.load(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/250.0_330.0GeV/pattern0/negative_250.0_330.0_pattern_0.npy')
    ##
    positive_211_259 = np.row_stack((positive_211_250,positive_250_330[np.where(positive_250_330[:,-1]<259)[0],:]))
    positive_259_330 = positive_250_330[np.where(positive_250_330[:,-1]>259)[0],:]
    negative_211_259 = np.row_stack((negative_211_250,negative_250_330[np.where(negative_250_330[:,-1]>-259)[0],:]))
    negative_259_330 = negative_250_330[np.where(negative_250_330[:,-1]<-259)[0],:]
    ##
    if os.path.isdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/211.0_259.0GeV/pattern0/'):
        pass
    else:
        mkdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/211.0_259.0GeV/pattern0/')
    if os.path.isdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/259.0_330.0GeV/pattern0/'):
        pass
    else:
        mkdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/259.0_330.0GeV/pattern0/')
    if os.path.isdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/211.0_259.0GeV/pattern0/'):
        pass
    else:
        mkdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/211.0_259.0GeV/pattern0/')
    if os.path.isdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/259.0_330.0GeV/pattern0/'):
        pass
    else:
        mkdir(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/259.0_330.0GeV/pattern0/')
    ##
    np.save(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/211.0_259.0GeV/pattern0/positive_211.0_259.0_pattern_0.npy',positive_211_259)
    np.save(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/positive/259.0_330.0GeV/pattern0/positive_259.0_330.0_pattern_0.npy',positive_259_330)
    np.save(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/211.0_259.0GeV/pattern0/negative_211.0_259.0_pattern_0.npy',negative_211_259)
    np.save(workpath+"/" + arguments.dataset + '/0.8_1130GeV_v2/results/rawdata/transferdata/negative/259.0_330.0GeV/pattern0/negative_259.0_330.0_pattern_0.npy',negative_259_330)    

##### rebinning_450_JUAMS #########################################################################



#########################################################################################
endtime = time.time()
print ((endtime - begintime)/60)



