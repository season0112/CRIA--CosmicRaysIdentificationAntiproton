import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-D','--dataset', help='which_dataset_you_choose')
parser.add_argument('-E','--energy_range', choices=['16.6-38.9GeV', '38.9-147GeV', '38.9-147GeV_v2','147-1000GeV','147-1000GeV_v2'], help='which_energy_you_choose')
arguments = parser.parse_args()

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)

if arguments.energy_range == '147-1000GeV_v2':
    binnings = np.array([147, 175, 211, 250, 330, 525, 1000])
elif arguments.energy_range == '38.9-147GeV_v2':
    binnings = np.array([38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0])
else:
    print("no energy range is given")
    os._exit(0)    

a0 = np.array([])

for binleft in range(binnings.shape[0]-1):
    b_tem = np.load("/p/scratch/cvsk10/li8/analysis_v7.0/" + arguments.dataset + "/" + arguments.energy_range + '/results/rawdata/transferdata/' + 'negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0'+'/'+'negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy')
    if a0.shape[0] == 0:
        a0 = b_tem
    else:
        a0 = np.row_stack((a0,b_tem)) 

plt.figure(figsize=(9,9))
plt.hist(a0[:,34],bins=900,log=True)
plt.savefig('/p/scratch/cvsk10/li8/analysis_v7.0/plot_true_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'.png')

plt.figure(figsize=(9,9))
plt.hist(a0[:,43],bins=900, log=True,)
plt.savefig('/p/scratch/cvsk10/li8/analysis_v7.0/plot_measure_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'.png')





