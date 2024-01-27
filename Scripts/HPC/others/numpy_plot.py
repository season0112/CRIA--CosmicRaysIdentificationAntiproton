#!/usr/bin/env python
from root_numpy import root2array, tree2array
from ROOT import TFile
import numpy as np
import os
import argparse
import math
import json
import collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import shutil


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)

parser = argparse.ArgumentParser()
parser.add_argument('--inputnumpyfile', help='inputnumpyfile')
parser.add_argument('--path', help='current path')
parser.add_argument('--variable', help='which variable you choose')
parser.add_argument('--low', help='low range')
parser.add_argument('--high', help='high range')
arguments = parser.parse_args()

os.chdir(arguments.path)
mkdir("plots")

inputfile = np.load(arguments.inputnumpyfile)

rigidity = inputfile['Rigidity']
tofbeta = inputfile['TofBeta']
mass2 = rigidity**2 * (1-tofbeta**2) / tofbeta**2

####
plt.figure(figsize=(18,18))
plt.hist(mass2, 500,  alpha=0.5, range=(-6,6), facecolor='blue',edgecolor='black')
plt.xlabel("mass2",fontsize=30)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.legend(loc='upper right',fontsize=30)
plt.savefig(arguments.inputnumpyfile+'mass2.png')

plt.figure(figsize=(18,18))
plt.hist(inputfile[arguments.variable], 500,  alpha=0.5, range=(float(arguments.low),float(arguments.high)), facecolor='blue',edgecolor='black')
plt.xlabel(arguments.variable,fontsize=30)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.legend(loc='upper right',fontsize=30)
plt.savefig(arguments.inputnumpyfile + arguments.variable+'.png')

shutil.move(arguments.inputnumpyfile+'mass2.png', 'plots/'+arguments.inputnumpyfile[0:-4]+'_'+'mass2.png')
shutil.move(arguments.inputnumpyfile+arguments.variable+'.png', 'plots/'+arguments.inputnumpyfile[0:-4]+'_'+arguments.variable+'.png')


