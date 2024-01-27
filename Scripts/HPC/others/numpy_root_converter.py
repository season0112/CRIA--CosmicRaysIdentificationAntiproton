#!/usr/bin/env python
from root_numpy import root2array, tree2array
from ROOT import TFile
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--mode', help='which transfer mode you want')
parser.add_argument('--inputrootfile', help='inputrootfile')
parser.add_argument('--treename', help='treename')
parser.add_argument('--inputnumpyfile', help='inputnumpyfile')
parser.add_argument('--path', help='current path')
arguments = parser.parse_args()


os.chdir(arguments.path)

if arguments.mode == "root2array":
    array = root2array(str(arguments.inputrootfile), arguments.treename)
    np.save(arguments.inputrootfile[0:-5]+".npy" , array)

elif arguments.mode == "tree2array":
    rootfile = TFile(arguments.inputrootfile)
    roottree = rootfile.Get(arguments.treename)
    array = tree2array(roottree)
    np.save(arguments.inputrootfile[0:-5]+".npy" , array) 

elif arguments.mode == "array2tree":
    array = np.load(arguments.inputnumpyfile)
##    tree = array2tree(array, name="tree")
    array2root(array, arguments.inputnumpyfile[0:-4]+'.root', 'tree', 'recreate')


