#!/bin/bash

#### Define binning and Rigidity Range
bin=(0.8 1.0 1.16 1.33 1.51 1.71 1.92 2.15 2.4 2.67 2.97 3.29 3.64 4.02 4.43 4.88 5.37 5.9 6.47 7.09 7.76 8.48 9.26 10.1 11.0 12.0 13.0 14.1 15.3 16.6 18.0 19.5 21.1 22.8 24.7 26.7 28.8 31.1 33.5 36.1 38.9 41.9 45.1 48.5 52.2 56.1 60.3 64.8 69.7 74.9 80.5 93 108 125 147 175 211 259 450)
binname=(0.8 1.0 1.16 1.33 1.51 1.71 1.92 2.15 2.4 2.67 2.97 3.29 3.64 4.02 4.43 4.88 5.37 5.9 6.47 7.09 7.76 8.48 9.26 10.1 11.0 12.0 13.0 14.1 15.3 16.6 18 19.5 21.1 22.8 24.7 26.7 28.8 31.1 33.5 36.1 38.9 41.9 45.1 48.5 52.2 56.1 60.3 64.8 69.7 74.9 80.5 93 108 125 147 175 211 259 450) # 18.0->18 etc.
#zsh shell index is from 1. [1] is 0.8, [2] is 1.0,  [11] is 2.97, [24] is 10.1, [30] is 16.6, [31] is 18.0.
#bash shell index is from 0. [0] is 0.8, [1] is 1.0,  [10] is 2.97, [23] is 10.1, [29] is 16.6, [30] is 18.0.

if [ $2 = intermediaterange ]; then
    Treename=AntiprotonIntermediateEnergy
elif [ $2 = highrange ]; then
    Treename=ExampleAnalysis
elif [ $2 = lowrange ]; then
    Treename=AntiprotonLowEnergy_Tree
fi


#### AC_Merge loop for all Rigidity bins
startnumber=$(($3))
endnumber=$(($4))

for (( i=${startnumber}; i<${endnumber}; i=i+1 ))
do
    echo ${i}
    echo ${bin[${i}]}_${bin[$((${i}+1))]}
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_May2015_${bin[${i}]}_${bin[$((${i}+1))]}.root --resultprefix $5_Tree_$6_May2015_${binname[${i}]}_${binname[$((${i}+1))]} --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_Nov2017_${bin[${i}]}_${bin[$((${i}+1))]}.root --resultprefix $5_Tree_$6_Nov2017_${binname[${i}]}_${binname[$((${i}+1))]} --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_all_${bin[${i}]}_${bin[$((${i}+1))]}.root --resultprefix $5_Tree_$6_${binname[${i}]}_${binname[$((${i}+1))]} --resultdir $1
done


#### Make a short version of large data tree as a test file (For low rigidity range) 
#### (Originally has 1000 root files, for test only merge 50 root files, therefore need to rescale 20 times later)
if [[ $2 = lowrange ]] && [[ $6 = positive ]]; then
    for (( i=${startnumber}; i<${endnumber}; i=i+1 ))
    do
        ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_00*5[0-4]_May2015_${bin[${i}]}_${bin[$((${i}+1))]}.root --resultprefix $5_Tree_$6_May2015_${binname[${i}]}_${binname[$((${i}+1))]}_test --resultdir $1
        ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_00*5[0-4]_Nov2017_${bin[${i}]}_${bin[$((${i}+1))]}.root --resultprefix $5_Tree_$6_Nov2017_${binname[${i}]}_${binname[$((${i}+1))]}_test --resultdir $1
        ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_00*5[0-4]_all_${bin[${i}]}_${bin[$((${i}+1))]}.root --resultprefix $5_Tree_$6_${binname[${i}]}_${binname[$((${i}+1))]}_test --resultdir $1
    done
fi


#### Extra step for High Rigidity range in last few bins:
if [ $2 = highrange ]; then
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_May2015_211_250.root --resultprefix $5_Tree_$6_May2015_211_250 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_May2015_250_330.root --resultprefix $5_Tree_$6_May2015_250_330 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_May2015_259_330.root --resultprefix $5_Tree_$6_May2015_259_330 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_May2015_330_525.root --resultprefix $5_Tree_$6_May2015_330_525 --resultdir $1

    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_Nov2017_211_250.root --resultprefix $5_Tree_$6_Nov2017_211_250 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_Nov2017_250_330.root --resultprefix $5_Tree_$6_Nov2017_250_330 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_Nov2017_259_330.root --resultprefix $5_Tree_$6_Nov2017_259_330 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_Nov2017_330_525.root --resultprefix $5_Tree_$6_Nov2017_330_525 --resultdir $1

    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_all_211_250.root --resultprefix $5_Tree_$6_211_250 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_all_250_330.root --resultprefix $5_Tree_$6_250_330 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_all_259_330.root --resultprefix $5_Tree_$6_259_330 --resultdir $1
    ac_merge --input $1/$5/$6/results/finebinroot/${Treename}_*_all_330_525.root --resultprefix $5_Tree_$6_330_525 --resultdir $1
fi


