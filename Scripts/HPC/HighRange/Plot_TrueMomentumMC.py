import numpy as np
import os
import matplotlib.pyplot as plt
import binning

resultpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis'
highpath = os.getenv('HPCHIGHENERGYDATADIR')

B1042all = np.array([])
for binname in binning.bins_525_zhili:
    a0 = np.load(highpath + "/ISS_anylsis/data/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.1800_7.6_all_Tree_negative_" + binname + "_Pattern_0_VGG16NN.npy")
    if B1042all.shape[0] == 0 :
        B1042all = a0
    else:
        B1042all = np.row_stack((B1042all,a0))


print(B1042all[:,3].shape)

fig, ax = plt.subplots(figsize=(30,18))
plt.hist(B1042all[:,3], 50, density=False, facecolor='g', alpha=0.75)
ax.tick_params(direction='in',length=10, width=3)
plt.xticks(fontsize=60)
plt.yticks(fontsize=80)
plt.ylabel("",fontsize=50)
#plt.xscale('log')
##plt.yscale('log')
plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + "R_True" + ".pdf")




