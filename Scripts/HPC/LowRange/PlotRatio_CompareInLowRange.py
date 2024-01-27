import asireader
import numpy as np
import matplotlib.pyplot as plt

# parameters for nice axis labels
plt.rcParams['font.size'] = 24.0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'sans\-serif'
plt.rcParams['mathtext.it'] = 'sans\-serif:italic'

# Load Published Data
reader = asireader.ASIReader_R()
readerPAMELA = asireader.ASIReader_R_PAMELAData()
AMSPaper2016Ratio = reader.readfluxratio("/home/bo791269/Software/antiprotonanalysis/ReferenceFiles/AntiprotonToProtonRatio/pbar_p_AMS_PRL2016_000.xml")
PAMELA2013Ratio = readerPAMELA.readfluxratio("/home/bo791269/Software/antiprotonanalysis/ReferenceFiles/AntiprotonToProtonRatio/pbar_p_PAM_JETPlett2013_R_000.xml") # Remove the first (0.5-1.01GV) and last bin (180-350GV) due to data information are different there.
AMSPaper2016Ratio.R = asireader.LaffertyWyatt(AMSPaper2016Ratio.rigidity_min, AMSPaper2016Ratio.rigidity_max, gamma=2.7)
PAMELA2013Ratio.R = asireader.LaffertyWyatt(PAMELA2013Ratio.rigidity_min, PAMELA2013Ratio.rigidity_max, gamma=2.7)

# Load My result
#f_pass78 = TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v10.0/" + "/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultspass7.8.root")
#Ratio_pass78 = f_pass78.Get("gRatio_unfolded")

#### Plot 
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16.0, 10.0))
# Draw 
plt.errorbar( AMSPaper2016Ratio.R, AMSPaper2016Ratio.fluxratio*10000, yerr=AMSPaper2016Ratio.err_total*10000, marker='o', linestyle="None", markersize=6, color='red', label=r' AMS2016Paper')
plt.errorbar( PAMELA2013Ratio.R, PAMELA2013Ratio.fluxratio*10000, yerr=(PAMELA2013Ratio.err_total_low*10000,PAMELA2013Ratio.err_total_high*10000), marker='o', linestyle="None", markersize=6, color='blue', label=r' PAMELA2013')

plt.xlabel("Rigidity (GV)",fontsize=40)
plt.ylabel("Antiproton to Proton Ratio",fontsize=40)
plt.legend(loc='best',fontsize=30)
#axes.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
#axes.tick_params(axis='y', labeltop=False, labelbottom=True,labelleft=False,labelright=True)
axes.tick_params(axis='both', which='both', direction='in',length=10, width=3)
plt.xticks(fontsize=30)
plt.yticks(fontsize=40)
axes.set_xlim(0.6,5)
axes.set_ylim(0,1.5)
#plt.xscale('log')
#plt.yscale('log')
plt.savefig("LowRatioCompare.pdf")



