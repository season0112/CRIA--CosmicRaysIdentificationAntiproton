import numpy as np
import matplotlib.pyplot as plt
# load data
Fisher_signal = np.load('application_results/mva_Fisher_signal_array.npy')
Fisher_background = np.load('application_results/mva_Fisher_background_array.npy')

Likelihood_signal = np.load('application_results/mva_Likelihood_signal_array.npy')
Likelihood_background = np.load('application_results/mva_Likelihood_background_array.npy')

LikelihoodKDE_signal = np.load('application_results/mva_LikelihoodKDE_signal_array.npy')
LikelihoodKDE_background = np.load('application_results/mva_LikelihoodKDE_background_array.npy')

SVM_signal = np.load('application_results/mva_SVM_signal_array.npy')
SVM_background = np.load('application_results/mva_SVM_background_array.npy')

MLPBNN_signal = np.load('application_results/mva_MLPBNN_signal_array.npy')
MLPBNN_background = np.load('application_results/mva_MLPBNN_background_array.npy')
'''
# rerange
Fisher_signal[np.where(Fisher_signal<-1)[0]] = -1
Fisher_background[np.where(Fisher_background<-1)[0]] = -1
Fisher_signal[np.where(Fisher_signal>1)[0]] = 1
Fisher_background[np.where(Fisher_background>1)[0]] = 1

Likelihood_signal[np.where(Likelihood_signal<-1)[0]] = -1
Likelihood_background[np.where(Likelihood_background<-1)[0]] = -1
Likelihood_signal[np.where(Likelihood_signal>1)[0]] = 1
Likelihood_background[np.where(Likelihood_background>1)[0]] = 1
'''
plt.figure(figsize=(18,18))
plt.hist(Fisher_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Fisher_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/Fisher.png')

plt.figure(figsize=(18,18))
plt.hist(Fisher_signal, bins=80, alpha=0.5,log=True, label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Fisher_background, bins=80, alpha=0.5,log=True, label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/Fisher_log.png')

plt.figure(figsize=(18,18))
plt.hist(Likelihood_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Likelihood_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/Likelihood.png')

plt.figure(figsize=(18,18))
plt.hist(Likelihood_signal, bins=80, alpha=0.5,log=True,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Likelihood_background, bins=80, alpha=0.5,log=True,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/Likelihood_log.png')

plt.figure(figsize=(18,18))
plt.hist(LikelihoodKDE_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(LikelihoodKDE_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/LikelihoodKDE.png')

plt.figure(figsize=(18,18))
plt.hist(LikelihoodKDE_signal, bins=80, alpha=0.5,log=True,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(LikelihoodKDE_background, bins=80, alpha=0.5,log=True,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/LikelihoodKDE_log.png')

plt.figure(figsize=(18,18))
plt.hist(SVM_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(SVM_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/SVM.png')

plt.figure(figsize=(18,18))
plt.hist(SVM_signal, bins=80, alpha=0.5,log=True, label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(SVM_background, bins=80, alpha=0.5,log=True, label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/SVM_log.png')

plt.figure(figsize=(18,18))
plt.hist(MLPBNN_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(MLPBNN_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/MLPBNN.png')

plt.figure(figsize=(18,18))
plt.hist(MLPBNN_signal, bins=80, alpha=0.5,log=True,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(MLPBNN_background, bins=80, alpha=0.5,log=True,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=35)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('plot/MLPBNN_log.png')


