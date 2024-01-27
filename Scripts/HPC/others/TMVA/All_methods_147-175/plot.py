import numpy as np
import matplotlib.pyplot as plt
# load data
Fisher_signal = np.load('mva_Fisher_signal_array.npy')
Fisher_background = np.load('mva_Fisher_background_array.npy')

Likelihood_signal = np.load('mva_Likelihood_signal_array.npy')
Likelihood_background = np.load('mva_Likelihood_background_array.npy')

LikelihoodKDE_signal = np.load('mva_LikelihoodKDE_signal_array.npy')
LikelihoodKDE_background = np.load('mva_LikelihoodKDE_background_array.npy')

SVM_signal = np.load('mva_SVM_signal_array.npy')
SVM_background = np.load('mva_SVM_background_array.npy')

MLPBNN_signal = np.load('mva_MLPBNN_signal_array.npy')
MLPBNN_background = np.load('mva_MLPBNN_background_array.npy')

# rerange
Fisher_signal[np.where(Fisher_signal<-1)[0]] = -1
Fisher_background[np.where(Fisher_background<-1)[0]] = -1
Fisher_signal[np.where(Fisher_signal>1)[0]] = 1
Fisher_background[np.where(Fisher_background>1)[0]] = 1

Likelihood_signal[np.where(Likelihood_signal<-1)[0]] = -1
Likelihood_background[np.where(Likelihood_background<-1)[0]] = -1
Likelihood_signal[np.where(Likelihood_signal>1)[0]] = 1
Likelihood_background[np.where(Likelihood_background>1)[0]] = 1


plt.figure(figsize=(18,18))
plt.hist(Fisher_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Fisher_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('Fisher.png')

plt.figure(figsize=(18,18))
plt.hist(Fisher_signal, bins=80, alpha=0.5,log=True, label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Fisher_background, bins=80, alpha=0.5,log=True, label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('Fisher_log.png')

plt.figure(figsize=(18,18))
plt.hist(Likelihood_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Likelihood_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('Likelihood.png')

plt.figure(figsize=(18,18))
plt.hist(Likelihood_signal, bins=80, alpha=0.5,log=True,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(Likelihood_background, bins=80, alpha=0.5,log=True,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('Likelihood_log.png')

plt.figure(figsize=(18,18))
plt.hist(LikelihoodKDE_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(LikelihoodKDE_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('LikelihoodKDE.png')

plt.figure(figsize=(18,18))
plt.hist(LikelihoodKDE_signal, bins=80, alpha=0.5,log=True,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(LikelihoodKDE_background, bins=80, alpha=0.5,log=True,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('LikelihoodKDE_log.png')

plt.figure(figsize=(18,18))
plt.hist(SVM_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(SVM_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('SVM.png')

plt.figure(figsize=(18,18))
plt.hist(SVM_signal, bins=80, alpha=0.5,log=True, label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(SVM_background, bins=80, alpha=0.5,log=True, label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('SVM_log.png')

plt.figure(figsize=(18,18))
plt.hist(MLPBNN_signal, bins=80, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(MLPBNN_background, bins=80, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('MLPBNN.png')

plt.figure(figsize=(18,18))
plt.hist(MLPBNN_signal, bins=80, alpha=0.5,log=True,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(MLPBNN_background, bins=80, alpha=0.5,log=True,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.legend(fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.savefig('MLPBNN_log.png')



