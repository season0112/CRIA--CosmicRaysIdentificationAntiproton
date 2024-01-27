import numpy as np
import matplotlib.pyplot as plt

array_all = []

for i in range(126):
    array  = np.load("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v8.0/totalall/numpyfiles/negative/test/negative_2.67_2.97_" + str(i) + ".npy")
    array_all.append(array.shape[0])
    
array_all = np.array(array_all)
print(array_all.shape)
print(array_all[0])
print(array_all[50])




plt.figure(figsize=(18,18))
plt.plot(np.arange(126), array_all, '*',lw=8)
plt.yscale('log')
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.savefig('test.png')



