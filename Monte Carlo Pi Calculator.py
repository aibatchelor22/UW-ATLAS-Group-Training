import numpy as np
import random
import math
import matplotlib.pyplot as plt
%matplotlib inline


#Calculate many values of pi and plot a histogram
dist = []
for i in range(1,100):
  Naccepted = 0
  Ntotal = 100000
  for j in range(1,Ntotal):
    x = random.random()
    y = random.random()
    if (x**2 + y**2) <= 1:
      Naccepted = Naccepted + 1
  pivalue = 4*Naccepted/Ntotal
  dist.append(pivalue)
mean = np.mean(dist)
std = np.std(dist)
print(mean)
print(std)
plt.hist(dist, bins=20)