import numpy as np
import random
import math
import matplotlib.pyplot as plt
%matplotlib inline

#Inverse transform distribution with falling exponential probability density
dist = []
for i in range(1,10000):
  y = (1-math.exp(-5))*random.random()
  x = math.log(1/(1-y))
  dist.append(x)

#Plot histogram with matplotlib
plt.hist(dist, bins=50)