import numpy as np
import random
import math
import matplotlib.pyplot as plt
%matplotlib inline

#Inverse transform distribution with gaussian probability density on range [0,10]
from scipy.special import erf
from scipy.special import erfinv
mu = 5
sigma = 2
x1 = 0
x2 = 10
y1 = 0.5*(1+erf((x1-mu)/(sigma*math.sqrt(2))))
y2 = 0.5*(1+erf((x2-mu)/(sigma*math.sqrt(2))))
dist = []
for i in range(1,100000):
  y = (y2-y1)*random.random()+y1
  x = sigma*math.sqrt(2)*erfinv(2*y-1) + mu
  dist.append(x)

#Plot histogram with matplotlib
plt.hist(dist, bins=50)