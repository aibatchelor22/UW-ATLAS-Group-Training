#Routine for 2D multiple splits until Ecrit is reached for all partons, initial parton along "x" direction
import random
import math
import matplotlib.pyplot as plt
%matplotlib inline
import pylab as pl
from matplotlib import collections as mc

zmin = 0.
zmax = 1.
thetamin = 0.
thetamax = math.pi / 2.
Ei = 100.
m = 1.
Pxi = math.sqrt(Ei**2-m**2)
Pyi = 0
Pi = [Ei, Pxi, Pyi]
Ecrit = 2.
nextLoc = [0.,0.]
momentumVectors = []
moreSplits = True

nextSplitData = [[Ei,Pxi,Pyi,0.,0.]]

#Graph initial parton as line from [-1,0] to [0,0]
lines = [[(-1, -0), (0, 0)]]
lcinit = mc.LineCollection(lines, linewidths=2)
fig, ax = pl.subplots()
ax.add_collection(lcinit)


#iterative block for each location
while moreSplits == True:

  moreSplits = False #This will become true if there is a split parton with E > Ecrit
  currentSplitData = nextSplitData
  #clear the split locations and data for the next iteration
  nextSplitData = []
    
  for currentSplitLoc in currentSplitData:

    #Retrieve current vector Pi from currentSplitData and assign variables from its components
    Ei = currentSplitLoc[0]
    Pxi = currentSplitLoc[1]
    Pyi = currentSplitLoc[2]
    Pi = [Ei, Pxi, Pyi]

    #determine thetaInitial from components of the current vector Pi
    thetaInitial = math.tan(Pyi/Pxi)

    #split stochastically and use total angles to calculate the components

    #Calculate z stochastically
    def zCDF(pos):
      zcdf = math.log(1+pos)
      return zcdf;
    def zCDFINV(pos):
      zcdfinv = math.exp(pos) - 1
      return zcdfinv;
    yzmin = zCDF(zmin)
    yzmax = zCDF(zmax)
    yz = (yzmax-yzmin)*random.random() + yzmin
    z = zCDFINV(yz)

    #Calculate theta stochastically
    def thetaCDF(pos):
      thetacdf = math.log(1+pos)
      return thetacdf;
    def thetaCDFINV(pos):
      thetacdfinv = math.exp(pos) - 1
      return thetacdfinv;
    ythetamin = thetaCDF(thetamin)
    ythetamax = thetaCDF(thetamax)
    ytheta = (ythetamax-ythetamin)*random.random() + ythetamin
    theta = thetaCDFINV(ytheta)


    #calculate two new momentum vectors for the split partons
   

    E1 = z*Ei
    E2 = (1-z)*Ei
    
    if (E1 > m) and (E2 > m):
      #Calculate P1 vector relative to the path of travel of the mother particle
      Px1Prime = math.sqrt(E1**2-m**2)*math.cos(theta)
      Py1Prime = math.sqrt(E1**2-m**2)*math.sin(theta)
  
      #Rotation operation R(thetaInitial) to set P1 vector at angle thetaInital w.r.t. the lab frame
      Px1 = Px1Prime*math.cos(thetaInitial) - Py1Prime*math.sin(thetaInitial)
      Py1 = Px1Prime*math.sin(thetaInitial) + Py1Prime*math.cos(thetaInitial)
      P1 = [E1, Px1, Py1]

      #Calculate P2 vector based on the P1 vector (Pi=P1+P2) using conservation of energy and momentum
   
      Px2 = Pxi-Px1
      Py2 = Pyi-Py1
      P2 = [E2, Px2, Py2]


      #store vectors for the two partons in the grand list of momentum vectors (always!)
      momentumVectors.append(P1)
      momentumVectors.append(P2)

      #define locations of next parton splittings (unitary vectors oriented at resepctive theta)
      currentLoc = [currentSplitLoc[3],currentSplitLoc[4]]
      nextLoc1 = [currentLoc[0] + Px1, currentLoc[1] + Py1]
      nextLoc2 = [currentLoc[0] + Px2, currentLoc[1] + Py2]


    
      #plot lines to the next two splits or end points (always!)
      lines = lines + [[(currentLoc[0], currentLoc[1]), (nextLoc1[0], nextLoc1[1])], [(currentLoc[0], currentLoc[1]), (nextLoc2[0], nextLoc2[1])]]
      lc = mc.LineCollection(lines, linewidths=2)
      ax.add_collection(lc)
      ax.autoscale()
    

      if E1 > Ecrit:

        #store vector and location for next split in nextSplitData list
        nextSplitData.append([E1,Px1,Py1,nextLoc1[0],nextLoc1[1]])

        #Set this to True so that the while loop repeats.  If there are no more splits to do (after all location iterations), the shower ends.
        moreSplits = True

      if E2 > Ecrit:

        #store vector and location for next split in nextSplitData list
        nextSplitData.append([E2,Px2,Py2,nextLoc2[0],nextLoc2[1]])
        
        #Set this to True so that the while loop repeats. If there are no more splits to do (after all location iterations), the shower ends.
        moreSplits = True
