#Routine for 3D multiple splits until Ecrit is reached for all partons, initial parton along "z" direction
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
phimin = 0.
phimax = math.pi
Ei = 100.
m = 1.
Pxi = 0
Pyi = 0
Pzi = math.sqrt(Ei**2-m**2)
Ecrit = 2.
nextLoc = [0.,0.,0.,0.]
momentumVectors = []
moreSplits = True

nextSplitData = [[Ei,Pxi,Pyi,Pzi]]

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
    Pzi = currentSplitLoc[3]
    Pi = [Ei, Pxi, Pyi, Pzi]

    #determine initial angular components of mother particle Cartesian components of the current vector Pi
    phiInitial = math.taninv(math.sqrt(Pxi**2+Pyi**2)/Pzi)
    thetaInitial = math.taninv(Pyi/Pxi)

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

    #Calculate phi stochastically
    phi = (phimax-phimin)*random.random() + phimin


    #calculate two new momentum vectors for the split partons
   

    E1 = z*Ei
    E2 = (1-z)*Ei
    
    if (E1 > m) and (E2 > m):

      #Calculate P1Prime magnitude relative to the path of travel of the mother particle
      P1Mag = math.sqrt(E1**2-m**2)   

      #Rotate angular components to align with path of mother particle
      phiTotal = phiInitial + phi
      thetaTotal = thetaInitial + theta
  
      #Convert from spherical to cartesian coordinates
      Px1 = P1Mag*math.sin(phiTotal)*math.cos(thetaTotal)
      Py1 = P1Mag*math.sin(phiTotal)*math.sine(thetaTotal)
      Pz1 = P1Mag*math.cos(phiTotal)
      P1 = [E1, Px1, Py1, Pz1]

      #Calculate P2 vector based on the P1 vector (Pi=P1+P2) using conservation of energy and momentum
   
      Px2 = Pxi-Px1
      Py2 = Pyi-Py1
      Pz2 = Pzi-Pz1
      P2 = [E2, Px2, Py2, Pz2]


      #store vectors for the two partons in the grand list of momentum vectors (always!)
      momentumVectors.append(P1)
      momentumVectors.append(P2)

      #define locations of next parton splittings
      nextLoc1 = [currentSplitLoc[1] + Px1, currentSplitLoc[2] + Py1, currentSplitLoc[3] + Pz2]
      nextLoc2 = [currentSplitLoc[1] + Px2, currentSplitLoc[2] + Py2, currentSplitLoc[3] + Pz2]

    

      if E1 > Ecrit:

        #store vector and location for next split in nextSplitData list
        nextSplitData.append([E1,Px1,Py1,Pz2])

        #Set this to True so that the while loop repeats.  If there are no more splits to do (after all location iterations), the shower ends.
        moreSplits = True

      if E2 > Ecrit:

        #store vector for next split in nextSplitData list
        nextSplitData.append([E2,Px2,Py2,Pz2])
        
        #Set this to True so that the while loop repeats. If there are no more splits to do (after all location iterations), the shower ends.
        moreSplits = True