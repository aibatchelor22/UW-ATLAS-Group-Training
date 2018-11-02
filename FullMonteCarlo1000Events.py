#Routine for 3D multiple splits until Ecrit is reached for all partons, initial parton along "z" direction, followed by clustering algorithm
import random
import math
import matplotlib.pyplot as plt
%matplotlib inline
import pylab as pl
import numpy as np
from matplotlib import collections as mc

#Simplified histogram method for number of jets in 1000 events
jetHistData = []
jets = []
for a in range(1,1001):

  #First parton
  zmin = 0.
  zmax = 1.
  thetamin = 0.
  thetamax = math.pi / 2.
  phimin = 0.
  phimax = math.pi
  Ei = 100.
  minit = 2.
  m = 1
  thetai = thetai = random.random()*math.pi/2
  phii = random.random()*math.pi

  Pmagi = math.sqrt(Ei**2-minit**2)
  Ecrit = 2.
  momentumVectors = []
  moreSplits = True

  nextSplitData = [[Ei,Pmagi,thetai,phii]]

  #iterative block for each location
  while moreSplits == True:

    moreSplits = False #This will become true if there is a split parton with E > Ecrit
    currentSplitData = nextSplitData
    #clear the split locations and data for the next iteration
    nextSplitData = []
    
    for currentSplitLoc in currentSplitData:

      #Retrieve current vector Pi from currentSplitData and assign variables from its components
      Ei = currentSplitLoc[0]
      Pmagi = currentSplitLoc[1]
      thetainit = currentSplitLoc[2]
      phiinit = currentSplitLoc[3]
      Pi = [Ei, Pmagi, thetai, phiinit]

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
        Pmag1 = math.sqrt(E1**2-m**2)

        #Rotate angular components to align with path of mother particle
        theta1 = thetainit + theta
        phi1 = (phiinit + phi) % (2 * math.pi)
        if theta1 > math.pi:
          theta1 = math.pi - theta1
          phi1 = -1*phi1

      
        P1 = [E1,Pmag1, theta1, phi1]
      

        #Calculate P2 vector based on the P1 vector (Pi=P1+P2) using conservation of energy and momentum
   
        Pmag2 = Pmagi-Pmag1
        theta2 = math.atan(Pmag1*math.sin(theta)/(Pmagi-Pmag1))
        phi2 = (phiinit - phi1) % (2 * math.pi)
        if theta2 > math.pi:
          theta2 = math.pi - theta2
          phi2 = -1*phi2
        P2 = [E2, Pmag2, theta2, phi2]


        #store vectors for the two partons in the grand list of momentum vectors (always!)
        momentumVectors.append(P1)
        momentumVectors.append(P2)

         

        if E1 > Ecrit:

          #store vector and location for next split in nextSplitData list
          nextSplitData.append(P1)

          #Set this to True so that the while loop repeats.  If there are no more splits to do (after all location iterations), the shower ends.
          moreSplits = True

        if E2 > Ecrit:

          #store vector for next split in nextSplitData list
          nextSplitData.append(P2)
        
          #Set this to True so that the while loop repeats. If there are no more splits to do (after all location iterations), the shower ends.
          moreSplits = True



  #Second parton
  zmin = 0.
  zmax = 1.
  thetamin = 0.
  thetamax = math.pi / 2.
  phimin = 0.
  phimax = math.pi
  Ei = 100.
  minit = 2.
  m = 1
  thetai = math.pi - thetai
  phii = math.pi + phii
  Pmagi = math.sqrt(Ei**2-minit**2)
  Ecrit = 2.

  moreSplits = True

  nextSplitData = [[Ei,Pmagi,thetai,phii]]

  #iterative block for each location
  while moreSplits == True:

    moreSplits = False #This will become true if there is a split parton with E > Ecrit
    currentSplitData = nextSplitData
    #clear the split locations and data for the next iteration
    nextSplitData = []
    
    for currentSplitLoc in currentSplitData:

      #Retrieve current vector Pi from currentSplitData and assign variables from its components
      Ei = currentSplitLoc[0]
      Pmagi = currentSplitLoc[1]
      thetai = currentSplitLoc[2]
      phii = currentSplitLoc[3]
      Pi = [Ei, Pmagi, thetai, phii]

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
        Pmag1 = math.sqrt(E1**2-m**2)

        #Rotate angular components to align with path of mother particle
        theta1 = thetainit + theta
        phi1 = (phiinit + phi) % (2 * math.pi)
        if theta1 > math.pi:
          theta1 = math.pi - theta1
          phi1 = -1*phi1 #Fix this for negative theta
        P1 = [E1,Pmag1, theta1, phi1]
      

        #Calculate P2 vector based on the P1 vector (Pi=P1+P2) using conservation of energy and momentum
   
        Pmag2 = Pmagi-Pmag1
        theta2 = math.atan(Pmag1*math.sin(theta)/(Pmagi-Pmag1))
        phi2 = (phiinit - phi1) % (2 * math.pi)
        if theta2 > math.pi:
          theta2 = math.pi - theta2
          phi2 = -1*phi2
        P2 = [E2, Pmag2, theta2, phi2]


        #store vectors for the two partons in the grand list of momentum vectors (always!)
        momentumVectors.append(P1)
        momentumVectors.append(P2)

         

        if E1 > Ecrit:

          #store vector and location for next split in nextSplitData list
          nextSplitData.append(P1)

          #Set this to True so that the while loop repeats.  If there are no more splits to do (after all location iterations), the shower ends.
          moreSplits = True

        if E2 > Ecrit:

          #store vector for next split in nextSplitData list
          nextSplitData.append(P2)
        
          #Set this to True so that the while loop repeats. If there are no more splits to do (after all location iterations), the shower ends.
          moreSplits = True

  import pandas as pd

  n = 0
  R = 0.5  #Usually between 0.4 and 0.7, see http://iopscience.iop.org/article/10.1088/1742-6596/645/1/012008/pdf

  pseudojets = []
  jets = []
  for vector in momentumVectors:
    E = vector[0]
    Pmag = vector[1]
    theta = vector[2]
    PT = Pmag*math.sin(theta)
    pseudorapidity = -1*math.log(math.tan(abs(theta/2)))
    phi = vector[3]
    pseudojets.append([E,PT,phi,pseudorapidity,])



  #iterate over all pseudojets until pseudojets list is empty
  while pseudojets != []:
    diBlist = []
    dijlist = []
    pseudojetscount = 0
    #Populate beam distance data diBData
    for pseudojetinstance in pseudojets:
      diB=(pseudojetinstance[1])**2
      diBlist.append([diB,pseudojetinstance[0],pseudojetinstance[1],pseudojetinstance[2],pseudojetinstance[3]])
      pseudojetscount = pseudojetscount + 1
    diBdata = pd.DataFrame(data=diBlist, columns = ['diB', 'E', 'PT', 'phi', 'pseudorapidity'])
    #Populate dijDdata with iteration for vectors (avoids duplicates) to 
    for i in range(1, pseudojetscount):
      for j in range(0, i):
        if (i != j):
          P1 = pseudojets[i]
          P2 = pseudojets[j]
          PT1 = P1[1]
          PT2 = P2[1]
          phi1 = P1[2]
          phi2 = P2[2]
          pseudorapidity1 = P1[3]
          pseudorapidity2 = P2[3]
  
          Delta12 = math.sqrt(((phi1-phi2) ** 2) + ((pseudorapidity1-pseudorapidity2) ** 2))
          dij = (min((PT1 ** (2*n)), (PT2 ** (2*n))) * Delta12) / R

          #Populate dataframe of dij (also include two four vectors associated with each dij and the corresponding indices)
          dijlist.append([dij,P1[0],P1[1],P1[2],P1[3],P2[0],P2[1],P2[2],P2[3]])       
    dijdata = pd.DataFrame(data=dijlist, columns = ['dij', 'E1', 'PT1', 'phi1', 'pseudorapidity1','E2', 'PT2', 'phi2', 'pseudorapidity2'])
    #Find minimum of diB and dij and label it dmin

  
    dmin1 = dijdata['dij'].min()
    dmin2 = diBdata['diB'].min()
    dmin = min(dmin1,dmin2)

    #If dmin is one of dij, create new pseudojet with four vector P = P1 + P2 (see page 5, arXiv:hep-ph/9305266)
    if dmin == dmin1:

      #find corresponding values for P1 and P2 in dataframe dijData df.loc[df['B'] == 3, 'A'].item()
      E1 = dijdata.loc[dijdata['dij']==dmin, 'E1'].item()
      E2 = dijdata.loc[dijdata['dij']==dmin, 'E2'].item()
      PT1 = dijdata.loc[dijdata['dij']==dmin, 'PT1'].item()
      PT2 = dijdata.loc[dijdata['dij']==dmin, 'PT2'].item()
      phi1  = dijdata.loc[dijdata['dij']==dmin, 'phi1'].item()
      phi2 = dijdata.loc[dijdata['dij']==dmin, 'phi2'].item()
      pseudorapidity1 = dijdata.loc[dijdata['dij']==dmin, 'pseudorapidity1'].item()
      pseudorapidity2 = dijdata.loc[dijdata['dij']==dmin, 'pseudorapidity2'].item()

      #combine pseudojets and put into the pseudojet list
      P = [E1+E2, PT1+PT2, ((phi1+phi2)/2), ((pseudorapidity1+pseudorapidity2)/2)]
      pseudojets.append(P)

      #remove these entries from list of pseudojets
      P1 = [E1, PT1, phi1, pseudorapidity1]
      P2 = [E2, PT2, phi2, pseudorapidity2]
      pseudojets.remove(P1)
      pseudojets.remove(P2)

    #If dmin is diB, then the pseudojet four vector associated with diB is stored as a jet and removed from the list of pseudojets
    elif dmin == dmin2:
      #find corresponding value in dataframe, retrieve momentum data, add to jets
      E = diBdata.loc[diBdata['diB']==dmin, 'E'].item()
      PT = diBdata.loc[diBdata['diB']==dmin, 'PT'].item()
      phi = diBdata.loc[diBdata['diB']==dmin, 'phi'].item()
      pseudorapidity = diBdata.loc[diBdata['diB']==dmin, 'pseudorapidity'].item()
      P = [E, PT, phi, pseudorapidity]
      jets.append(P)
      #remove this entry from list of pseudojets
      pseudojets.remove(P)
    else:
      for pseudojetleft in pseudojets:
        jets.append(pseudojetleft)
      pseudojets=[]
      
  jetDataFrame = pd.DataFrame(data=jets, columns = ['E', 'PT', 'phi', 'pseudorapidity'])
  jetHistData.append(jetDataFrame['E'].count())
plt.hist(jetHistData, bins = 50)