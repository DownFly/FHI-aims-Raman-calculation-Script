#!/usr/bin/python

#-------2014.02.06: download from http://www.kuliklab.org--------------------------
#-------2014.02.06: shanghui edit to plot IR spectoscopy for fhi-aism--------------
#how to use: python IR_plot.py 10 data  (here 10 is mean 10cm-1 tp broaden

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import sys
import os
from scipy.stats import norm

#---------------read data-------------------
broaden=8 #cm^-1 
rawdata=open(sys.argv[1],'r').readlines()    

mymin=100000.0
mymax=-100000.0
iy=0.0
peakcent=[]
peakintens=[]


for lines in range(0,len(rawdata)):
   freq=float(rawdata[lines].split()[0])
    
   if freq < mymin:
     mymin=freq
   if freq > mymax:
     mymax=freq

   peakcent.append(freq)
   peakintens.append(float(rawdata[lines].split()[1].strip('\n')))


if mymin < 0.0:
   mymin=0.0

points=mymax-mymin+200.0

#-----determines the x axis of the spectrum, and number of grid points
ix = np.linspace(mymin-100.0,mymax+100.0,int(points))

#-----peak center positions
for peaks in range(0,len(peakcent)):
#    iy+=2.51225*broaden*peakintens[peaks]*mlab.normpdf(ix,peakcent[peaks],broaden)
    iy+=2.51225*broaden*peakintens[peaks]*norm.pdf(ix,peakcent[peaks],broaden)


finwrite=open('Raman.dat','w')
for lines in range(0,len(ix)):
    finwrite.write('%s %s\n' %(ix[lines],0.004796*iy[lines]))

#-----plot------------------------
fig,ax = plt.subplots()
plt.xlabel("Raman shift (cm$^{-1}$)")
plt.ylabel('Intensity')
plt.xlim(200.0, 600)
#plt.ylim(-0.05, 2.0)
#ax.plot(ix,iy,color='blue',label="harmonic",lw=3)
#ax.plot(ix,iy,color='red',label=(sys.argv[2]+' harmonic Raman'),lw=3)
ax.plot(ix,iy,color='red',label=(sys.argv[2]),lw=3)
ax.legend(loc=9)
fig = plt.gcf()
#plt.show()
#fig.savefig(sys.argv[2] + '.pdf')
fig.savefig(sys.argv[2] + '.png')
