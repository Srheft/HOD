from numpy import *
import numpy as np
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
import scipy.optimize as op
from scipy import stats
import emcee
import sys
import corner
import matplotlib as plt

filename = 'FIT1_NOV2018-4KDE+KO12+eBOSS.dat'
#filename = 'FIT2_NOV2018-allKO12+eBOSS.dat'
#filename = 'FIT3_NOV2018-4KDE+eBOSS.dat'

dat = loadtxt(filename)

naccept = dat[:,0]

fsat = dat[:,1]

Mm = dat[:,2]
delm = dat[:,3]
chi2 = dat[:,4]

wgood = np.where(Mm <0.24e13) # FIT1 TRIM
#wgood= np.where(Mm < 0.22e13)  # FIT2 TRIM
#wgood= np.where(fsat > 0.045) # FIT3 TRIM

Mm = Mm[wgood]
fsat = fsat[wgood]
chi2 = chi2[wgood]
delm = delm[wgood]


w=np.where(chi2 == min(chi2))
#w=np.where(fsat <= 0.011 and fsat >= 0.009)

Mm= np.asarray(Mm)

samples2=np.zeros((len(fsat),2))
samples2[:,0]= fsat
#samples2[:,1]= delm
samples2[:,1]= Mm

fsat_true, mass_true = np.mean(samples2[:,0]),np.mean(samples2[:,1]) #FIT1 

#fsat_true, mass_true = stats.mode(samples2[:,0])[0][0],stats.mode(samples2[:,1])[0][0]

#fsat_true, mass_true =fsat[w][0], Mm[w][0]
#fsat_true, mass_true =0.010, 1.3e12 # FIT2 cross hair

#print len(w),len(fsat),min(fsat),max(chi2), max(chi2),len(np.where(chi2))
print('F_sat=',fsat_true, 'M_m=',mass_true/1e12,'x e12 M_sun/h')

fig = corner.corner(samples2, bins=15, labels=[r"$f_{sat}$", r"$M_{m}$"],truths=[fsat_true, mass_true],color='blue',fontsize=25,truth_color='red',scale_hist=False,space=0)

#fig.show()
fig.savefig("CornerFIT1_NOV2018_KDE+KO12_eBOSS.eps")
