from numpy import *
import numpy as np
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
import scipy.optimize as op
from scipy import stats
import emcee
import sys
import corner
import matplotlib as plt

#filename = 'chain2_2pars_new_ONLY1KDE_ALL.dat'
filename = 'chain2_2pars_newdatapoints_ALL.dat'
dat = loadtxt(filename)

naccept = dat[:,0]

fsat = dat[:,1]

Mm = dat[:,2]
delm = dat[:,3]
chi2 = dat[:,4]

w=np.where(chi2 == min(chi2))

Mm= np.asarray(Mm)

samples2=np.zeros((len(fsat),2))
samples2[:,0]= fsat
#samples2[:,1]= delm
samples2[:,1]= Mm

#fsat_true, mass_true = np.mean(samples2[:,0]),np.mean(samples2[:,1])

#fsat_true, mass_true = stats.mode(samples2[:,0])[0][0],stats.mode(samples2[:,1])[0][0]

fsat_true, mass_true =fsat[w][0], Mm[w][0]

#print len(w),len(fsat),min(fsat),max(chi2), max(chi2),len(np.where(chi2))
print 'F_sat=',fsat_true, 'M_m=',mass_true/1e12,'x e12 M_sun/h'

fig = corner.corner(samples2, bins=15, labels=[r"$f_{sat}$", r"$M_{m}$"],truths=[fsat_true, mass_true],color='blue',fontsize=25,truth_color='red',scale_hist=False,space=0)

#fig.show()
fig.savefig("cornerplot_SEP18_4KDE.eps")
