import astropy.io.fits as fits
import numpy as np
import fitsio, os, sys
from astropy.table import Column
from astropy.table import Table
from math import cos,pi,log10
import matplotlib.pylab as plt
import healpy as hp
import pylab as py
from desitarget.geomask import radec_match_to 

path='/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/Mar2020/'
plotpath='/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/Mar2020/plots/'

#### MY BINNING: Binningrp: 30kpc/h - ~200 Mpc/h
#powmax=2.5
#powmin=-1.5
#delpow=0.19
#edges=10.**(np.arange(powmin,powmax,delpow))
#mids=[]
#for i in range(len(edges)-1):
#        mids.append((10**((powmin)+((i+1)*(delpow)))+10**((powmin)+(i*(delpow))))/2)
#nbins = len(mids)-1
#print(nbins,mids, edges,len(mids))

##############################################
####### FAIZAN'S BINNING  from DD_rp-pi_QSO_NGC_v7_2.dat
edges = [0.089125,0.112202,0.141254,0.177828,0.223872,0.281838,0.354813,0.446684,0.562341,0.707946,0.891251,1.122018,1.412538,1.778279,2.238721,2.818383,3.548134,4.466836,5.623413,7.079458,8.912509,11.220185,14.125375,17.782794,22.387211,28.183829,35.481339,44.668359,56.234133,70.794578,89.125094,112.201845,141.253754,177.827941]
mids = [0.1,0.125893,0.158489,0.199526,0.251189,0.316228,0.398107,0.501187,0.630957,0.794328,1.0,1.258925,1.584893,1.995262,2.511886,3.162278,3.981072,5.011872,6.309573,7.943282,10,12.589254,15.848932,19.952623,25.118864,31.622777,39.810717,50.118723,63.095734,79.432823,100.0,125.892541,158.489319]
nbins = len(mids)
##############################################
import glob
ddfiles = glob.glob(path+'/DD_QSO_NGC_v7_2/*_DD.fits') 

qso,h = fitsio.read(path+'eBOSS_QSO_clustering_NGC_v7_2_withS_with60PIPs.fits',header=True)


sumarr_rp_DD = np.zeros((10000,nbins))
sumarr_rp_wDD = np.zeros((10000,nbins))

for f in ddfiles:
    print('working on ',f)
    dd,h = fitsio.read(f,header=True)
    
    ind1 = dd['index1']
    ind2 = dd['index2']
    w1 = dd['w1']
    w2 = dd['w2']
    weights = w1*w2

    pips = []
    for p in range(len(dd)):
        c = (qso['PIPs'][ind1[p]] & qso['PIPs'][ind2[p]])
        cbit = [bin(ci) for ci in c]
        nones = sum([cb.count('1') for cb in cbit])
        if nones >0:
            pips.append(1860/nones)
        else:
            pips.append(0)
        
 
    pips = np.asarray(pips)
    
    pipnweights = pips * weights
    
    s1 =  qso[ind1]['s']
    s2 = qso[ind2]['s']
    
    d = dd['dist']
    rp = (s1+s2)*d*pi/180./2.
    pii = abs(s1-s2)
    arr_rp = np.zeros((10000,nbins))
    warr_rp = np.zeros((10000,nbins))
    
    
    print('start ...')
    
    for j in range(nbins):

        wrp  = ((rp < edges[j+1]) & (rp > edges[j]))
        arr = np.zeros(10000)
        warr = np.zeros(10000)

        if np.sum(wrp)>0:
            w = pipnweights[wrp]
            a = pii[wrp]
            
            for k in range(len(a)):
                arr[int(a[k])] += 1
                warr[int(a[k])] += 1*w[k]

        arr_rp[:,j] = arr   
        warr_rp[:,j] = warr    
        #print('bin', j,'DD', np.sum(arr))
        #print('bin', j,'wDD', np.sum(warr))

    sumarr_rp_DD = sumarr_rp_DD + arr_rp   
    sumarr_rp_wDD = sumarr_rp_wDD + warr_rp   
    
wDDp = []
DDp = []
for i in range(nbins):

    
    DDp.append(sum(sumarr_rp_DD[:,i]))
    wDDp.append(sum(sumarr_rp_wDD[:,i]))
    print(edges[i],'\t',edges[i+1],'\t',mids[i],'\t','DD: ',sum(sumarr_rp_DD[:,i]),'\t','norm_DD: ',sum(sumarr_rp_wDD[:,i]))
    

np.savetxt("normDD_2Darray.txt", sumarr_rp_wDD)
np.savetxt("DD_2Darray.txt",sumarr_rp_DD)  

print('*******************************************************************') 

