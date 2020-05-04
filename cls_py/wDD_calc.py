import astropy.io.fits as fits
import numpy as np
from time import time
import fitsio, os, sys
from astropy.table import Column
from astropy.table import Table
from math import cos,pi,log10
import matplotlib.pylab as plt
import healpy as hp
import pylab as py
from desitarget.geomask import radec_match_to 

path='/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/'


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
####### BINNING  from DD_rp-pi_QSO_NGC_v7_2.dat
edges = [0.089125,0.112202,0.141254,0.177828,0.223872,0.281838,0.354813,0.446684,0.562341,0.707946,0.891251,1.122018,1.412538,1.778279,2.238721,2.818383,3.548134,4.466836,5.623413,7.079458,8.912509,11.220185,14.125375,17.782794,22.387211,28.183829,35.481339,44.668359,56.234133,70.794578,89.125094,112.201845,141.253754,177.827941]
mids = [0.1,0.125893,0.158489,0.199526,0.251189,0.316228,0.398107,0.501187,0.630957,0.794328,1.0,1.258925,1.584893,1.995262,2.511886,3.162278,3.981072,5.011872,6.309573,7.943282,10,12.589254,15.848932,19.952623,25.118864,31.622777,39.810717,50.118723,63.095734,79.432823,100.0,125.892541,158.489319]
nbins = len(mids)
##############################################
import glob

ddfiles = glob.glob(path+'/Apr2020/DD_QSO_NGC_v7_2/*_DD.fits') 

qso,h = fitsio.read(path+'/Mar2020/eBOSS_QSO_NGC_pip_v7_2.dat_z0.8_z2.2_withS_withPIX.fits',header=True)

sumarr_rp_DD = np.zeros((10000,nbins))

rp_bin_angup = np.zeros(nbins)

rp_bin_weights = np.zeros(nbins)


time0 = time()

for f in ddfiles[0:1]:
    
    print('working on ',f)

    dd,h = fitsio.read(f,header=True)
    
    ind1 = dd['index1']
    ind2 = dd['index2']
    w1 = dd['w1']
    w2 = dd['w2']
    pix1 = dd['pix1']
    pix2 = dd['pix2']
    stat1 = qso['STATUS'][ind1]
    stat2 = qso['STATUS'][ind2]
    weights = w1*w2
    pips = []
    angup = []

    print('settint up took {} min ...'.format((time()-time0)/60., 'min' ))
    time0 = time()
    for p in range(len(dd)):

        c = (qso['WEIGHT_BW'][ind1[p]] & qso['WEIGHT_BW'][ind2[p]])
        cbit = [bin(ci) for ci in c]
        nones = sum([cb.count('1') for cb in cbit])
        if nones >0:
            pips.append(1860/nones)
        else:
            pips.append(0)
       
    print('Getting PIPs took {} min ...'.format((time()-time0)/60., 'min' ))
    time0 = time()

    pips = np.asarray(pips)
    
    pipnweights = pips * weights
    
    s1 = dd['s1']
    s2 = dd['s2']
    d = dd['dist']
    rp = (s1+s2)*d*pi/180./2.
    pii = abs(s1-s2)

    arr_rp = np.zeros((10000,nbins))
    
    #----------------------------------  intermediate ile for calculation verification
    write_pair_info = False
    
    if write_pair_info:
    
        t = Table(dd)

        col = Column(qso['RA'][ind1],name='RA1')
        t.add_column(col,index=0)

        col = Column(qso['DEC'][ind1],name='DEC1')
        t.add_column(col,index=0)

        col = Column(qso['RA'][ind2],name='RA2')
        t.add_column(col,index=0)

        col = Column(qso['DEC'][ind2],name='DEC2')
        t.add_column(col,index=0)

        col = Column(qso['EBOSS_TARGET_ID'][ind1],name='eBOSS_TARGET_ID1')
        t.add_column(col,index=0)

        col = Column(qso['EBOSS_TARGET_ID'][ind2],name='eBOSS_TARGET_ID2')
        t.add_column(col,index=0)

        col = Column(weights,name='SYSTOTxNOZxFKP_1x2')
        t.add_column(col,index=0)

        col = Column(pips,name='PIPs=1860/sumbitwise_bit1bit2')
        t.add_column(col,index=0)

        t.write(path+'DD_batch0_infotable.fits')
    #---------------------------------------------------------------------------------
    
    print('starting the rp-pi binning...')
    
    for j in range(nbins):
  
        time0 = time()
        wrp  = ((rp < edges[j+1]) & (rp > edges[j]))
        dd_pr = np.sum(wrp)
        w_fb = (((2**0 & stat2[wrp]) !=0) & ((stat2[wrp] & 2**0) !=0))
        dd_fb = np.sum(w_fb)

        angsep = d[wrp]
        #sep =  np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)+np.sin(dec1)*np.sin(dec2)
        #theta = np.arccos(sep)*180./np.pi
        #tmin = np.min(theta)
        #st  =  np.max(theta)-tmin
        #bint = int(log10(theta/tmin)/st)
        print('bin',j,edges[j],'<rp [Mpc/h] <', edges[j+1],np.round(np.min(angsep),3),'< theta [deg] <' ,np.round(np.max(angsep),3), 'Angup:', dd_pr/dd_fb)#, np,min(theta), np.max(theta)),bint)

        arr = np.zeros(10000)
  
        if np.sum(wrp)>0 :
                       
            rp_bin_weights[j] = np.sum(pipnweights[wrp]) #total weights of the pairs in this bin of rp
            a = pii[wrp]
            arr,b = np.histogram(a,bins=np.arange(10001))
            rp_bin_angup[j] = dd_pr*1.0/dd_fb 

        arr_rp[:,j] = arr   
        
        print('Got xi_rp_pii took {} min ...'.format((time()-time0)/60., 'min' ))

    sumarr_rp_DD = sumarr_rp_DD + arr_rp     
    
    JK = False
    #######  JACKKNIF:  getting RR pairs in pixels ##############################
    if JK:
        print([ind1])
        pixx1 = pix[ind1]
        pixx2 = pix[ind2]
        inpix = (pixx1 == pixx2)
        #import healpy as hp
        #nside = 16 
        #NPIX = hp.nside2npix(nside)
        npix = 3072  # nside =16 -> total number of pixels NPIX = 3072

        pixes = np.unique(pixx1)

        weights_DD_pix = np.zeros((npix,nbins))

        for p in list(pixes):

          ind = (pixx1 == p)
          pind = pixx2[ind]
          rp_pix = rp[pind]
          pii_pix = pii[pind]
          weights_pix = weights[pind]
          inp = (pind == p)
          rp_jk = rp_pix[inp]
          pii_jk = pii_pix[inp]

          weights_jk = weights_pix[inp]

          arr_rp_pi_jk = np.zeros((npix,10000,nbins))

          rp_bin_weights_jk = np.zeros((npix,nbins))

          for j in range(nbins):

              wrp  = ((rp_jk < edges[j+1]) & (rp_jk > edges[j]))

              rp_bin_weights_jk[p,j] = np.sum(weights_jk[wrp])

              arr = np.zeros(10000)
              if np.sum(wrp)>0:
                  a = pii_jk[wrp]
                  arr,b = np.histogram(a,bins=np.arange(10001))
                  arr_rp_pi_jk[p,:,j] = arr  

        np.savez_compressed(path+"/Apr2020/DD_JK_batch_"+str(int(n))+".npz", WEIGHTS =rp_bin_weights,XI_RP_PII=arr_rp_pi,PIX_WEIGHTS=rp_bin_weights_jk ,PIX_XI_RP_PI=arr_rp_pi_jk)

DDp = []
wDD_angup = []

for i in range(nbins):

    DDp.append(sum(sumarr_rp_DD[:,i]))
   
    wDD_angup.append(sum(sumarr_rp_DD[:,i])*rp_bin_angup[i]*rp_bin_weights[i])

    print(edges[i],'\t',edges[i+1],'\t',mids[i],'\t','DD: ',sum(sumarr_rp_DD[:,i]),'\t','norm_DD[PIP+NOZ,SYS,FKP]: ',sum(sumarr_rp_DD[:,i])*rp_bin_weights[i], '\t','norm_DD[all+angular]: ',sum(sumarr_rp_DD[:,i])*rp_bin_angup[i]*rp_bin_weights[i])
    
np.savez_compressed(path+'DD_QSO_NGC_v7_2_rp_pi_2darray.npz',DD_RP_PII=sumarr_rp_DD, PIPnWEIGHTS=rp_bin_weights, PIPnWEIGHTSnANGUP=rp_bin_angup)

print('*******************************************************************') 

