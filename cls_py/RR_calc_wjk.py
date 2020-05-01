# -*- coding: utf-8 -*-
from time import time
from time import gmtime, strftime
import sys
import os
import numpy as np
from math import *
from astrometry.util.starutil_numpy import * 
from astrometry.libkd import spherematch
import random, fitsio
from numpy import radians, degrees, sin, cos, arctan2, hypot
from astropy.io import fits

#python data_data_srh.py 1 0 0 for a single cpu 


## If you are going to use 
## k3match, uncomment this line 
#from k3match import *
#################################################################


####### BINNING  from DD_rp-pi_QSO_NGC_v7_2.dat
edges = [0.089125,0.112202,0.141254,0.177828,0.223872,0.281838,0.354813,0.446684,0.562341,0.707946,0.891251,1.122018,1.412538,1.778279,2.238721,2.818383,3.548134,4.466836,5.623413,7.079458,8.912509,11.220185,14.125375,17.782794,22.387211,28.183829,35.481339,44.668359,56.234133,70.794578,89.125094,112.201845,141.253754,177.827941]
mids = [0.1,0.125893,0.158489,0.199526,0.251189,0.316228,0.398107,0.501187,0.630957,0.794328,1.0,1.258925,1.584893,1.995262,2.511886,3.162278,3.981072,5.011872,6.309573,7.943282,10,12.589254,15.848932,19.952623,25.118864,31.622777,39.810717,50.118723,63.095734,79.432823,100.0,125.892541,158.489319]
nbins = len(mids)

#################################################################
      
if __name__ == '__main__':
  
  path='/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/'

  filename = path+'/Mar2020/eBOSS_QSO_NGC_pip_v7_2.ran_z0.8_z2.2_withS_withPIX.fits'
  a=fitsio.read(filename)
  ra = a['RA']  
  dec = a['DEC']
  pix = a['PIXEL']
  ### w=w_SYSTOT*w_FKP*w_NOZ
  w = a['WEIGHT_SYSTOT']*a['WEIGHT_FKP']*a['WEIGHT_NOZ']

  # Python equivalent of matchlength
  radius_in_deg = 10.  # degrees
  
  size = len(ra)
  doJob = True
  p = 1000
  
  #n_cpu =  int(sys.argv[1])
  batch_no = int(sys.argv[1])
  #cpu_ID = int(sys.argv[3]) 
  
  n = batch_no #* n_cpu + cpu_ID
  low = p*n
  up = p*(n+1)
  if low <= size and up > size:
      up = size
  if low >= size:
      doJob = False
  if up<low or up<0 or low<0:
      doJob = False
  
  
  if doJob:

      ra_small = ra[low:up]
      dec_small = dec[low:up]
      w_small = w[low:up]
      pix_small = pix[low:up]
      print("No.: ", n, "Go ... ")
      
      time0 = time()
      
      (m1,m2,distance) = spherematch.match_radec(ra_small, dec_small, ra, dec, radius_in_deg,notself=True)

      writeFits = False
      if writeFits:
	
        m1_2 = m1 + low
	# Creating the output fits table
	
	# Format J for 32 bit Integer
	# Format D for 64 bit floating point 
	# L: Logical (Boolean)
	# B: Unsigned Byte
	# I: 16-bit Integer
	# J: 32-bit Integer
	# K: 64-bit Integer
	# E: Single-precision Floating Point
	# D: Double-precision Floating Point
	# C: Single-precision Complex
	# M: Double-precision Complex
	# A: Character
        
        #col1 = fits.Column(name='ra1', format='D', array=ra[m1_2])
        #col2 = fits.Column(name='dec1', format='D', array=dec[m1_2])
        col1 = fits.Column(name='index1', format='K', array=m1_2)
        #col4 = fits.Column(name='ra2', format='D', array=ra[m2])
        #col5 = fits.Column(name='dec2', format='D', array=dec[m2])
        col2 = fits.Column(name='index2', format='K', array=m2)
        col3 = fits.Column(name='dist', format='D', array=distance)
        col4 = fits.Column(name='w1', format='D', array=w[m1_2])
        col5 = fits.Column(name='w2', format='D', array=w[m2])
        col6 = fits.Column(name='pix1',format='K',array=pix[m1_2])
        col7 = fits.Column(name='pix2',format='K',array=pix[m2]) 
        col8 = fits.Column(name='s1',format='D',array=a['s'][m1_2])
        col9 = fits.Column(name='s2',format='D',array=a['s'][m2])
        col10 = fits.Column(name='w1w2',format='D',array=w[m2]*w[m1_2])
        cols = fits.ColDefs([col1, col2, col3,col4,col5,col6,col7,col8,col9,col10])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        # clobber = True  (it overwrites/updates)
        tbhdu.writeto(path+'/RR_QSO_NGC_v7_2/spherematch.'+str(int(n))+'_RR.fits', clobber=True)
        
      #######################
      ind1 = m1 + low
      ind2 = m2
      s1 = a['s'][ind1]
      s2 = a['s'][ind2]
      d = distance
      rp = (s1+s2)*d*pi/180./2.
      pii = abs(s1-s2)

      weights = w[ind1]*w[ind2]
      
      arr_rp_pi = np.zeros((10000,nbins))
      
      rp_bin_weights = np.zeros(nbins)
               
      print('start binning ...')

      for j in range(nbins):

        wrp  = ((rp < edges[j+1]) & (rp > edges[j]))
        arr = np.zeros(10000)
        
        if np.sum(wrp)>0:

            rp_bin_weights[j] = np.sum(weights[wrp]) #total weights of the pairs in this bin of rp
            a = pii[wrp]
            arr,b = np.histogram(a,bins=np.arange(10001))
             
        arr_rp_pi[:,j] = arr
             
      
      ####### JACKKNIF:  getting RR pairs in pixels:
      print([ind1])
      pixx1 = pix[ind1]
      pixx2 = pix[ind2]
      inpix = (pixx1 == pixx2)

      #import healpy as hp
      #nside = 16 
      #NPIX = hp.nside2npix(nside)
      npix = 3072  # nside =16 -> total number of pixels NPIX = 3072

      RR_jk = np.zeros((npix,nbins))

      pixes = np.unique(pixx1)

      weights_RR_pix = np.zeros((npix,nbins))

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
                  arr_rp_pi_jk[p,:,j] = arr  # [:,j] is the j-th column of the xi   
      #np.savetxt(path+"/RR_QSO_NGC_v7_2/wRR_2Darray_batch_"+str(int(n))+".txt", sumarr_rp_wRR)
      #np.savetxt(path+"/RR_QSO_NGC_v7_2/RR_2Darray_batch_"+str(int(n))+".txt",sumarr_rp_RR)
      np.savez_compressed(path+"/Apr2020/RR_QSO_NGC_v7_2/wRR_2Darray_batch_"+str(int(n))+".npz", WEIGHTS =rp_bin_weights,XI_RP_PII=arr_rp_pi,PIX_WEIGHTS=rp_bin_weights_jk ,PIX_XI_RP_PI=arr_rp_pi_jk) 

      #grid = np.load(path+'test.npz')

      #grid['X'] = sumarr_rp_wRR

      print('')
      print("run_time = ", (time()-time0)/60., 'min')
      print("# of total matches:", len(distance))
      print("Catalog size:", size)
      print("lower index: ", low)
      print("upper index: ", up)
      print('')
  
  
  if not doJob:
    print
    print("Catalog size:", size)
    print("lower index: ", low)
    print("upper index: ", up)
    print("No Job was done ...")
    print("[Error] Index out of bound")
    print('') 
  print("Finished at:", strftime("%Y-%m-%d %H:%M:%S", gmtime()))
  print("*********************************")
  print('')
  print('')
  
  
    
