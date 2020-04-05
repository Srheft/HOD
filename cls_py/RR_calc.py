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


####### FAIZAN'S BINNING  from DD_rp-pi_QSO_NGC_v7_2.dat
edges = [0.089125,0.112202,0.141254,0.177828,0.223872,0.281838,0.354813,0.446684,0.562341,0.707946,0.891251,1.122018,1.412538,1.778279,2.238721,2.818383,3.548134,4.466836,5.623413,7.079458,8.912509,11.220185,14.125375,17.782794,22.387211,28.183829,35.481339,44.668359,56.234133,70.794578,89.125094,112.201845,141.253754,177.827941]
mids = [0.1,0.125893,0.158489,0.199526,0.251189,0.316228,0.398107,0.501187,0.630957,0.794328,1.0,1.258925,1.584893,1.995262,2.511886,3.162278,3.981072,5.011872,6.309573,7.943282,10,12.589254,15.848932,19.952623,25.118864,31.622777,39.810717,50.118723,63.095734,79.432823,100.0,125.892541,158.489319]
nbins = len(mids)

#################################################################
      
if __name__ == '__main__':
  
  path='/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/Mar2020/'

  filename = path+'eBOSS_QSO_clustering_NGC_v7_1_withS.ran.fits'
  a=fitsio.read(filename)
  ra = a['RA']  
  dec = a['DEC']
  ### w=w_SYSTOT*w_FKP*w_NOZ
  w = a['WEIGHT_SYSTOT']*a['WEIGHT_FKP']*a['WEIGHT_NOZ']

  # Python equivalent of matchlength
  radius_in_deg = 10.  # degrees
  
  size = len(ra)
  doJob = True
  p = 1000
  
  n_cpu =  int(sys.argv[1])
  batch_no = int(sys.argv[2])
  cpu_ID = int(sys.argv[3]) 
  
  n = batch_no * n_cpu + cpu_ID
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
        cols = fits.ColDefs([col1, col2, col3,col4,col5])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        # clobber = True  (it overwrites/updates)
        tbhdu.writeto('./RR_QSO_NGC_v7_1/spherematch.'+str(int(n))+'_RR.fits', clobber=True)
        
      #######################
      ind1 = m1 + low
      ind2 = m2
      s1 = a['s'][ind1]
      s2 = a['s'][ind2]
      d = distance
      rp = (s1+s2)*d*pi/180./2.
      pii = abs(s1-s2)
      weights = w[ind1]*w[ind2]
      arr_rp = np.zeros((10000,nbins))
      warr_rp = np.zeros((10000,nbins))
      
      print('start ...')

      for j in range(nbins):

        wrp  = ((rp < edges[j+1]) & (rp > edges[j]))
        arr = np.zeros(10000)
        warr = np.zeros(10000)

        if np.sum(wrp)>0:
            w = weights[wrp]
            a = pii[wrp]

            for k in range(len(a)):
                arr[int(a[k])] += 1
                warr[int(a[k])] += 1*w[k]

        arr_rp[:,j] = arr
        warr_rp[:,j] = warr
       

      sumarr_rp_RR = arr_rp
      sumarr_rp_wRR = warr_rp


      np.savetxt(path+"/RR_QSO_NGC_v7_1/wRR_2Darray_batch_"+str(int(n))+".txt", sumarr_rp_wRR)
      np.savetxt(path+"/RR_QSO_NGC_v7_1/RR_2Darray_batch_"+str(int(n))+".txt",sumarr_rp_RR)
     
      
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
  
  
    
