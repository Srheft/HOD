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
import gmpy
### to run on a single cpu:
#python DD_Angular_paircount.py 1 0 0 


#################################################################

      
if __name__ == '__main__':
  
  path='/utahpath/July2020/'

  filename = 'eBOSS_LRG_NGC_pip_v7_2_new_withS_z0.6_z1.0.fits'
  a = fitsio.read(filename)
  ra = a['RA']  
  dec = a['DEC']
  #pix = a['PIXEL']
  s = a['s']
  ### for DD and DR pairs, only three of the weights are going to be use[see PIP paper]: w=w_SYSTOT*w_FKP*w_NOZ
  w = a['WEIGHT_SYSTOT']*a['WEIGHT_FKP']*a['WEIGHT_NOZ']

  indices = np.arange(len(ra))
 
  ra = ra[indices]
  dec = dec[indices]
  w = w[indices]
 
  # Python equivalent of matchlength
  radius_in_deg = 7.0  # degrees
  
  size = len(ra)
  doJob = True
  p = 10000
  
  batch_no = int(sys.argv[1])

  n = batch_no 
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

      writeFits = True
      if writeFits:
	
        m1_2 = m1 + low


        PIP_fib = np.ones(len(m1))*1860.

        ll = (distance < 1.0) ####  the angular scale over which PIP and ANGULAR upweghts are effective is way smaller than 1 degrees at the redshift range of eboss samples
        PIP_fib[ll] = 0         

        num = np.sum(ll)
        
        ind = np.where(ll)[0]
        print(num)
        for i in ind:

            for j in range(60):

                 wei = (a[m1_2[i]]['WEIGHT_BW'][j]) & (a[m2[i]]['WEIGHT_BW'][j])

                 PIP_fib[i] += gmpy.popcount(int(wei))
                 #print(PIP_fib[i])
        print('Done with PIP weight calculation')

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
        
        col1 = fits.Column(name='index1', format='K', array=m1_2)
        col2 = fits.Column(name='index2', format='K', array=m2)
        col3 = fits.Column(name='dist', format='D', array=distance)
        col4 = fits.Column(name='w1', format='D', array=w[m1_2])
        col5 = fits.Column(name='w2', format='D', array=w[m2])
        col6 = fits.Column(name='s1',format='D',array=s[m1_2])
        col7 = fits.Column(name='s2',format='D',array=s[m2])
        col8 = fits.Column(name='w1w2',format='D',array=w[m1_2]*w[m2])
        col9 = fits.Column(name='PIP',format='D',array=PIP_fib)
        cols = fits.ColDefs([col1, col2, col3,col4,col5,col6,col7,col8,col9])

        tbhdu = fits.BinTableHDU.from_columns(cols)
	
	# clobber = True  (it overwrites/updates)
        tbhdu.writeto(path+'/DD_NGC_LRG_z0.6_z1.0/spherematch.'+str(int(n))+'_DD.fits', clobber=True)
      
      #######################     
      
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
  
  
    
