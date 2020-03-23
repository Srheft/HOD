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


#################################################################
def my_shuffle(array):
        random.seed(0)
        random.shuffle(array)
        return array
#################################################################

def angleto3Dradius(angle, isDegree=True):
  
  if isDegree:
    angle = angle*pi/180.
  
  
  return sqrt((sin(angle))**2 + (1-cos(angle))**2)

#################################################################

      
if __name__ == '__main__':
  
  path='/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/Mar2020/'

  filename = path+'eBOSS_QSO_clustering_NGC_v7_2_ids.dat.fits'
  a=fitsio.read(filename)
  ra = a['RA']  
  dec = a['DEC']
  ### w=w_SYSTOT*w_FKP*w_NOZ
  w = a['WEIGHT_SYSTOT']*a['WEIGHT_FKP']*a['WEIGHT_NOZ']

  indices = np.arange(len(ra))
 
  indices = my_shuffle(indices)
 
  ra = ra[indices]
  dec = dec[indices]
  w = w[indices]
 
  # Python equivalent of matchlength
  radius_in_deg = 5.  # degrees
  
  size = len(ra)
  doJob = True
  p = 10000
  
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


      writeFits = True
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
        tbhdu.writeto('./DD_QSO_NGC_v7_2/spherematch.'+str(int(n))+'_DD.fits', clobber=True)
      
      #######################     
      
      print('')
      print("run_time = ", (time()-time0)/60., 'sec')
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
  
  
    
