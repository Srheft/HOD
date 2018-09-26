# -*- coding: utf-8 -*-
from time import time
from time import gmtime, strftime
import sys
import os
import numpy as np
from math import *
from astrometry.util.starutil_numpy import * 
from astrometry.libkd import spherematch
import random 
import pyfits
from pyfits import *
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
  
  
  filename = 'CORE_87300_2.20000_3.39979_WITHsyscuts_RADECzs.fits'
  a = pyfits.open(filename)
  hdr = getheader(filename, 1) # header
  d = a[1].data
  
  ra = d.field('ra')
  dec = d.field('dec')


  #indices = np.arange(len(ra))
  ##print indices
  #indices = my_shuffle(indices)
  ##print indices
  #ra = ra[indices]
  #dec = dec[indices]
  
 
  # Python equivalent of matchlength
  radius_in_deg = 25.  # degrees
  r = angleto3Dradius(radius_in_deg)
  
  size = len(ra)
  doJob = True
  p = 100
  
  #n_cpu =  int(sys.argv[1])
  #batch_no = int(sys.argv[2])
  #cpu_ID = int(sys.argv[3])
      
  #n = batch_no * n_cpu + cpu_ID
  nn=ceil(len(ra)/p)
  
  for k in range(int(nn)):


   low = p*k
   up = p*(k+1)
   if low <= size and up > size:
      up = size
   if low >= size:
      doJob = False
   if up<low or up<0 or low<0:
      doJob = False
   print low, up
  
   if doJob:

      ra_small = ra[low:up]
      dec_small = dec[low:up]
      xyz_small = radectoxyz(ra_small, dec_small)
      xyz = radectoxyz(ra, dec)

    
      #print "No.: ", n, "Go ... "
      
  
      time0 = time()
      
      #############

      (inds,dists) = spherematch.match(xyz_small, xyz, r)
      dist_in_deg = dist2deg(dists)
      m1, m2 = inds[:,0], inds[:,1]
      d = dist_in_deg[:,0]
      w=np.where(d > 0.)
      print 'Out of ',len(d),' pairs ',len(d[w]), ' of them are real pairs'


      ####################
      ## For k3match
      ## Uncomment k3match-import section at the top of this file 
      ## All tested for peorformance and accuracy by E.Kourkchi
      ## (August 03 2015)
      #####################
      ##(m1, m2, d) =  celestial(ra_small, dec_small, ra, dec, radius_in_deg)



      #############
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
	col1 = fits.Column(name='index1', format='K', array=m1_2[w])
	col2 = fits.Column(name='index2', format='K', array=m2[w])
	col3 = fits.Column(name='dist', format='D', array=d[w])
	cols = fits.ColDefs([col1, col2, col3])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	
	# clobber = True  (it overwrites/updates)
	tbhdu.writeto('pairs/python.out.'+str(int(k))+'_DD.100.fits', clobber=True)
      
      #############     
      writeHist = True
      if writeHist:
	bin_delta = 0.25
	bins = np.arange(0, radius_in_deg + bin_delta, bin_delta)
	hist = np.histogram(d[w], bins = bins)
	col1 = fits.Column(name='Hist_TOT', format='J', array=hist[0])
	col2 = col1   # If we don't calcualte NGC separately

	

	### Calculating North
	NGC = True
	if NGC:
	  # NGC: ra > 90 and ra < 270
	  
	  size_pair = len(d)
	  
	  ra_small_pair = ra_small[m1]
	  #dec_small_pair = dec_small[m1]
	  ra_pair = ra[m2]
	  #dec_pair = dec[m2]
	  
	  q1 = np.zeros((size_pair,), dtype=np.int)
	  q2 = np.zeros((size_pair,), dtype=np.int)
	  q3 = np.zeros((size_pair,), dtype=np.int)
	  q4 = np.zeros((size_pair,), dtype=np.int)
	  q1[np.where(ra_small_pair>90)] = 1
	  q2[np.where(ra_small_pair<270)] = 1
	  q3[np.where(ra_pair>90)] = 1
	  q4[np.where(ra_pair<270)] = 1
	  qq = q1+q2+q3+q4
	  indices = np.where(qq==4)
	  
	  d_NGC = d[indices]
	  
	  w=where(d_NGC > 0)
	  hist = np.histogram(d_NGC[w], bins = bins)
	  col2 = fits.Column(name='Hist_NGC', format='J', array=hist[0])

          indices_sgc = np.where(qq!=4)
	  d_SGC = d[indices_sgc]
	  w=where(d_SGC > 0)
	  hist = np.histogram(d_SGC[w], bins = bins)
	  col3 = fits.Column(name='Hist_SGC', format='J', array=hist[0])
	### END - Calculating North


	cols = fits.ColDefs([col1, col2, col3])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto('DD_withsyscus/python.out.'+str(int(k))+'.DD_100.NGC_SGC_hist.fits', clobber=True)
      ### END - Write Histogram     
      
      print
      print "run_time = ", (time()-time0)/60., 'sec'
      print "# of total matches:", len(d)
      print "Catalog size:", size
      print "lower index: ", low
      print "upper index: ", up
      print
  
      #### test
      #ra_max = ra[m1[where(d == max(d))]]
      #dec_max = dec[m1[where(d == max(d))]]
      #print "Max Distance:", max(d), ra_max, dec_max
  
  if not doJob:
    print
    print "Catalog size:", size
    print "lower index: ", low
    print "upper index: ", up
    print "No Job was done ..."
    print "[Error] Index out of bound"
    print 
  print "Finished at:", strftime("%Y-%m-%d %H:%M:%S", gmtime())
  print "*********************************"
  print
  print
  
  
    
