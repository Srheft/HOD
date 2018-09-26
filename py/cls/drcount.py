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
  
  
  datfile = 'CORE_87300_2.20000_3.39979_WITHsyscuts_RADECzs.fits'
  ranfile = 'rancat_1623434_2.2000052_3.3998463_WITHsyscuts_RADECzs.fits'
  dat = pyfits.open(datfile)
  ran = pyfits.open(ranfile)
  ####hdr = getheader(filename, 1) # header
  ddat = dat[1].data
  dran = ran[1].data
  
  ra_dat = ddat.field('ra')
  dec_dat = ddat.field('dec')
  s_dat = ddat.field('s')

  ra_ran = dran.field('ra')
  dec_ran = dran.field('dec')
  s_ran = dran.field('s')

  #indices = np.arange(len(ra))
  ##print indices
  ##indices = my_shuffle(indices)
  ##print indices
  #ra_dat = ra_dat[indices]
  #dec_dat = dec_dat[indices]
  #ra_ran= ra_ran[indices]
  #dec_ran= dec_ran[indices]
  
 
  # Python equivalent of matchlength
  radius_in_deg = 25.  # degrees
  r = angleto3Dradius(radius_in_deg)
  
  size = len(ra_dat)
  doJob = True
  p =100
  
  #n_cpu =  int(sys.argv[1])
  #batch_no = int(sys.argv[2])
  #cpu_ID = int(sys.argv[3])
  
  #n = batch_no * n_cpu + cpu_ID
  
  nn=ceil(len(ra_dat)/p)
  
  for k in range(17,int(nn)):
   
   low = p*k
   up = p*(k+1)
   print k,low, up
   if low <= size and up > size:
      up = size
   if low >= size:
      doJob = False
   if up<low or up<0 or low<0:
      doJob = False
  
  
   if doJob:

      ra_small_dat = ra_dat[low:up]
      dec_small_dat = dec_dat[low:up]
      xyz_small_dat = radectoxyz(ra_small_dat, dec_small_dat)
      xyz_ran = radectoxyz(ra_ran, dec_ran)

    
      #print "No.: ", n, "Go ... "
      #print "patch No.: ", nn, "Go ... "
      
  
      time0 = time()
      
      #############

      (inds,dists) = spherematch.match(xyz_small_dat, xyz_ran, r)
      dist_in_deg = dist2deg(dists)
      m1, m2 = inds[:,0], inds[:,1]
      d = dist_in_deg[:,0]
      

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
	
	
	s1=s_dat[m1_2]
	s2=s_ran[m2]
	#print shape(s1),shape(s2),shape(d),shape(s1+s2)
	pi=abs(s1-s2)
	rp=(s1+s2)*d*pi/180./2.
	ss=sqrt(pi**2.+rp**2.)
	col1 = fits.Column(name='index1', format='K', array=m1_2)
	col2 = fits.Column(name='index2', format='K', array=m2)
	col3 = fits.Column(name='dist', format='D', array=d)
	col4 = fits.Column(name='pi', format='D', array=pi)
	col5 = fits.Column(name='rp', format='D', array=rp)
	col6 = fits.Column(name='s', format='D', array=ss)

        cols = fits.ColDefs([col1, col2, col3, col4, col5, col6])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	
	# clobber = True  (it overwrites/updates)
	tbhdu.writeto('DR_withsyscus/python.out.'+str(int(k))+'.100_d_pi_rp_s.fits', clobber=True)
        print 
      #############     
      writeHist = True
      if writeHist:
        m1_2 = m1 + low
        s1=s_dat[m1_2]
	s2=s_ran[m2]
	pi=abs(s1-s2)
	rp=(s1+s2)*d*pi/180./2.
	ss=sqrt(pi**2.+rp**2.)

	bin_delta = 0.25
	powmax=3.1
	powmin=-0.22
	delpow=0.082732
	bins = np.arange(0, radius_in_deg + bin_delta, bin_delta)
	midtheta=np.arange(bin_delta/2, radius_in_deg + bin_delta, bin_delta)
	    
	edges=10.**(np.arange(powmin,powmax,delpow))
	mids=[]
	for i in range(len(edges)):
            mids.append((10**((powmin)+((i+1)*(delpow)))+10**((powmin)+(i*(delpow))))/2)
        
        midbins=[]
        midbins[0:int(len(mids))-2]=mids[0:int(len(mids))-1]
        midbins[int(len(mids))-1:int(len(bins))-2]=np.zeros(len(bins)-len(mids))
	hist = np.histogram(d, bins = bins)
	hist_s = np.histogram(ss, bins = edges)
	col1 = fits.Column(name='midbins_deg', format='J', array=midtheta[0:100])
	col2 = fits.Column(name='Hist_TOT', format='J', array=hist[0])
	col3 = fits.Column(name='midbins_s', format='J', array=midbins)
        col4 = fits.Column(name='Hist_TOT_s', format='J', array=hist_s[0])
	

	### Calculating North
	NGC = True
	if NGC:
	  # NGC: ra > 90 and ra < 270
	  
	  size_pair = len(d)
	  
	  ra_small_pair = ra_small_dat[m1]
	  #dec_small_pair = dec_small[m1]
	  ra_pair = ra_ran[m2]
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
	  s_NGC = ss[indices]
	  w=np.where(d_NGC > 0.)
	  hist = np.histogram(d_NGC[w], bins = bins)
	  hist_s_ngc = np.histogram(s_NGC[w], bins = edges)
	  col5 = fits.Column(name='Hist_NGC', format='J', array=hist[0])
	  col6 = fits.Column(name='Hist_NGC_s', format='J', array=hist_s_ngc[0])
	### END - Calculating North

	### Calculating South
          indices = np.where(qq!=4)
          d_SGC = d[indices]
	  s_SGC = ss[indices]
	  w=np.where(d_SGC > 0.)
	  hist = np.histogram(d_SGC[w], bins = bins)
          hist_s_sgc = np.histogram(s_SGC[w], bins = edges)
          col7 = fits.Column(name='Hist_SGC', format='J', array=hist[0])
          col8 = fits.Column(name='Hist_SGC_s', format='J', array=hist_s_sgc[0])

	### END - Calculating North

	cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
	tbhdu = fits.BinTableHDU.from_columns(cols)
	tbhdu.writeto('DR_withsyscus_hist/python.out.'+str(int(k))+'.100.hist_d_s_NGC_SGC.fits', clobber=True)
      ### END - Write Histogram     
      
      #print
      #print "run_time = ", (time()-time0)/60., 'sec'
      #print "# of total matches:", len(d)
      #print "Random Catalog size:", len(ra_ran)
      #print "Data Catalog size:", len(ra_dat)
      #print
  
      #### test
      #ra_max = ra[m1[where(d == max(d))]]
      #dec_max = dec[m1[where(d == max(d))]]
      #print "Max Distance:", max(d), ra_max, dec_max
  
  if not doJob:
    print
    print "Random Catalog size:", len(ra_ran)
    print "lower index: ", low
    print "upper index: ", up
    print "No Job was done ..."
    print "[Error] Index out of bound"
    print 
  print "Finished at:", strftime("%Y-%m-%d %H:%M:%S", gmtime())
  print "*********************************"
  print
  print
  
  
    
