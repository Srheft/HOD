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
#import pyfits
#from pyfits import *
from numpy import radians, degrees, sin, cos, arctan2, hypot
from astropy.io import fits

#python data_data_srh.py 1 0 0 for a single cpu 


## If you are going to use 
## k3match, uncomment this line 
#from k3match import *
#################################################################

def angleto3Dradius(angle, isDegree=True):
  
  if isDegree:
    angle = angle*pi/180.
  
  
  return sqrt((sin(angle))**2 + (1-cos(angle))**2)

#################################################################

      
if __name__ == '__main__':
  
  rpath='redacted for security'
  dpath = 'redacted for security'
  dfilename = dpath+'eBOSS_QSO_NGC_pip_v7_2.dat_z0.8_z2.2_withS_withPIX.fits'
  rfilename = rpath+'eBOSS_QSO_NGC_pip_v7_2.ran_z0.8_z2.2_withS_withPIX.fits'
  
  d=fitsio.read(dfilename)
  dra = d['RA']  
  ddec = d['DEC']
  dpix = d['PIXEL']
  dw = d['WEIGHT_SYSTOT']*d['WEIGHT_FKP']*d['WEIGHT_NOZ']
  ds = d['s']

  r = fitsio.read(rfilename)
  rra = r['RA']
  rdec = r['DEC']
  rpix = r['PIXEL']
  rw = r['WEIGHT_SYSTOT']*r['WEIGHT_FKP']*r['WEIGHT_NOZ']
  rs = r['s']
  # Python equivalent of matchlength
  radius_in_deg = 10.  # degrees
  
  
  size = len(dra)
  doJob = True
  p = 1000
 
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

      ra_small = dra[low:up]
      dec_small = ddec[low:up]
      w_small = dw[low:up]
      

    
      print("No.: ", n, "Go ... ")
      
      
  
      time0 = time()
      
      #############

      
      (m1,m2,distance) = spherematch.match_radec(ra_small, dec_small, rra, rdec, radius_in_deg,notself=True)


      ####################
      ## For k3match
      ## Uncomment k3match-import section at the top of this file 
      ## All tested for peorformance and accuracy by E.Kourkchi
      ## (August 03 2015)
      #####################
      ##(m1, m2, d) =  celestial(ra_small, dec_small, ra, dec, radius_in_deg)



      #############
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
        col4 = fits.Column(name='w1', format='D', array=dw[m1_2])
        col5 = fits.Column(name='w2', format='D', array=rw[m2])
        col6 = fits.Column(name='pix1',format='K',array=dpix[m1_2])
        col7 = fits.Column(name='pix2',format='K',array=rpix[m2]) 
        col8 = fits.Column(name='s1',format='D',array=ds[m1_2])
        col9 = fits.Column(name='s2',format='D',array=rs[m2])
        col10 = fits.Column(name='w1w2',format='D',array=dw[m1_2]*rw[m2])
        cols = fits.ColDefs([col1, col2, col3,col4,col5,col6,col7,col8,col9,col10])


        tbhdu = fits.BinTableHDU.from_columns(cols)
	
	# clobber = True  (it overwrites/updates)
        tbhdu.writeto('./DR_QSO_NGC_v7_2/spherematch.'+str(int(n))+'_DR.fits', clobber=True)
      
      #############          
      
      print('')
      print("run_time = ", (time()-time0)/60., 'min')
      print("# of total matches:", len(d))
      print("Catalog size:", size)
      print("lower index: ", low)
      print("upper index: ", up)
      print('')
  
      #### test
      #ra_max = ra[m1[where(d == max(d))]]
      #dec_max = dec[m1[where(d == max(d))]]
      #print "Max Distance:", max(d), ra_max, dec_max
  
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
  
plot =True

if plot: 

    import matplotlib.pylab as plt

    a,h=fitsio.read(outpath+'NGC_QSO_v7_2_AngUpWeight_lookupTable_upto'+str(max_ang)+'_'+str(binsize)+'_ALLTARGETS.fits',header=True)

    mids = (a['theta_min']+a['theta_max'])/2


    plt.plot(mids,a['dd_par']/(a['dd_fib']),color='blue')
    plt.plot(mids,a['AngUpWeight_wpip'],color='orange')

    plt.legend(['AngUpWeight (without PIP)','AngUpWeight_wpip (with PIP)'])
    plt.xlabel('Angular separation [degrees]')
    plt.ylabel('Angular upweight of DD pairs')

    plt.savefig(outpath+'Angular_upweight'+str(max_ang)+'_'+str(binsize)+'_ALLTARGETS.png')
    plt.show()
  
    
