# -*- coding: utf-8 -*-
from time import time
import sys
import os
import numpy as np
from math import *
import pyfits
from pyfits import *
from astropy.io import fits
import glob
from astropy.table import Table, Column 



#################################################################

      
if __name__ == '__main__':
  
  # 32-bit signed integer
  bin_ang = np.zeros((100,), dtype=np.dtype('float64'))
  hist_all_ang = np.zeros((100,), dtype=np.dtype('int64'))
  bin_s = np.zeros((100,), dtype=np.dtype('float64'))
  hist_all_s = np.zeros((100,), dtype=np.dtype('int64'))
  hist_NGC_ang = np.zeros((100,), dtype=np.dtype('int64'))
  hist_NGC_s = np.zeros((100,), dtype=np.dtype('int64'))
  hist_SGC_ang = np.zeros((100,), dtype=np.dtype('int64'))
  hist_SGC_s = np.zeros((100,), dtype=np.dtype('int64')) 
  #for n in range(873):
  for n in range(16235):
  #for filename in glob.glob('./*fits'):

    
    #print filename
    print n
    #n = 16230
    #filename = 'DR_withsyscus_hist/python.out.'+str(int(n))+'.100.hist_d_s_NGC_SGC.fits'
    #filename = 'DD_withsyscus_hist/python.out.'+str(int(n))+'.DD_100.hist_d_s_NGC_SGC.fits'
    filename = 'RR_withsyscus_hist/python.out.'+str(int(n))+'.100.hist_d_s_NGC_SGC.fits'
    a = pyfits.open(filename)
    d = a[1].data
    #mid_ang = d.field('midbins_deg')
    #ALL_ang = d.field('Hist_d_TOT')
    #mid_s = d.field('midbins_s')
    #ALL_s = d.field('Hist_s_TOT')
    #NGC_ang = d.field('Hist_d_NGC')
    #NGC_s = d.field('Hist_s_NGC')
    #SGC_ang = d.field('Hist_SGC')
    #SGC_s = d.field('Hist_SGC_s')    

    mid_ang = d.field('midbins_deg')
    ALL_ang = d.field('Hist_TOT')
    mid_s = d.field('midbins_s')
    ALL_s = d.field('Hist_TOT_s')
    NGC_ang = d.field('Hist_NGC')
    NGC_s = d.field('Hist_NGC_s')
    SGC_ang = d.field('Hist_SGC')
    SGC_s = d.field('Hist_SGC_s')    

    bin_ang += mid_ang
    hist_all_ang += ALL_ang
    bin_s += mid_s
    hist_all_s += ALL_s
    hist_NGC_ang += NGC_ang
    hist_NGC_s += NGC_s
    hist_SGC_ang += SGC_ang
    hist_SGC_s += SGC_s
    
  col1 = fits.Column(name='midbins_deg', format='D', array=bin_ang)  
  col2 = fits.Column(name='AllHist_TOT_ang', format='K', array=hist_all_ang)
  col3 = fits.Column(name='midbins_s', format='D', array=mid_s)
  col4 = fits.Column(name='AllHist_TOT_s', format='K', array=hist_all_s)
  col5 = fits.Column(name='AllHist_NGC_ang', format='K', array=hist_NGC_ang)
  col6 = fits.Column(name='AllHist_NGC_s', format='K', array=hist_NGC_s)
  col7 = fits.Column(name='AllHist_SGC_ang', format='K', array=hist_SGC_ang)
  col8 = fits.Column(name='AllHist_SGC_s', format='K', array=hist_SGC_s)
  cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8])
  tbhdu = fits.BinTableHDU.from_columns(cols)
  #tbhdu.writeto('python.out.ALL.DD.hist.fits', clobber=True)
  #tbhdu.writeto('python.out.ALL.DR.hist.fits', clobber=True)
  tbhdu.writeto('python.out.ALL.RR.hist.fits', clobber=True)
  #print hist_all.dtype
  #print hist_NGC.dtype
  myTable = Table()
  empty = []

  hist_all_ang += ALL_ang
  hist_all_s += ALL_s
  hist_NGC_ang += NGC_ang
  hist_NGC_s += NGC_s
  hist_SGC_ang += SGC_ang 
  hist_SGC_s += SGC_s
  myTable.add_column(Column(data=hist_all_ang,name='AllHist_TOT_ang', dtype=np.dtype('int64')))
  myTable.add_column(Column(data=hist_all_s,name='AllHist_TOT_s', dtype=np.dtype('int64')))
  myTable.add_column(Column(data=hist_NGC_ang,name='AllHist_NGC_ang', dtype=np.dtype('int64')))
  myTable.add_column(Column(data=hist_NGC_s,name='AllHist_NGC_s', dtype=np.dtype('int64')))
  myTable.add_column(Column(data=hist_SGC_ang,name='AllHist_SGC_ang', dtype=np.dtype('int64')))
  myTable.add_column(Column(data=hist_SGC_s,name='AllHist_SGC_s', dtype=np.dtype('int64')))
  myTable.write('python.out.ALL.RR.hist.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
  #myTable.write('python.out.ALL.DR.hist.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
  #myTable.write('python.out.ALL.DD.hist.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
  
