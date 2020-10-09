import numpy as np
from time import time
import fitsio, os, sys
from astropy.table import Column
from astropy.table import Table
from math import cos,pi,log10
import matplotlib.pylab as plt
import pylab as py
from astrometry.libkd import spherematch
#gmpy is A C-coded Python extension module that wraps the GMP library.
import gmpy

''' Example of runninf syntax: python angupweight.py 'eBOSS_LRG_SGC_pip_v7_2.dat.fits' 'log' 0.3

    NOTE: - Make sure the target catalog that is used is not cut to any specofic redshift range 
          - Make sure that the catalog is not limited to targets from the clustering set or the
          targets that have been spectroscopically identifed as the desired target class (i.e., ELG, LRG, QSO). 
          Keep everything given in the parent catalog.
          
          examples: 
              'eBOSS_LRG_SGC_pip_v7_2.dat.fits'
              'eBOSS_ELG_NGC_pip_v7.dat.fits'
              'eBOSS_QSO_SGC_pip_v7_2_new.dat.fits'
'''


path='/utahpath/'

outpath = path

inputfile = sys.argv[1]

filename = path+inputfile

binning = sys.argv[2]

tgt = filename.split(path)[1].split('_pip')[0]

catalog,h = fitsio.read(filename, header = True)

max_ang = sys.argv[3]

if binning=='log':
     binsize=''
     hedges = np.logspace(-2,np.log10(max_ang),30)
     ledges = np.asarray([0]+list(hedges[0:len(hedges)-1]))
     seps = hedges 

elif binning=='linear':
     binsize=0.003
     ledges = np.round(np.arange(0,max_ang,binsize),3)
     hedges = np.round(np.arange(binsize,binsize+max_ang,binsize),3)
     mids = np.round((ledges+hedges)/2,3)

avg_pip_bin = np.zeros(len(hedges))

DD_fib = np.zeros(len(hedges))

DD_par = np.zeros(len(hedges))

wDD_angup = np.zeros(len(hedges))

time0 = time()

for k,sep in enumerate(seps):
      
      print('----------------------------------------')
      print('working on {},{}'.format(k,sep))

      (match1,match2,distance12) = spherematch.match_radec(catalog['RA'], catalog['DEC'], catalog['RA'], catalog['DEC'], sep, notself=True)

      a = (match1 != match2) & (distance12 > ledges[k])

      match1=match1[a]

      match2=match2[a]
      
      num_par = len(match1)
      
      print('number of pairs in target catalog within ', sep, ' degrees = ',num_par)

      a = (catalog['FIBER'] | catalog['CLUSTERING']) !=0 
      
      catalog = catalog[a]
      
      print('working with {} matches ...'.format(np.sum(a)))
      
      (match1,match2,distance12) = spherematch.match_radec(catalog['RA'], catalog['DEC'], catalog['RA'], catalog['DEC'], sep, notself=True)

      a = (match1 != match2) & (distance12 > ledges[k])
      
      match1=match1[a]

      match2=match2[a]

      num_fib = len(match1)

      print('number of DDfib counts within ', sep, ' degrees = ',num_fib)
      
      PIP_fib=np.zeros(num_fib)

      for i in range(num_fib):

          for j in range(60):

              a = (catalog[match1[i]]['WEIGHT_BW'][j]) & (catalog[match2[i]]['WEIGHT_BW'][j]) 

              PIP_fib[i] += gmpy.popcount(int(a))

      
      DD_fib_PIP = np.sum(1860./PIP_fib)

      print('wDD_ang= DD_par/DD_fib_PIP= ',num_par/DD_fib_PIP)  ### definition in equation 9 of PIP paper

      avg_pip_bin[k] = np.sum(1860./PIP_fib)

      DD_fib[k] = num_fib
      
      DD_par[k] = num_par
      
      wDD_angup[k] = num_par/DD_fib_PIP 
     
write_table = True

if write_table:

    a = Table()

    aa = Column(seps,name='theta_in_deg')
    a.add_column(aa,index=0)

    aa = Column(DD_par,name='DD_par')
    a.add_column(aa,index=0)

    bb = Column(DD_fib,name='DD_fib')
    a.add_column(bb,index=0)

    cc = Column(avg_pip_bin,name='avg_pip_inbin')
    a.add_column(cc,index=0)

    dd = Column(wDD_angup,name='wDD_angup')
    a.add_column(dd,index=0)

    outname = tgt+'_AngUpWeight_lookupTable_upto'+str(max_ang)+'_'+str(binsize)+'_ALLTARGETS_inPython.fits'

    if os.path.exists(outpath+outname):

          os.system('rm '+outpath+outname)

    a.write(outpath+outname)

  

print("run_time = ", (time()-time0)/60., 'min')

from time import gmtime, strftime

print("Finished at:", strftime("%Y-%m-%d %H:%M:%S", gmtime()))

