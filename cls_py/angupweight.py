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
from functools import reduce
from astrometry.libkd import spherematch

#gmpy is A C-coded Python extension module that wraps the GMP library.
import gmpy

#### To RUN: python angupweight.py 'eBOSS_LRG_SGC_pip_v7_2.dat.fits' 'log'

path='/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/July2020/'

outpath = path

inputfile = sys.argv[1]


''' target file should not be cut to any redshift range 
    examples: 'eBOSS_LRG_SGC_pip_v7_2.dat.fits'
              'eBOSS_ELG_NGC_pip_v7.dat.fits'
              'eBOSS_QSO_SGC_pip_v7_2_new.dat.fits'
'''

filename = path+inputfile

binning = sys.argv[2]

tgt = filename.split(path)[1].split('_pip')[0]

catalog,h = fitsio.read(filename, header = True)

max_ang = 0.3#3.0

if binning=='log':
     binsize=''
elif binning=='linear':
     binsize=0.003

#ledges = np.round(np.arange(0,max_ang,binsize),3)
#hedges = np.round(np.arange(binsize,binsize+max_ang,binsize),3)
#mids = np.round((ledges+hedges)/2,3)

hedges = np.logspace(-2,np.log10(max_ang),30)
ledges = np.asarray([0]+list(hedges[0:len(hedges)-1]))
seps = hedges #[0.003]

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

      #print('PIP-weighted pair counts is ', np.sum(1860./PIP_fib))

      #print('Average PIP weight for this bin:',np.sum(1860./PIP_fib)/num_fib )# ????? thi sis not what equation 9 of PIP paper 

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

