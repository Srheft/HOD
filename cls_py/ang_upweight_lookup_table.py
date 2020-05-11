import numpy as np
import fitsio
import matplotlib.pylab as plt
import pylab as      pl
from astropy.table import Table
import os, sys, glob
from astropy.table import Table,Column
from time import time
#### this directory the pairs of targets from 0 to 10 degress for the quasars in  NGC_QSO v7_2  at 0.8<z<2.2

path = 'for security reasons this is redacted in public codes'

DDpath = path+'/DD_QSO_NGC_v7_2/'


max_ang = 1.0
binsize=0.003

start = time()
#give_arr=True
write_table = False # default

#if give_arr: 
if not os.path.exists(path+'NGC_QSO_v7_2_AngUpWeight_lookupTable_upto'+str(max_ang)+'_'+str(binsize)+'.fits'):    # True
    
    write_table = True

    cat = path+'/eBOSS_QSO_NGC_pip_v7_2.dat_z0.8_z2.2_withS_withPIX.fits'
    qso = fitsio.read(cat)

    status = qso['STATUS']
    files = glob.glob(DDpath+'*.fits')

    
    ledges = np.round(np.arange(0,max_ang,binsize),3)
    hedges = np.round(np.arange(binsize,binsize+max_ang,binsize),3)
    mids = np.round((ledges+hedges)/2,3)
    nbins = len(mids)

    dd_par = np.zeros(len(mids))
    dd_fb = np.zeros(len(mids))
   
    for f in files:
    
        print('Working on ',f)
        a = fitsio.read(f)
        ind1 = a['index1']
        ind2 = a['index2']
        
        stat1 = status[ind1]
        stat2 = status[ind2]
        angsep = a['dist']

        for i in range(nbins):

            wt = (ledges[i] < angsep) & (angsep <= hedges[i])

            dd_par[i] += np.sum(wt)
    
            # both members of the pair have to have gotten a fiber --> i.e., the 0 bit has to be flipped
            w_fb = (((2**0 & stat1[wt]) !=0) & ((stat2[wt] & 2**0) !=0))
          
            dd_fb[i] += np.sum(w_fb)

            print('Finished bin #{}:  {} < theta[deg] < {}: DD_pair/DD_fib = {}, total upweight for the bin so far: {}'.format(i,ledges[i],hedges[i],np.sum(wt)/np.sum(w_fb),dd_par[i]/dd_fb[i]))
       
   
angup = dd_par/dd_fb      



if write_table:

    a = Table()

    aa = Column(ledges,name='theta_min')
    a.add_column(aa,index=0)

    aa = Column(hedges,name='theta_max')
    a.add_column(aa,index=0)

    aa = Column(dd_par,name='dd_par')
    a.add_column(aa,index=0)

    bb = Column(dd_fb,name='dd_fib')
    a.add_column(bb,index=0)
    
    cc = Column(angup,name='AngUpWeight')
    a.add_column(cc,index=0)

    a.write(path+'NGC_QSO_v7_2_AngUpWeight_lookupTable_upto'+str(max_ang)+'_'+str(binsize)+'.fits')
    
    print('Finished writing the lookup Table of the angular upweights! took {} min'.format((time()-start)/60.)) 

 else:
    print('The requested look up table already exist here: {}'.format(path+'NGC_QSO_v7_2_AngUpWeight_lookupTable_upto'+str(max_ang)+'_'+str(binsize)+'.fits'))
   
