import numpy as np
import fitsio
import matplotlib.pylab as plt
import pylab as      pl
from astropy.table import Table
import os, sys, glob
from astropy.table import Table,Column
from time import time
#### this directory the pairs of targets from 0 to 10 degress for the quasars in  NGC_QSO v7_2  at 0.8<z<2.2
path = 'redacted for security'
DDpath = path+'/DD_QSO_NGC_v7_2_319636/'


max_ang = 0.1
binsize=0.003

start = time()
give_arr=True

if give_arr: #not os.path.exists(path+'/NGC_QSO_v7_2_pairs_stat_sep.fits'):    # True

    cat = path+'/eBOSS_QSO_NGC_pip_v7_2.dat_319636_withS.fits'
    qso = fitsio.read(cat)
    z = qso['Z']
    ebid = qso['EBOSS_TARGET_ID']
    ra = qso['RA']
    dec = qso['DEC']
    status = qso['STATUS']
    files = glob.glob(DDpath+'*.fits')

    
    ledges = np.round(np.arange(0,max_ang,binsize),3)
    hedges = np.round(np.arange(binsize,binsize+max_ang,binsize),3)
    mids = np.round((ledges+hedges)/2,3)
    #print(list(zip(list(ledges),list(hedges),list(mids))))
    nbins = len(mids)

    dd_par = np.zeros(len(mids))
    dd_fb = np.zeros(len(mids))
    pip_fb = np.zeros(len(mids))

    for f in files:
    
        print('Working on ',f)
        batch = f.split('spherematch.')[1].split('_')[0]
        a = fitsio.read(f)
        ind1 = a['index1']
        ind2 = a['index2']
        
        stat1 = status[ind1]
        stat2 = status[ind2]
        angsep = a['dist']
        pips = []

        for i in range(nbins):
            
            wt = (ledges[i] < angsep) & (angsep <= hedges[i])
            i1 = ind1[wt]
            i2 = ind2[wt]
            
            dd_par[i] += np.sum(wt)
            ## both members are in clustering catalog
            #w_fb = (((2**1 & stat1[wt]) !=0) & ((stat2[wt] & 2**1) !=0))
            # both members got a fiber
            w_fb = (((2**0 & stat1[wt]) !=0) & ((stat2[wt] & 2**0) !=0))
            
            if np.sum(w_fb) >0 :

                ind1w = i1[w_fb]
                ind2w = i2[w_fb]

                for p in range(np.sum(w_fb)):

                    c = (qso['WEIGHT_BW'][ind1w[p]] & qso['WEIGHT_BW'][ind2w[p]])
                    cbit = [bin(ci) for ci in c]
                    nones = sum([cb.count('1') for cb in cbit])

                    if nones >0:
                        pips.append(1860/nones)
                    else:
                        pips.append(0)

            #### to get a taste of the measurement
            #if i==0:
                 
            #     idradecz1= list(zip(i1,ra[i1],dec[i1],z[i1],status[i1]))
            #     idradecz2= list(zip(i2,ra[i2],dec[i2],z[i2],status[i2]))

            #     for i in range(11):

            #          print(idradecz1[i],'\t', idradecz2[i],'\t','{:>8}'.format(angsep[i]),'\t','{:>8}'.format(np.round(pips[i],4)))

 
            dd_fb[i] += np.sum(w_fb)
            pip_fb[i] += np.sum(pips) 
            print('Finished bin #{} in batch {}:  {} < theta[deg] < {}: DD_pair/DD_fib = {}, total upweight for the bin so far: {},{}'.format(i,batch,ledges[i],hedges[i],np.sum(wt)/np.sum(w_fb),dd_par[i]/dd_fb[i], dd_par[i]/pip_fb[i]))
       
       # sephist,b = np.histogram(angsep,bins = [0,0.1,.3,.4,.5])
       # print(sephist,b)
   
angup = dd_par/dd_fb      
angupw = dd_par/pip_fb

write_table = True

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

    cc = Column(angupw,name='AngUpWeight_wpip')
    a.add_column(cc,index=0)

    cc = Column(pip_fb,name='pipsforDD_fib')
    a.add_column(cc,index=0)

    a.write(path+'NGC_QSO_v7_2_AngUpWeight_lookupTable_upto'+str(max_ang)+'_'+str(binsize)+'_319636.fits')
    
print('Finished writing the lookup Table of the angular upweights! took {} min'.format((time()-start)/60.)) 
    
   
