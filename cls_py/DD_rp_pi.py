from glob import glob
import time, numpy as np
import fitsio, matplotlib.pyplot as plt


##### Shadab's binning 
#rper (rp)  in log bin with   bin edges as 
logrper=np.linspace(-1.5,1.477,25)
#rpar (pi) in linear bin with bin edges by 
rpar=np.linspace(0,picut,picut+1)



### Getting the table with stored angular upweights for bins in angular separations 
reg = 'NGC'
tgt = 'LRG'
angupfile = 'eBOSS_'+tgt+'_'+reg+'_AngUpWeight_lookupTable_upto0.0100000_1.12210_ALLTARGETS_new.fits' 
path = '/utahpath/July2020/'

wei,h = fitsio.read(path+angupfile, header=True)

wdd = wei['WDD_ANGUP']
lang = wei['LANG']
hang = wei['HANG']


chunks = glob(path+'/DD_NGC_LRG_z0.6_z1.0/*.fits')

nbin = len(logrper)-1
edges = 10**logrper

picut = 1000
dd_rppi = np.zeros((picut,nbin),float)
wdd_rppi = np.zeros((picut,nbin),float)
wpipdd_rppi = np.zeros((picut,nbin),float)

for f in chunks:
    print(f)
    dd,h = fitsio.read(f,header=True)    
    pii = abs(dd['s1']-dd['s2'])

    ''' NOTE: To avoid unnecessary looping, we can limit the pairs in each chunk to only pairs with
       pi < 40Mpc/h that is the upper cut on comoving separations along line of sight:
    pimax = (pii <= picut)
    dd = dd[pimax]
    '''
    theta = dd['dist']*np.pi/180.
    pii = abs(dd['s1']-dd['s2'])
    rp = ((dd['s1']+dd['s2'])/2.)*theta
    
    ########  Weights 
#     print(np.min(dd['PIP'], np.max(dd['']))
    PIP = 1860./dd['PIP']  ### equation 6 of Faizan's paper
    wsys = dd['w1w2']   ### (w1_tot*w1_sys*w1_FKP) x (w2_tot*w2_sys*w2_FKP)
    
    weight = wsys*PIP
        
    DD = np.zeros((picut,nbin), float)
    wDD = np.zeros((picut,nbin), float)
    wpipDD = np.zeros((picut,nbin), float)
    
    DDcount = np.zeros(nbin)
    
    for j in range(nbin):
        
        t0 = time.time()

        #print('working on bin #{}'.format(j))
        wrp = (rp < edges[j+1]) & (rp > edges[j])
        cnt = np.sum(wrp)
        DDcount[j] += cnt

        if cnt > 0:
             
            a = pii[wrp]
            allwei = weight[wrp] 
            syswei = wsys[wrp] 
            upws = np.ones(len(a))
            
            for i in range(len(lang)):
   
                ii = (dd['dist'][wrp] < hang[i]) & (dd['dist'][wrp] > lang[i])
                upws[ii] *= wdd[i] 

            warr = np.zeros(picut)
            wpiparr = np.zeros(picut)
            arr = np.zeros(picut)
                  
            for k in range(len(a)): 
            
                warr[int(a[k])] += 1*syswei[k]*upws[k]
                wpiparr[int(a[k])] += 1*allwei[k]*upws[k]
                arr[int(a[k])] += 1

            DD[:,j] = arr
            wDD[:,j] = warr
            wpipDD[:,j] = wpiparr

    dd_rppi = dd_rppi + DD
    wdd_rppi = wdd_rppi + wDD
    wpipdd_rppi = wpipdd_rppi + wpipDD
    
np.savez_compressed(path+'DD_rp_pi_NGC_LRG_raw.npz', array = dd_rppi)
np.savez_compressed(path+'DD_rp_pi_NGC_LRG_angup.npz', array = wdd_rppi)
np.savez_compressed(path+'DD_rp_pi_NGC_LRG_angupPIP.npz', array = wpipdd_rppi)
 
