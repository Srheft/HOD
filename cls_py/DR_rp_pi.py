from glob import glob
import time, fitsio, numpy as np

picut = 40

##### Shadab's binning 
#rper (rp)  in log bin with bin edges as 
logrper=np.linspace(-1.5,1.477,25)
#rpar (pi) in linear bin with bin edges by 
rpar=np.linspace(0,picut,picut+1)

### Getting the table with stored angular upweights for bins in angular separations 

reg = 'NGC'
tgt = 'LRG'
angupfile = 'eBOSS_'+tgt+'_'+reg+'_AngUpWeight_lookupTable_upto0.0100000_0.105929_ALLTARGETS_DR.fits' 
path = '/utahpath/July2020/'

wei,h = fitsio.read(path+angupfile, header=True)

wdr = wei['WDR_ANGUP']
lang = wei['LANG']
hang = wei['HANG']

##########################################################################

#chunks = glob(path+'/DR_NGC_LRG_z0.6_z1.0/*.fits')
chunks = glob(path+'/DR_NGC_LRG_z0.6_z1.0_downsampled/*.fits')

nbin = len(logrper)-1
edges = 10**logrper

dr_rppi = np.zeros((picut,nbin),float)
wdr_rppi = np.zeros((picut,nbin),float)
wsysdr_rppi = np.zeros((picut,nbin),float)


for f in chunks:
    print(f)
    dr,h = fitsio.read(f,header=True)    
    pii = abs(dr['s1']-dr['s2'])
    ''' NOTE: To avoid unnecessary looping, we can limit the pairs in each chunk to only pairs with
       pi < 40Mpc/h that is the upper cut on comoving separations along line of sight:

     pimax = (pii <= picut)
     dr = dr[pimax]
    '''    
    theta = dr['dist']*np.pi/180.
    pii = abs(dr['s1']-dr['s2'])
    rp = ((dr['s1']+dr['s2'])/2.)*theta
    PIP = 1860./dr['PIP']  ### equation 6 of Faizan's paper
    wsys = dr['w1w2']   ### (w1_noz*w1_systot*w1_FKP) x (w2_noz*w2_systot*w2_FKP)
    weight = wsys*PIP
    
    
    DR = np.zeros((picut,nbin), float)
    wDR = np.zeros((picut,nbin), float)
    wsysDR = np.zeros((picut,nbin), float)

    for j in range(nbin):

#         print('working on bin #{}'.format(j))
        wrp = (rp < edges[j+1]) & (rp > edges[j])
        cnt = np.sum(wrp)
        
        if cnt > 0: 
            a = pii[wrp]
            upws = np.ones(len(a))
            allwei = weight[wrp] 
            syswei = wsys[wrp] 

            for i in range(len(lang)):
   
                 ii = (dr['dist'][wrp] < hang[i]) & (dr['dist'][wrp] > lang[i])
                 upws[ii] *= wdr[i]

            warr = np.zeros(picut)
            arr = np.zeros(picut)
            wsysarr = np.zeros(picut)
                  
            for k in range(len(a)): 
            
                warr[int(a[k])] += 1*syswei[k]*upws[k]   ### AngUp & systematics
                wsysarr[int(a[k])] += 1*syswei[k]   ### only systematic weights were applied w_tot x w_tot'
                arr[int(a[k])] += 1   ### raw counts

            DR[:,j] = arr
            wDR[:,j] = warr
            wsysDR[:,j] = wsysarr
          

    dr_rppi = dr_rppi + DR
    wdr_rppi = wdr_rppi + wDR
    wsysdr_rppi = wsysdr_rppi + wsysDR
    
    

np.savez_compressed(path+'DR_rp_pi_NGC_LRG_raw_downsampled.npz', array = dr_rppi)
np.savez_compressed(path+'DR_rp_pi_NGC_LRG_angupPIP_downsampled.npz', array = wdr_rppi)
np.savez_compressed(path+'DR_rp_pi_NGC_LRG_sysraw_downsampled.npz', array = wsysdr_rppi)
        
