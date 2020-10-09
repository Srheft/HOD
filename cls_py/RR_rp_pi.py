from glob import glob
import time, fitsio, numpy as np
 

picut = 100
path = '/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/July2020/'
##### Shadab's binning 
#rper (rp)  in log bin with   bin edges as 
logrper=np.linspace(-1.5,1.477,25)
#rpar (pi) in linear bin with bin edges by 
rpar=np.linspace(0,picut,picut+1)


chunks = glob(path+'/RR_NGC_LRG_z0.6_z1.0_downsampled/*.fits')

nbin = len(logrper)-1
edges = 10**logrper

rr_rppi = np.zeros((picut,nbin),float)
wrr_rppi = np.zeros((picut,nbin),float)



for n,f in enumerate(chunks):
    
    print(n,f)
    rr,h = fitsio.read(f,header=True)    
    pii = abs(rr['s1']-rr['s2'])
    
    pimax = (pii <= picut)
    rr = rr[pimax]
    
    theta = rr['dist']*np.pi/180.
    pii = abs(rr['s1']-rr['s2'])
    rp = ((rr['s1']+rr['s2'])/2.)*theta
    wsys = rr['w1w2']   ### (w1_tot*w1_sys*w1_FKP) x (w2_tot*w2_sys*w2_FKP)
    weight = wsys
    RR = np.zeros((picut,nbin), float)
    wRR = np.zeros((picut,nbin), float)

    for j in range(nbin):
        
        t0 = time.time()

        #print('working on bin #{}'.format(j))
        wrp = (rp < edges[j+1]) & (rp > edges[j])
        cnt = np.sum(wrp)
        
        if cnt > 0: 
            a = pii[wrp]
            allwei = weight[wrp] 

            arr = np.zeros(picut)
            warr = np.zeros(picut)

            for k in range(len(a)): 
            
                warr[int(a[k])] += 1*allwei[k]
                arr[int(a[k])] += 1
                
            RR[:,j] = arr
            wRR[:,j] = warr

        #print('Took {} min'.format((time.time()-t0)/60.))

    rr_rppi = rr_rppi + RR
    wrr_rppi = rr_rppi + wRR

    
np.savez_compressed(path+'RR_rp_pi_NGC_LRG_raw_10deg.npz', array = rr_rppi)
np.savez_compressed(path+'RR_rp_pi_NGC_LRG_10deg.npz', array = wrr_rppi)
    
