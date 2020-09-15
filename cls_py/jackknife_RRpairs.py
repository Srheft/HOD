from glob import glob
import healpy as hp
import fitsio,numpy as np


##########################################################

def get_pixnum(ra,dec,nside):
    
    import numpy as np, healpy as hp 
    
    NPIX = hp.nside2npix(nside)

    theta = 0.5*np.pi-dec*np.pi/180.
    pix=hp.ang2pix(nside,theta,phi=ra*np.pi/180.)

    print('Number of unique pixels the data covers out of {} total pixels with nisde {} : {}'.format(NPIX,nside,len(np.unique(pix))))

    histpix,edges=np.histogram(pix, bins=np.arange(np.max(pix)+1))

    p = (histpix > 0)
    print('Average and median population of the filled pixels: {}, {}'.format(np.round(np.mean(histpix[p]),3), np.median(histpix[p])))
    xpop = (histpix == np.max(histpix)) 
    mpop = (histpix == np.min(histpix[p])) # [p] to ensure the min is not zero

    print('Most populated pixel is {} and has {} random points'.format(np.where(xpop)[0],np.max(histpix[xpop])))
    print('Least poulated pixel is {} and has {} random points'.format(np.where(mpop)[0],np.min(histpix[mpop])))

    avgpop = np.round(np.mean(histpix[p]),3)
    low = (histpix[p] < (2/3.)*np.mean(histpix[p])) 
    print('Number of pixels below the 2-3rd of average population:{}'.format(np.sum(low)))
    wl = np.where(low)[0]
    weights = np.ones(len(np.unique(pix)))
    weights[wl] = np.round(histpix[p][wl]/avgpop, 2)
#     print(len(weights),len(low),len(histpix),weights)
    upix = np.unique(pix)
    return pix,upix, weights


###########################################################
nside = 8
path = '/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/July2020/'
data = 'eBOSS_LRG_NGC_pip_v7_2_new.dat.fits'
random = 'eBOSS_LRG_NGC_pip_v7_2_ran_withS.fits'

d = fitsio.read(path+data)
dra = d['RA']
ddec = d['DEC']

r = fitsio.read(path+random)
rra = r['RA']
rdec = r['DEC']

pix, uniqpix, weights = get_pixnum(rra,rdec,nside)

picut, pimax = 60, 60



###################################
chunks = glob(path+'/RR_NGC_LRG_z0.6_z1.0/*.fits')

logrper=np.linspace(-1.5,1.477,25)
rpar=np.linspace(0,picut,picut)

nbin = len(logrper)-1
edges = 10**logrper



for i,p in enumerate(uniqpix):
    
    rr_rppi_pix = np.zeros((len(rpar),nbin),float)
    wrr_rppi_pix = np.zeros((len(rpar),nbin),float)
    

    pixweit = weights[i]

    print('Excluding pixel {}, with weight {} '.format(p,pixweit))
        
    for f in chunks:

        rr,h = fitsio.read(f,header=True)

        m1 = rr['index1']
        m2 = rr['index2']
        pixnum1 = pix[m1]
        pixnum2 = pix[m2]
#         print(len(pixnum1),len(dr))
       
        
        ind = (pixnum2 == p) | (pixnum1 == p)
        
#         print(np.sum(ind), np.sum(~ind))
        
        if np.sum(~ind) > 0:
            
#             print('Excluding pixel {}, with weight {} and population {}'.format(p,pixweit,np.sum(ind)))
        
            rr_pix = rr[~ind]

            pii = abs(rr_pix['s1']-rr_pix['s2'])

            pimax = (pii <= picut)
            rr_pix = rr_pix[pimax]

            theta = rr_pix['dist']*np.pi/180.
            pii = abs(rr_pix['s1']-rr_pix['s2'])
            rp = ((rr_pix['s1']+rr_pix['s2'])/2.)*theta
            wsys = rr_pix['w1w2']   ### (w1_tot*w1_sys*w1_FKP) x (w2_tot*w2_sys*w2_FKP)
            weight = wsys
    
            RR = np.zeros((picut,nbin), float)
            wRR = np.zeros((picut,nbin), float)
            

            for j in range(nbin):
        #         print('working on bin #{}'.format(j))
                wrp = (rp < edges[j+1]) & (rp > edges[j])
                cnt = np.sum(wrp)

                if cnt > 0: 
        #             print(cnt)
                    a = pii[wrp]
                    syswei = wsys[wrp] 

                    warr = np.zeros(picut)
                    arr = np.zeros(picut)

                    for k in range(len(a)): 

                        warr[int(a[k])] += 1*syswei[k]   ### systematics
                        arr[int(a[k])] += 1   ### raw counts

                    RR[:,j] = arr
                    wRR[:,j] = warr
                   

            rr_rppi_pix = rr_rppi_pix + RR
            wrr_rppi_pix = wrr_rppi_pix + wRR
            


    np.savez_compressed(path+'jkpix_LRG_NGC/RR_rp_pi_NGC_LRG_raw_pix'+str(p)+'.npz', array = rr_rppi_pix)
    np.savez_compressed(path+'jkpix_LRG_NGC/RR_rp_pi_NGC_LRG_angup_pix'+str(p)+'.npz', array = wrr_rppi_pix)
    
    print('finished this pix')
                    
                    
 
