import healpy as hp
import numpy as np
import fitsio
from glob import glob
import sys

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

def getXi_wp(ddfile, drfile, rrfile, picut, verbose=False):

    
    dd_rppi = np.load(ddfile)['array']
    dr_rppi = np.load(drfile)['array']
    rr_rppi = np.load(rrfile)['array']
    
    n = np.sqrt(np.sum(rr_rppi)/np.sum(dd_rppi))  # r/d
    m = np.sqrt(np.sum(dr_rppi)/np.sum(dd_rppi))  # r/d
    z = np.sqrt(np.sum(rr_rppi)/np.sum(dr_rppi))  # r/d
    
    xi = (dd_rppi/rr_rppi)*n**2 - 1
    xi_LS = (dd_rppi/rr_rppi)*n**2  - 2*z*(dr_rppi/rr_rppi) + 1.

    wp = np.zeros_like(xi[0,:])
    wp_LS = np.zeros_like(xi_LS[0,:])
    rpbinnedRRcount = np.zeros_like(rr_rppi[0,:])

    for i in range(xi.shape[0]):

        wp += xi[i,:]
        wp_LS += xi_LS[i,:]
        rpbinnedRRcount += rr_rppi[i,:]
        
    return wp, wp_LS, xi, xi_LS, rpbinnedRRcount 
########################################################################### Getting Cij

def get_cij(i, j, pimax, npix, fullwp, pixeled_wp, rpbinnedRRcount, rpbinpixRRcount):
    
    import sys
    
    # sanity check
    if (i > len(fullwp) | j > len(fullwp)):
        
        print('Error: the i or j index are beyond the size of measured correlation!')
        sys.exit()
    cij = 0
    
    for L in range(npix):
        
        cij += (rpbinpixRRcount[L][i]/rpbinnedRRcount[i])*(pixeled_wp[L][i]-fullwp[i])*(rpbinpixRRcount[L][j]/rpbinnedRRcount[j])*(pixeled_wp[L][j]-fullwp[j])
        
    return cij    
       

######################################################################### Getting pixel number

nside = 8
path = '/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/July2020/'
random = 'eBOSS_LRG_NGC_pip_v7_2_ran_withS.fits'

r = fitsio.read(path+random)
rra = r['RA']
rdec = r['DEC']

pix, uniqpix, weights = get_pixnum(rra,rdec,nside)

picut = 40


RRpath = '/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/July2020/jkpix_LRG_NGC_downsampled/'

rrs = glob(RRpath+'/RR*_raw_pix*.npz')

DDpath = '/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/July2020/jkpix_LRG_NGC/'

dds = glob(DDpath+'/DD*_angup_pix*.npz')


#######################################################################   Getting pixelated  wp 
def get_pixeled_wp(RRpath,DDpath,uniqpix):

    logrper=np.linspace(-1.5,1.477,25)

    nbin = len(logrper)-1

    pixeled_wp = np.zeros((len(uniqpix),nbin),float)


    for j,p in enumerate(uniqpix):
    
        rrfile = RRpath+'/RR_rp_pi_NGC_LRG_angup_pix'+str(int(p))+'.npz'
        ddfile = DDpath+'/DD_rp_pi_NGC_LRG_raw_pix'+str(int(p))+'.npz' 
        drfile = rrfile
        wp, wp_LS, xi, xi_LS,rpbinnedRRcount = getXi_wp(ddfile, drfile, rrfile, picut, verbose=False)
    
        pixeled_wp[j,:] = wp
    return pixeled_wp
######################################################################  Getting full wp before pixelation

path = '/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/July2020/'
ddfile = path + 'DD_rp_pi_NGC_LRG_angup.npz' 
drfile = path + 'DR_rp_pi_NGC_LRG_raw_downsampled.npz'
rrfile = path + 'RR_rp_pi_NGC_LRG_downsampled.npz'

fullwp, wp_LS, xi, xi_LS, rpbinnedRRcount = getXi_wp(ddfile, drfile, rrfile, picut, verbose=False)


##########################################################################   Getting RR pair count for all rp bins of each pixel
rpbinpixRRcount = np.load(path + "Pixeled_RRcounts_LRG_NGC_86pix24rpbins.npz")['array']

npix = len(uniqpix)

cov = np.zeros((nbin,nbin),float)

for i in range(nbin):
   
    for j in range(nbin):

        cov[i][j] = get_cij(i, j, picut, npix, fullwp, pixeled_wp, rpbinnedRRcount, rpbinpixRRcount)

np.savez_compressed(path+'Covarianve_LRG_NGC_'+str(len(uniqpix))+'pix'+str(nbin)+'rpbins.npz', array = cov)



chunks = glob(path+'/RR_NGC_LRG_z0.6_z1.0_downsampled/*.fits')

logrper=np.linspace(-1.5,1.477,25)
rpar=np.linspace(0,picut,picut)

nbin = len(logrper)-1
edges = 10**logrper

picut = 40

rpbinnedRRcount = np.zeros((len(uniqpix),nbin))


for i,p in enumerate(uniqpix[70:86]):

    rpPIpixRRcount = np.zeros((picut,nbin))

    pixweit = weights[i]

    print('Excluding pixel {}, with weight {} '.format(p,pixweit))

    for f in chunks:

        rr,h = fitsio.read(f,header=True)

        m1 = rr['index1']
        m2 = rr['index2']
        pixnum1 = pix[m1]
        pixnum2 = pix[m2]

        ind = (pixnum2 == p) | (pixnum1 == p)
        
        RR = np.zeros((picut,nbin))

        if np.sum(~ind) > 0:

            rr_pix = rr[~ind]

            pii = abs(rr_pix['s1']-rr_pix['s2'])

            pimax = (pii <= picut)

            rr_pix = rr_pix[pimax]

            theta = rr_pix['dist']*np.pi/180.

            pii = abs(rr_pix['s1']-rr_pix['s2'])

            rp = ((rr_pix['s1']+rr_pix['s2'])/2.)*theta

            for j in range(nbin):

                wrp = (rp < edges[j+1]) & (rp > edges[j])

                cnt = np.sum(wrp)

                rpbinnedRRcount[i,j] += cnt

                if cnt > 0:
 
                    a = pii[wrp]

                    arr = np.zeros(picut)

                    for k in range(len(a)): 
            
                        arr[int(a[k])] += 1
                
                RR[:,j] = arr

        rpPIpixRRcount = rpPIpixRRcount + RR

        np.savez_compressed(path+'Pixeled_RRcounts_LRG_NGC_pix'+str(p)+'_'+str(picut)+'pi_'+str(nbin)+'rpbins.npz', array = rpPIpixRRcount)

    np.savez_compressed(path+'Pixeled_RRcounts_LRG_NGC_'+str(len(uniqpix))+'pix'+str(nbin)+'rpbins.npz', array = rpbinnedRRcount)

    print('finished this pix')


##############################################

def make_covariance(picut, npix, fullwp, pixeled_wp, rpbinnedRRcount, pixRRcount):

     # It needs to be checked that the number of pixels (npix) and the maximum pi reach (pimax)
     # matches the size of the pixeled_wp, rpbinnedRRcount, pixRRcount arrays
     # pixeled_wp, a 86 x 24 size array  -->made by get_pixeled_wp() (exists in this code )
     # fullwp --> a 24 size array --> using  getXi_wp() function
     # rpbinnedRRcount, a 24 size array --> using  getXi_wp() function 
     # pixRRcount a 86 x 24 array ---> the code above this unit 
     # npix = len(uniqpix)
     # 

     cov = np.zeros((nbin,nbin))

    for i in range(nbin):

        for j in range(nbin):

            cov[i][j] = get_cij(i, j, picut, npix, fullwp, pixeled_wp, rpbinnedRRcount, pixRRcount)


    ncov = np.zeros((nbin,nbin))

    for i in range(nbin):

        for k in range(nbin):
                        
            ncov[i,k] = cov[i,k]/np.sqrt(cov[i,i])/np.sqrt(cov[k,k])


    np.savez_compressed(path+'Covarianve_LRG_NGC_'+str(len(uniqpix))+'pix'+str(nbin)+'rpbins.npz', array = cov)
    np.savez_compressed(path+'NormCovarianve_LRG_NGC_'+str(len(uniqpix))+'pix'+str(nbin)+'rpbins.npz', array = ncov)

################################################


def plot_cov(input_cov):


    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm
    import pylab as py


    fig = plt.figure(figsize=(5, 5), dpi=100)
    ax = fig.add_axes([0.1,0.1,0.8,0.8])  

    from matplotlib.colors import LogNorm

    ax.imshow(ncov)#, norm=LogNorm(vmin=0.01, vmax=1000000))

    ax.savefig(path+'Normalized_covar_LRG_NGC.png')



