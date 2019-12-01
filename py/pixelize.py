from time import time
import os, sys, fitsio, datetime
from astropy import units
from astropy.io import fits
from astropy.coordinates import SkyCoord
import numpy as np
from math import pi,cos,sin
from astropy.table import Table
import matplotlib.pyplot as plt
from glob import glob
from collections import Counter
import desimodel.focalplane
import desimodel.footprint
from desitarget.targetmask import desi_mask, obsconditions
import random,fitsio
from desimodel.io import (load_desiparams, load_fiberpos, load_platescale,load_tiles, load_deviceloc)
import pylab as py
from matplotlib import gridspec
from matplotlib.patches import Circle
#from desiutil.plots import init_sky, plot_healpix_map, p lot_grid_map, plot_sky_circles, plot_sky_binned, prepare_data


path= '/uufs/chpc.utah.edu/common/home/bolton-group1/bolton_data2/kdawson/sarahE/eboss/cat/'
qncat='eBOSS_QSO_hod_utah_eboss_SGC_v7_trimm_legacy_no-obs.dat.fits'
qscat='eBOSS_QSO_hod_utah_eboss_NGC_v7_trimm_legacy_no-obs.dat.fits'

nq =fitsio.read(path+qncat)
ns = fitsio.read(path+qscat)
q = np.concatenate((nq,ns),axis=0)


import healpy as hp

### set pixeled map resolution 
nside = 4
print("Approximate resolution at nside {} is {:.2} deg".format(
        nside, hp.nside2resol(nside, arcmin=True) / 60))

NPIX = hp.nside2npix(nside)
print('Number of pixels given the resolution:',NPIX)
print((hp.ang2pix(32,45*pi/180,90*pi/180.)))
