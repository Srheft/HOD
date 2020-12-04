I - go to work area:

/uufs/chpc.utah.edu/common/home/astro/dawson/sarahE/eboss/MultiHODFitter/

II - In your .bashrc,  export and then echo these env variable:

     $datadir
     $jnroot
     $xi2droot
     $xi02root
     $pcroot

III - edit the config.ini file (config_sarah.ini):
   in the last line, edit the jnfile name format for the redshift limits of the sample (zlim is 0.6 and 1.0) 
_________________________________________________________________________
1- Set up the jacknife regions

python Jacknife_region_2D_keys.py -r1 eBOSS_LRG_NGC_pip_v7_2_ran_withS_z0.6_z1.0_downsampled.fits -NjackRA 15 -NjackDEC 6 -ax1 RA -ax2 DE
C -odir $jnroot/eBOSS_LRG_NGC

OUTPUT:

XI/PairCount/lrg-ngc-logrp-pi-JNdir/*.dat ---> 4 sets of 90 pair count files 


_________________________________________________________________________
2- Run the pair count

python Runme_Correlation.py -config_file ./config/config_sarah.ini

OUTPUT:

XI/PairCount/lrg-ngc-All-logrp-pi-RR.dat
XI/PairCount/lrg-ngc-All-logrp-pi-DD.dat
XI/PairCount/lrg-ngc-All-logrp-pi-DR.dat

XI/PairCount/lrg-ngc-All-logrp-pi-norm.dat ---> ?

_________________________________________________________________________
3- Convert pair count to desired measurements

python PairCountTOxi.py -sampmode 3 -plots 1 -pcroot ./XI/PairCount/lrg-ngc -xi2droot ./XI/XI2D/lrg-ngc -xi02root ./XI/WP/lrg-ngc -njn 90


OUTPUT:

wp(rp) for all 90 jacknife resamplings:        XI/WP/lrg-ngc-wp-logrp-pi-NJN-90.txt
xi(rp,pi) for all 90 jacknife resamplings:     XI/XI2D/lrg-ngc-logrp-pi-NJN-90.txt



