# HOD
Halo Occupation playground



# Catalogs
  moved to version 7.2 of randoms and data for quasars;  


   QSO => 0.8  < z < 2.2
   
   LRG => 0.6  < z < 1.0
   
   ELG => 0.6  < z < 1.1

   NGC_QSOs    191,430
   
   NGC_LRG     149,416
   
   NGC_ELG     83,769

   SGC_QSOs    125,899
   
   SGC_LRG     115,969
   
   SGC_ELG     89,967

LRGs and QSO overlap in 0.8 < z < 1.0 

 #hod_utah_eboss_NGC_v7_trimm_legacy_no-obs
**i)** the 'eBOSS_QSO_hod_utah_eboss_*_v7_trimm_legacy_no-obs.*.fits' contain all objects within the geometry of the clustering catalogues that do not fall in one of the veto masks (except the collision priority veto mask). These are the catalogues to be used to compute the 'angular weights'. The angular weights are the ratio between the raw DD(DR) pair counts (do not use any weight here) of all targets and the PIP-corrected DD(DR) pairs (again do not use any weight other than PIP) from objects that got a fiber in the actual eBOSS targeting (for each target this information is contained in the first bit of the bitwise weights);
 
#301,680 NGC-QSQ, 3,377,232 Random
#203,194 SGC-QSO, 2,712,291 Random
#149,416 NGC-LRG, 4,397,705 Random
#115,969 SGC-LRG, 4,018,148 Random


#eBOSS_LRG_clustering_NGC_v7_COMP_BOSS.ran  --> 4 files NGC or SGC lrgs and qsos
**ii)** the 'eBOSS_QSO_clustering_*_v7_1_COMP_BOSS.ran.fits' are the randoms to be used for the 2PCF measurements of the eBOSS clustering catalogues. These are exactly the catalogues you can find on wiki except the two SEQUELS chunks that are removed and the WEIGHT_SYSTOT is rescaled by COMP_BOSS since the PIP weights already account for the completeness;

**iii)** 'eBOSS_QSO_clustering_*GC_v7_1_ids.dat.fits': These are slightly different than the catalogues on wiki because I removed the SEQUELS area in NGC, added IMATCH=8 targets and added the EBOSS_TARGET_ID column. [Use EBOSS_TARGET_ID to match objects between the catalogues and the files containing the bitwise weights (the IDs are int64 type integers)].


 #weights_tiled_QSO_eBOSS_NGC_31x60_trimm_promoted_v7_no-obs.dat
 #weights_tiled_QSO_eBOSS_SGC_31x60_trimm_promoted_v7_no-obs.dat
 
 # My products
  DD_s_eBOSS_QSO_NGC_v7_1.fits ---> 19 bins from 40 kpc/h till ~200 Mpc/h

  'eBOSS_QSO_clustering_NGC_v7_1_ids_zerosones.fits' ---> converted the base-2 binaries in weights*.dat files to the sequence of 1860 long sequence of zeos and ones and stored it as a string in the last column of this file.---> the problem is resoved: use firmat(bit,"031b") NOT format(bit,"b"). This would ensure equal length for all bits.

 
  These are ASCII files with columns EBOSS_TARGET_ID, RA, DEC, Z (neglect Z) and a sequence of 60 columns that contain the bitwise weights resulting from the 1860 fiber assignment runs. The first bit in the first of these 60 columns is the outcome of the actual eBOSS tiling run while the remaining bits are random realizations of the tiling. 


see the catalog analysis here /git_repo/HOD/catalogs.ipynb
/git_repo/HOD/catalogs.ipynb
/git_repo/HOD/weighted_DD.ipynb

# Findings
Not all quasars in NGC_v7_1 have the same length of number of fiberassign trials (I assumed to get 1860 length after converting the PIP weight into sequence of zeros and ones) yet I get a range of 1768 to 1860. then the multiplication of two of the weights, would have to be cut to the shared length.



