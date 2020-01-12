# HOD
Halo Occupation playground
 random catalogs are at version 7_1 for the quasars, 7 for the LRGs and ELGs (6 only had elg changes).
 
 
 
 #hod_utah_eboss_NGC_v7_trimm_legacy_no-obs
**i)** the 'eBOSS_QSO_hod_utah_eboss_*_v7_trimm_legacy_no-obs.*.fits' contain all objects within the geometry of the clustering catalogues that do not fall in one of the veto masks (except the collision priority veto mask). These are the catalogues to be used to compute the 'angular weights'. The angular weights are the ratio between the raw DD(DR) pair counts (do not use any weight here) of all targets and the PIP-corrected DD(DR) pairs (again do not use any weight other than PIP) from objects that got a fiber in the actual eBOSS targeting (for each target this information is contained in the first bit of the bitwise weights);
 
# 301,680 NGC-QSQ, 3,377,232 Random
# 203,194 SGC-QSO, 2,712,291 Random
# 149,416 NGC-LRG, 4,397,705 Random
# 115,969 SGC-LRG, 4,018,148 Random



#eBOSS_LRG_clustering_NGC_v7_COMP_BOSS.ran  --> 4 files NGC or SGC lrgs and qsos
**ii)** the 'eBOSS_QSO_clustering_*_v7_1_COMP_BOSS.ran.fits' are the randoms to be used for the 2PCF measurements of the eBOSS clustering catalogues. These are exactly the catalogues you can find on wiki except the two SEQUELS chunks that are removed and the WEIGHT_SYSTOT is rescaled by COMP_BOSS since the PIP weights already account for the completeness;

**iii)** 'eBOSS_QSO_clustering_*GC_v7_1_ids.dat.fits': These are slightly different than the catalogues on wiki because I removed the SEQUELS area in NGC, added IMATCH=8 targets and added the EBOSS_TARGET_ID column. [Use EBOSS_TARGET_ID to match objects between the catalogues and the files containing the bitwise weights (the IDs are int64 type integers)].


 #weights_tiled_QSO_eBOSS_NGC_31x60_trimm_promoted_v7_no-obs.dat
 #weights_tiled_QSO_eBOSS_SGC_31x60_trimm_promoted_v7_no-obs.dat
 
  These are ASCII files with columns EBOSS_TARGET_ID, RA, DEC, Z (neglect Z) and a sequence of 60 columns that contain the bitwise weights resulting from the 1860 fiber assignment runs. The first bit in the first of these 60 columns is the outcome of the actual eBOSS tiling run while the remaining bits are random realizations of the tiling. 
