# -*- coding: utf-8 -*-
import sys
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
import numpy as np
from time import time, sleep
from numpy import loadtxt, zeros
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import seaborn as sns
import pandas as pd
from scipy.stats import spearmanr
from wpfunc4hod import wpfunc,wpfunc10c
from wp_mcmc_2pars_CORRECTED import median_error
from scipy.interpolate import interp1d

def onesig_envelope(arr):
    
    arr = np.sort(arr)
    size = len(arr)
 
    cen = np.median(arr)
    up = arr[int(round(0.84*size))]
    lo =  arr[int(round(0.16*size))] 
    return (lo,cen,up)

def twosig_envelope(arr):
    
    arr = np.sort(arr)
    size = len(arr)
 
    cen = np.median(arr)
    up = arr[int(round(0.979*size))]
    lo =  arr[int(round(0.021*size))] 
    return (lo,cen,up)

#============================================================================
def meanNM(M,delm,Mm1):
     
    return (1./(sqrt(2.*pi)*delm))*exp(-(ln(M/Mm1)*ln(M/Mm1))/(2.*delm*delm))
#============================================================================

def stats(seed,nchain,pars,plots=False):
  
  filename1 ='FIT1_NOV2018-4KDE+KO12+eBOSS.dat'
  filename2 ='FIT2_NOV2018-allKO12+eBOSS_TRIMMED.dat'
  filename3 ='FIT3_NOV2018-4KDE+eBOSS_TRIMMED.dat'
  filename4 ='10cFIT3_OCT2018-4KDE+eBOSS.dat'

  #print('Working with ', filename1)
  dropouts = 0
  cols1=["n","fsat1", "Mm1","dm","chi2","wp0","wp1","wp2","wp3","wp4","wp5","wp6","wp7","wp8","wp9","wp10","wp11","wp12","wp13","wp14","wp15","wp16","wp17","wp18","wp19","wp20","wp21","wp22","wp23","wp24","wp25","wp26"]
  dat1 = loadtxt(filename1)

  #-----------------------Del m has been the same for al three fits 
  delm = dat1[:,3]
  delm = delm[dropouts:len(dat1)]
  bdelm_err = median_error(delm)
  #---------------------------------------
  
  fsat1 = dat1[:,1]
  Mm1 = dat1[:,2]
  fsat1 = fsat1[dropouts:len(dat1)]
  bfsat1_err1 = median_error(fsat1)
  Mm1 = Mm1[dropouts:len(dat1)]
  bMm1_err1 = median_error(Mm1)
  Mm1 = np.asarray(Mm1)
  fsat1 = np.asarray(fsat1)
  chi2_1 = dat1[:,4]  
  chi2_1 = chi2_1[dropouts:len(dat1)]
  chi2_1 = np.asarray(chi2_1)
  bferr1=median_error(fsat1)
  merr1= median_error(Mm1)   
  w1=np.where(chi2_1 == min(chi2_1))
  bMm1 = np.mean(Mm1)
  bfsat1 = np.mean(fsat1)

  cols2=["n","fsat1", "Mm1","dm","chi2","wp0","wp1","wp2","wp3","wp4","wp5","wp6","wp7","wp8","wp9","wp10","wp11","wp12","wp13","wp14","wp15","wp16","wp17","wp18","wp19","wp20","wp21","wp22","wp23","wp24"]
  dat2 = loadtxt(filename2)
  fsat2 = dat2[:,1]
  Mm2 = dat2[:,2]
  fsat2 = fsat2[dropouts:len(dat2)]
  bfsat2_err2 = median_error(fsat2)
  Mm2 = Mm2[dropouts:len(dat2)]
  bMm2_err2 = median_error(Mm2)
  Mm2 = np.asarray(Mm2)
  fsat2 = np.asarray(fsat2)
  chi2_2 = dat2[:,4]  
  chi2_2 = chi2_2[dropouts:len(dat2)]
  chi2_2 = np.asarray(chi2_2)
  bferr2=median_error(fsat2)
  merr2= median_error(Mm2)   
  w2=np.where(chi2_2 == min(chi2_2))
  bMm2 = np.mean(Mm2)
  bfsat2 = np.mean(fsat2)

  cols3=["n","fsat1", "Mm1","dm","chi2","wp0","wp1","wp2","wp3","wp4","wp5","wp6","wp7","wp8","wp9","wp10","wp11","wp12","wp13","wp14","wp15","wp16","wp17","wp18"]
  dat3 = loadtxt(filename3)
  fsat3 = dat3[:,1]
  Mm3 = dat3[:,2]
  fsat3 = fsat3[dropouts:len(dat3)]
  bfsat3_err3 = median_error(fsat3)
  Mm3 = Mm3[dropouts:len(dat3)]
  bMm3_err3 = median_error(Mm3)
  Mm3 = np.asarray(Mm3)
  fsat3 = np.asarray(fsat3)
  chi2_3 = dat3[:,4]  
  chi2_3 = chi2_3[dropouts:len(dat3)]
  chi2_3 = np.asarray(chi2_3)
  bferr3=median_error(fsat3)
  merr3= median_error(Mm3)   
  w3=np.where(chi2_3 == min(chi2_3))
  bMm3 = np.mean(Mm3)
  bfsat3 = np.mean(fsat3)

  cols4=["n","fsat1", "Mm1","dm","chi2","wp0","wp1","wp2","wp3","wp4","wp5","wp6","wp7","wp8","wp9","wp10","wp11","wp12","wp13","wp14","wp15","wp16","wp17","wp18"]
  dat4 = loadtxt(filename4)
  fsat4 = dat4[:,1]
  Mm4 = dat4[:,2]
  fsat4 = fsat4[dropouts:len(dat4)]
  bfsat4_err4 = median_error(fsat4)
  Mm4 = Mm4[dropouts:len(dat4)]
  bMm4_err4 = median_error(Mm4)
  Mm4 = np.asarray(Mm4)
  fsat4 = np.asarray(fsat4)
  chi2_4 = dat4[:,4]  
  chi2_4 = chi2_4[dropouts:len(dat4)]
  chi2_4 = np.asarray(chi2_4)
  bferr4=median_error(fsat4)
  merr4= median_error(Mm4)   
  w4=np.where(chi2_4 == min(chi2_4))
  bMm4 = np.mean(Mm4)
  bfsat4 = np.mean(fsat4)

  # Reporting .....
  pars = 2
  print('  ********************************  ')
  print('This chain has reached to a chisq of ', min(chi2_1), '(reduced chisq of ',min(chi2_1)/(len(cols1)-5-pars-1), ') for dof= ',(len(cols1)-5-pars-1),':')
  print('Mm1 = ',Mm1[w1][0]/1e12, 'e12 M_sun/h  ')
  print('fsat1 = ', fsat1[w1][0])
  print('mean fsat1= ',np.mean(fsat1),'+',bferr1[1],'-',bferr1[2])
  print('mean Mm1= ' ,np.mean(Mm1)/1e12,' e12+',merr1[1]/1e12,' e12-',merr1[2]/1e12)
  print('mean chi2= ',np.mean(chi2_1))

  print('  ********************************  ')
  print('This chain has reached to a chisq of ', min(chi2_2), '(reduced chisq of ',min(chi2_2)/(len(cols2)-5-pars-1), ') for dof= ',(len(cols2)-5-pars-1),':')
  print('Mm2 = ',Mm2[w2][0]/1e12, 'e12 M_sun/h  ')
  print('fsat2 = ', fsat2[w2][0])
  print('mean fsat2= ',np.mean(fsat2),'+',bferr2[1],'-',bferr2[2])
  print('mean Mm2= ' ,np.mean(Mm2)/1e12,' e12+',merr2[1]/1e12,' e12-',merr2[2]/1e12)
  print('mean chi2= ',np.mean(chi2_2))

  print('  ********************************  ')
  print('This chain has reached to a chisq of ', min(chi2_3), '(reduced chisq of ',min(chi2_3)/(len(cols3)-5-pars-1), ') for dof= ',(len(cols3)-5-pars-1),':')
  print('Mm3 = ',Mm3[w3][0]/1e12, 'e12 M_sun/h  ')
  print('fsat3 = ', fsat3[w3][0])
  print('mean fsat3= ',np.mean(fsat3),'+',bferr3[1],'-',bferr3[2])
  print('mean Mm3= ' ,np.mean(Mm3)/1e12,' e12+',merr3[1]/1e12,' e12-',merr3[2]/1e12)
  print('mean chi2= ',np.mean(chi2_3))

  print('  ********************************  ')
  print('This chain has reached to a chisq of ', min(chi2_4), '(reduced chisq of ',min(chi2_4)/(len(cols4)-5-pars-1), ') for dof= ',(len(cols4)-5-pars-1),':')
  print('Mm4 = ',Mm4[w4][0]/1e12, 'e12 M_sun/h  ')
  print('fsat4 = ', fsat4[w4][0])
  print('mean fsat4= ',np.mean(fsat4),'+',bferr4[1],'-',bferr4[2])
  print('mean Mm4= ' ,np.mean(Mm4)/1e12,' e12+',merr4[1]/1e12,' e12-',merr4[2]/1e12)
  print('mean chi2= ',np.mean(chi2_4))
  print('  ********************************  ')


  
  
  #--------------------------------- plotting 
  
  nmod = 35
  rparr = np.logspace(-1.95, log10(100),num=nmod, endpoint=True )
      
  rp = [1.88000000e-02,   2.28000000e-02,   2.76000000e-02,3.34000000e-02, 0.034 ,  0.06069197,  0.10833867,  0.1933908 ,  0.34521376,0.61622653,  1.1, 1.10100000e+00,   1.91300000e+00,2.76600000e+00,   3.99800000e+00,   5.77900000e+00,8.35300000e+00,   1.20740000e+01,   1.74520000e+01,2.52250000e+01,   3.64620000e+01,   4.38370000e+01,5.27030000e+01,   6.33630000e+01,   7.61800000e+01, 9.15880000e+01]
      
  #wparr= wpfunc(1.55,6.46178e-6,rparr,np.mean(Mm1),np.mean(fsat1),np.mean(delm), 0.01, 1000, 0.5, B=0)      
      
  # SE: creating the envelope --------------------------------------------------
  wup1 = []
  wlo1 = []
  wup2 = []
  wlo2 = []
  wcen = []
  
      
  #ALL FOUR KDE POINTS-----------------------------------------------------------------------------------
  nwp=27  # how many total number of wp data points exits in the data file that had been used for the fit
  wKDE=3  # what element of the wp array is the KDE data point(s)?
  weBOSS = 12 #where the eBOSS stars
  wst,wen,wst2,wen2 = 0,3,7,12  #where the KO12 set starts and ends - before and after the KDE points
  #-------------------------------------------------------------------------------------------------------
    
  for i in range(nwp):
          
          w= dat1[:,i+5]
          en1 = onesig_envelope(w)
          #en2 = twosig_envelope(w)          
          wup1.append(en1[2])
          wcen.append(en1[1])
          wlo1.append(en1[0])
          #print('wup1= ', wup1[i],'wlo1= ', wlo1[i], 'wcen=',wcen[i])
          #wup2.append(en2[2])
          #wlo2.append(en2[0])
          
  wlo1 = np.asarray(wlo1)   
  wup1 = np.asarray(wup1) 
  #wlo2 = np.asarray(wlo2)   
  #wup2 = np.asarray(wup2)   

  wcen = np.asarray(wcen) 
#-------------------------------------------------------------------      
  rparr_kpc = np.logspace(-1.95, log10(2.0),num=10, endpoint=True )
      
  rparr_mpc = np.logspace(0.301, log10(100),num=25, endpoint=True )
      
  wp_kpc= wpfunc(1.55,6.46178e-6,rparr_kpc,np.mean(Mm1),np.mean(fsat1),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_up= wpfunc(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm1)+merr1[1]),(np.mean(fsat1)+bferr1[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_lo= wpfunc(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm1)-merr1[2]),(np.mean(fsat1)-bferr1[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc= wpfunc(1.55,6.46178e-6,rparr_mpc,np.mean(Mm1),np.mean(fsat1),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc_up= wpfunc(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm1)+merr1[1]),(np.mean(fsat1)+bferr1[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)

  wp_mpc_lo= wpfunc(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm1)-merr1[2]),(np.mean(fsat1)-bferr1[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)


   
     
  obs = loadtxt('FIT1_rp_wp_err_8KO12+4KDE+eBOSS.dat')
  rpd = obs[:,0]
  wpd = obs[:,1]
  err = obs[:,2]
    
  #------------------------------------------Wp best fit PLOTTING ----------------------------------------
  fig=plt.figure()
  ax= fig.add_subplot(111)
      
  # axes limits
  x1=0.01; x2=90
  y1=0.2; y2=3.e5

  # Genral axes
  ax.set_xlim(x1, x2)
  ax.set_ylim(y1, y2)
  ax.minorticks_on()


  # additional Y-axis (on the right)
  y_ax = ax.twinx()
  y_ax.set_yscale('log')
  y_ax.set_ylim(y1, y2)
  y_ax.set_yticklabels([])
  y_ax.minorticks_on()

  # additional X-axis (on the top)
  x_ax = ax.twiny()
  x_ax.set_xscale('log')
  x_ax.set_xlim(x1, x2)
  x_ax.set_xticklabels([])
  x_ax.minorticks_on()

  ax.set_xscale('log')
  ax.set_yscale('log') 
  ax.set_ylabel(r'$\rm w_p[h^{-1}Mpc]$', fontsize=12)
  ax.set_xlabel(r'$\rm r_p[h^{-1}Mpc]$', fontsize=12)
      
  # Grid
  ax.grid(False)

  # Plots    
  p1,= ax.plot(rparr_kpc,wp_kpc,'b-', label=r'$\rm KDE+KO12+eBOSS$')

  ax.fill_between(rparr_kpc,wp_kpc_up,wp_kpc_lo,alpha=0.6, edgecolor='#00EEEE', facecolor='#00EEEE') #for kpc 
  ax.fill_between(rparr_mpc,wp_mpc_up,wp_mpc_lo,alpha=0.6, edgecolor='#00EEEE', facecolor='#00EEEE') #for Mpc   
  
  ko = ax.errorbar(rpd[wst2:wen2],wpd[wst2:wen2],yerr=err[wst2:wen2],linestyle='None',marker='s',markersize=4.,color='purple',markerfacecolor='white',elinewidth=0.6,ecolor='purple',capsize=2.0,capthick=1.,label='KO12') 
 

  kdeb = ax.errorbar(rpd[wen2:len(cols1)-1],wpd[wen2:len(cols1)-1],yerr=err[wen2:len(cols1)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.6,ecolor='black',capsize=2.0,capthick=1.,label='This work: KDE+eBOSS') 
  ax.errorbar(rpd[wen:wst2],wpd[wen:wst2],yerr=err[wen:wst2],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.6,ecolor='black',capsize=2.0,capthick=1.) 

  p12,= ax.plot(rparr_mpc,wp_mpc,'b-')#,color="blue",label=r'$\rm Best fit$')



  #------------------------------------------- SECOND BEST FIT  ---------------------------------------------

  #ALL KO12+eBOSS POINTS-----------------------------------------------------------------------------------
  nwp=25  # how many total number of wp data points exits in the data file that had been used for the fit
  weBOSS = 10 #where the eBOSS stars
  wst,wen = 0,10 #where the KO12 set starts and ends - before and after the KDE points
  #-------------------------------------------------------------------------------------------------------
  wup1 = []
  wlo1 = []
  wcen = []
  
  for i in range(nwp):
          
          w= dat2[:,i+5]
          en1 = onesig_envelope(w)
          wup1.append(en1[2])
          wcen.append(en1[1])
          wlo1.append(en1[0])

  wlo1 = np.asarray(wlo1)   
  wup1 = np.asarray(wup1)    
  wcen = np.asarray(wcen) 
  
  
  rparr_kpc = np.logspace(-1.95, log10(2.0),num=10, endpoint=True )
      
  rparr_mpc = np.logspace(0.301, log10(100),num=25, endpoint=True )
      
  wp_kpc= wpfunc(1.55,6.46178e-6,rparr_kpc,np.mean(Mm2),np.mean(fsat2),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_up= wpfunc(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm2)+merr2[1]),(np.mean(fsat2)+bferr2[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_lo= wpfunc(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm2)-merr2[2]),(np.mean(fsat2)-bferr2[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc= wpfunc(1.55,6.46178e-6,rparr_mpc,np.mean(Mm2),np.mean(fsat2),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc_up= wpfunc(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm2)+merr2[1]),(np.mean(fsat2)+bferr2[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)

  wp_mpc_lo= wpfunc(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm2)-merr2[2]),(np.mean(fsat2)-bferr2[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)
   
     
  obs = loadtxt('FIT2_rp_wp_err_10KO12+eBOSS.dat')
  rpd = obs[:,0]
  wpd = obs[:,1]
  err = obs[:,2]

  p2,= ax.plot(rparr_kpc,wp_kpc,'g-', label=r'$\rm KO12+eBOSS$')

  ax.fill_between(rparr_kpc,wp_kpc_up,wp_kpc_lo,alpha=0.6, edgecolor='#ffff00', facecolor='#ffff00') #for kpc 
  ax.fill_between(rparr_mpc,wp_mpc_up,wp_mpc_lo,alpha=0.6, edgecolor='#ffff00', facecolor='#ffff00') #for Mpc   
  
  ko = ax.errorbar(rpd[wst:wen],wpd[wst:wen],yerr=err[wst:wen],linestyle='None',marker='s',markersize=4.,color='purple',markerfacecolor='white',elinewidth=0.6,ecolor='purple',capsize=2.0,capthick=1.,label='KO12') 
  
  kdeb = ax.errorbar(rpd[wen:len(cols2)-1],wpd[wen:len(cols2)-1],yerr=err[wen:len(cols2)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.,label='This work: KDE+eBOSS') 
  ax.errorbar(rpd[wen:len(cols2)-1],wpd[wen:len(cols2)-1],yerr=err[wen:len(cols2)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.) 

  p22,= ax.plot(rparr_mpc,wp_mpc,'g-')#,color="blue",label=r'$\rm Best fit$')

  #-------------------------------------------------------------------------------------------
  #------------------------------------------- THIRD BEST FIT  ---------------------------------------------

  #ALL KO12+eBOSS POINTS-----------------------------------------------------------------------------------
  nwp=19  # how many total number of wp data points exits in the data file that had been used for the fit
  weBOSS = 4 #where the eBOSS stars
  wst,wen = 0,4 #where the KO12 set starts and ends - before and after the KDE points
  #-------------------------------------------------------------------------------------------------------
  wup1 = []
  wlo1 = []
  wcen = []
  
  for i in range(nwp):
          
          w= dat3[:,i+5]
          en1 = onesig_envelope(w)
          wup1.append(en1[2])
          wcen.append(en1[1])
          wlo1.append(en1[0])

  wlo1 = np.asarray(wlo1)   
  wup1 = np.asarray(wup1)    
  wcen = np.asarray(wcen) 
  
  
  rparr_kpc = np.logspace(-1.95, log10(2.0),num=10, endpoint=True )
      
  rparr_mpc = np.logspace(0.301, log10(100),num=25, endpoint=True )
      
  wp_kpc= wpfunc(1.55,6.46178e-6,rparr_kpc,np.mean(Mm3),np.mean(fsat3),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_up= wpfunc(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm3)+merr3[1]),(np.mean(fsat3)+bferr3[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_lo= wpfunc(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm3)-merr3[2]),(np.mean(fsat3)-bferr3[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc= wpfunc(1.55,6.46178e-6,rparr_mpc,np.mean(Mm3),np.mean(fsat3),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc_up= wpfunc(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm3)+merr3[1]),(np.mean(fsat3)+bferr3[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)

  wp_mpc_lo= wpfunc(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm3)-merr3[2]),(np.mean(fsat3)-bferr3[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)
   
     
  obs = loadtxt('FIT3_rp_wp_err_4KDE+eBOSS.dat')
  rpd = obs[:,0]
  wpd = obs[:,1]
  err = obs[:,2]

  p3,= ax.plot(rparr_kpc,wp_kpc,'r-', label=r'$\rm KDE+eBOSS$')

  ax.fill_between(rparr_kpc,wp_kpc_up,wp_kpc_lo,alpha=0.6, edgecolor='#ff00ff', facecolor='#ff00ff') #for kpc 
  ax.fill_between(rparr_mpc,wp_mpc_up,wp_mpc_lo,alpha=0.6, edgecolor='#ff00ff', facecolor='#ff00ff') #for Mpc   
  
  kde = ax.errorbar(rpd[wst:wen],wpd[wst:wen],yerr=err[wst:wen],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='k',capsize=2.0,capthick=1.,label='KDE') 
  
  kdeb = ax.errorbar(rpd[wen:len(cols3)-1],wpd[wen:len(cols3)-1],yerr=err[wen:len(cols3)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.,label='This work: KDE+eBOSS') 
  ax.errorbar(rpd[wen:len(cols3)-1],wpd[wen:len(cols3)-1],yerr=err[wen:len(cols3)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.) 

  p23,= ax.plot(rparr_mpc,wp_mpc,'r-')#,color="blue",label=r'$\rm Best fit$')

  #-------------------------------------------------------------------------------------------
   #------------------------------------------- Fourth BEST FIT  ---------------------------------------------

  #ALL KO12+eBOSS POINTS 10c -----------------------------------------------------------------------------------
  nwp=19  # how many total number of wp data points exits in the data file that had been used for the fit
  weBOSS = 4 #where the eBOSS stars
  wst,wen = 0,4 #where the KO12 set starts and ends - before and after the KDE points
  #-------------------------------------------------------------------------------------------------------
  wup1 = []
  wlo1 = []
  wcen = []
  
  for i in range(nwp):
          
          w= dat4[:,i+5]
          en1 = onesig_envelope(w)
          wup1.append(en1[2])
          wcen.append(en1[1])
          wlo1.append(en1[0])

  wlo1 = np.asarray(wlo1)   
  wup1 = np.asarray(wup1)    
  wcen = np.asarray(wcen) 
  
  
  rparr_kpc = np.logspace(-1.95, log10(2.0),num=10, endpoint=True )
      
  rparr_mpc = np.logspace(0.301, log10(100),num=25, endpoint=True )
      
  wp_kpc= wpfunc10c(1.55,6.46178e-6,rparr_kpc,np.mean(Mm4),np.mean(fsat4),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_up= wpfunc10c(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm4)+merr4[1]),(np.mean(fsat4)+bferr4[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_kpc_lo= wpfunc10c(1.55,6.46178e-6,rparr_kpc,(np.mean(Mm4)-merr4[2]),(np.mean(fsat4)-bferr4[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc= wpfunc10c(1.55,6.46178e-6,rparr_mpc,np.mean(Mm4),np.mean(fsat4),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
  wp_mpc_up= wpfunc10c(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm4)+merr4[1]),(np.mean(fsat4)+bferr4[1]),np.mean(delm), 0.01, 1000, 0.5, B=0)

  wp_mpc_lo= wpfunc10c(1.55,6.46178e-6,rparr_mpc,(np.mean(Mm4)-merr4[2]),(np.mean(fsat4)-bferr4[2]),np.mean(delm), 0.01, 1000, 0.5, B=0)
   
     
  obs = loadtxt('FIT3_rp_wp_err_4KDE+eBOSS.dat')
  rpd = obs[:,0]
  wpd = obs[:,1]
  err = obs[:,2]

  p4,= ax.plot(rparr_kpc,wp_kpc,'-',color='orange', label=r'$\rm (c=10c)\, KDE+eBOSS$')

  ax.fill_between(rparr_kpc,wp_kpc_up,wp_kpc_lo,alpha=0.6, edgecolor='#CC4F1B', facecolor='#CC4F1B') #for kpc 
  ax.fill_between(rparr_mpc,wp_mpc_up,wp_mpc_lo,alpha=0.6, edgecolor='#CC4F1B', facecolor='#CC4F1B') #for Mpc   
  
  kde = ax.errorbar(rpd[wst:wen],wpd[wst:wen],yerr=err[wst:wen],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='k',capsize=2.0,capthick=1.,label='KDE') 
  
  kdeb = ax.errorbar(rpd[wen:len(cols4)-1],wpd[wen:len(cols4)-1],yerr=err[wen:len(cols4)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.,label='This work: KDE+eBOSS') 
  ax.errorbar(rpd[wen:len(cols4)-1],wpd[wen:len(cols4)-1],yerr=err[wen:len(cols4)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.) 

  p24,= ax.plot(rparr_mpc,wp_mpc,'-',color='orange')#,color="blue",label=r'$\rm Best fit$')

  #-------------------------------------------------------------------------------------------

  lns = [p12, p1, p22, p2, p23, p3, p24, p4]
  ax.legend(handles=lns, loc='best')
  fig.savefig('HOD_4BESTFITS.eps',bbox_inches='tight')

  plt.show()
      
  sys.exit()

def main():
    
       stats(73,3000,2,plots=True)

       plt.show()
       
if __name__ == '__main__':
  main()
      
      
