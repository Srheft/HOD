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
from wp_mcmc_2pars import median_error
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
def meanNM(M,delm,Mm):
     
    return (1./(sqrt(2.*pi)*delm))*exp(-(ln(M/Mm)*ln(M/Mm))/(2.*delm*delm))
#============================================================================

def stats(seed,nchain,pars,plots=False):
  
  filename ='chain2_2pars_new_ONLY1KDE_ALL.dat'#'chain2_2pars_newdatapoints_ALL.dat' #'chain2_2pars_Jul18.txt'#24 wp #'results_seed73_nchain3000_3pars.txt' #'results_seed'+str(seed)+'_nchain'+str(nchain)+'_'+str(pars)+'pars.txt'
  print('Working with ', filename)
  dropouts = 0#40
  cols=["n","fsat", "Mm","dm","chi2","wp0","wp1","wp2","wp3","wp4","wp5","wp6","wp7","wp8","wp9","wp10","wp11","wp12","wp13","wp14","wp15","wp16","wp17","wp18","wp19","wp20","wp21","wp22"]#,"wp23","wp24"]
  dat = loadtxt(filename)#,dtype=<type 'float'>, comments='#',delimiter=',')#wpobs_err_4mcmc.dat
  naccept = dat[:,0]
  fsat = dat[:,1]
  Mm = dat[:,2]
  
  fsat = fsat[dropouts:len(dat)]
  bfsat_err = median_error(fsat)
  
  Mm = Mm[dropouts:len(dat)]
  bMm_err = median_error(Mm)
  Mm = np.asarray(Mm)
  fsat = np.asarray(fsat)
  
  
  delm = dat[:,3]
  delm = delm[dropouts:len(dat)]
  bdelm_err = median_error(delm)
  bdm_uerr = bdelm_err[1]
  bdm_lerr = bdelm_err[2]
  chi2 = dat[:,4]
  delm = np.asarray(delm)
  
  chi2 = chi2[dropouts:len(dat)]
  chi2 = np.asarray(chi2)

  bferr=median_error(fsat)
  print( 'mean fsat= ',np.mean(fsat),'+',bferr[1],'-',bferr[2])
    
  merr= median_error(Mm)
  print('mean Mm= ' ,np.mean(Mm)/1e12,' e12+',merr[1]/1e12,' e12-',merr[2]/1e12)
  
  print('mean delm= ',np.mean(delm),'+',bdm_uerr, '-',bdm_lerr)
  bdelm = np.mean(delm)
   
            
  print('mean chi2= ',np.mean(chi2))
  
  w=np.where(chi2 == min(chi2))
  bMm = np.mean(Mm)
  bfsat = np.mean(fsat)

  
  print('  ********************************  ')
  print('This chain has reached to a chisq of ', min(chi2), '(reduced chisq of ',min(chi2)/(len(cols)-pars-1), ') for:')
  print('Mm = ',Mm[w][0]/1e12, 'e12 M_sun/h  ')
  print('fsat = ', fsat[w][0])
  
  if (pars ==3): 
      print('delm =', delm[w][0]) 
  
  print('  ********************************')
  
  
  if plots:
      
      nmod = 35
      rparr = np.logspace(-1.95, log10(100),num=nmod, endpoint=True )
      
      rp = [1.88000000e-02,   2.28000000e-02,   2.76000000e-02,3.34000000e-02, 0.034 ,  0.06069197,  0.10833867,  0.1933908 ,  0.34521376,0.61622653,  1.1, 1.10100000e+00,   1.91300000e+00,2.76600000e+00,   3.99800000e+00,   5.77900000e+00,8.35300000e+00,   1.20740000e+01,   1.74520000e+01,2.52250000e+01,   3.64620000e+01,   4.38370000e+01,5.27030000e+01,   6.33630000e+01,   7.61800000e+01, 9.15880000e+01]
      
      wparr= wpfunc(1.55,6.46178e-6,rparr,np.mean(Mm),np.mean(fsat),np.mean(delm), 0.01, 1000, 0.5, B=0)      
      
      # SE: creating the envelope --------------------------------------------------
      wup1 = []
      wlo1 = []
      wup2 = []
      wlo2 = []
      wcen = []
      
      ##ONLY ONE KDE POINT-----------------------------------------------------------------------------------
      nwp=23  # how many total number of wp data points exits in the data file that had been used for the fit
      wKDE=4  # what element of the wp array is the KDE data point(s)?
      weBOSS = 9 #where the eBOSS stars
      wst,wen,wst2,wen2 = 0,4,5,9  #where the KO12 set starts and ends - before and after the KDE points
      
      #ALL FOUR KDE POINTS-----------------------------------------------------------------------------------
      #nwp=25  # how many total number of wp data points exits in the data file that had been used for the fit
      #wKDE=3  # what element of the wp array is the KDE data point(s)?
      #weBOSS = 12 #where the eBOSS stars
      #wst,wen,wst2,wen2 = 0,3,7,11  #where the KO12 set starts and ends - before and after the KDE points
     #-------------------------------------------------------------------------------------------------------
      
      
      for i in range(nwp):
          
          w= dat[:,i+5]
          en1 = onesig_envelope(w)
          en2 = twosig_envelope(w)          
          wup1.append(en1[2])
          wcen.append(en1[1])
          wlo1.append(en1[0])
          
          wup2.append(en2[2])
          wlo2.append(en2[0])
          
      wlo1 = np.asarray(wlo1)   
      wup1 = np.asarray(wup1) 
      wlo2 = np.asarray(wlo2)   
      wup2 = np.asarray(wup2)   

      wcen = np.asarray(wcen) 
      
      
      
      print 'wup1= ', wup1
      print 'wlo1= ',wlo1
      print 'wcen= ',wcen
      print 'wup2= ', wup2
      print 'wlo2= ',wlo2
      
      #sys.exit()
      #--------------------------------------------------  PLOTTING ---------------------------------------------------------  
      rparr_10c = np.logspace(-1.95, log10(2.0),num=10, endpoint=True )
      
      rparr_mpc = np.logspace(2.0, log10(100),num=10, endpoint=True )
      
      wp_10c = wpfunc10c(1.55,6.46178e-6,rparr_10c,np.mean(Mm),np.mean(fsat),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
      wp_10c_mpc = wpfunc10c(1.55,6.46178e-6,rparr_mpc,np.mean(Mm),np.mean(fsat),np.mean(delm), 0.01, 1000, 0.5, B=0)

      wparr3= wpfunc(1.55,6.46178e-6,rparr,np.mean(Mm),0.01,np.mean(delm), 0.01, 1000, 0.5, B=0)
      
      wparr4= wpfunc(1.55,6.46178e-6,rparr,np.mean(Mm),0.9,np.mean(delm), 0.01, 1000, 0.5, B=0)
      
      #wp_10c_3pars = wpfunc10c(1.55,6.46178e-6,rparr,np.mean(Mm),np.mean(fsat),np.mean(delm), 0.01, 1000, 0.5, B=0)
      
      dat2 = loadtxt('chain2_2pars_new_ONLY1KDE_ALL.dat')#loadtxt('results_seed73_nchain2000_2pars.txt')#'results_seed'+str(seed)+'_nchain'+str(nchain)+'_2pars.txt'
      cols=["n","fsat", "Mm","dm","chi2","wp0","wp1","wp2","wp3","wp4","wp5","wp6","wp7","wp8","wp9","wp10","wp11","wp12","wp13","wp14","wp15","wp16","wp17","wp18","wp19","wp20","wp21","wp22"]#,"wp23","wp24"]
      naccept2 = dat2[:,0]
      fsat2 = dat2[:,1]
      Mm2 = dat2[:,2]
      chi22 = dat2[:,3]
      fsat2 = fsat[dropouts:len(dat)]
      Mm2 = Mm[dropouts:len(dat)]

      #wparr2= wpfunc(1.55,6.46178e-6,rparr,np.mean(Mm2),np.mean(fsat2),0.75, 0.01, 1000, 0.5, B=0) # fitting parameters are Mm,fsat,delm

      
      #wp_10c_2pars_highn = wpfunc10c(1.55,6.46178e-6,rparr,np.mean(Mm2),np.mean(fsat2),0.75, 0.01, 1000, 0.5, B=0)
      
      obs = loadtxt('rp_wp_err_ONLY1KDE.dat')#'rp_wp_err_NOTALL_KO12.dat')#
      rpd = obs[:,0]
      wpd = obs[:,1]
      err = obs[:,2]
    
    
    #------------------------------------------Wp best fit PLOTTING ----------------------------------------
      fig=plt.figure()
      ax= fig.add_subplot(111)
      
      # axes limits
      x1=0.01; x2=100
      y1=0.1; y2=3.e5

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
      print(len(rpd[0:len(cols)-1]),len(wup1),len(wlo1))
      
      #plots
      p1,= ax.plot(rparr,wparr,'b-', label="Best fit")
      #ax.fill_between(rpd,wup1,wlo1,alpha=0.6, edgecolor='#00EEEE', facecolor='#00EEEE')
      
      
      #ax.fill_between(rpd[0:len(cols)-1],wup2,wlo2,alpha=0.6, edgecolor='#CC4F1B', facecolor='#FF9848')

      d1= ax.errorbar(rpd[wst:wen],wpd[wst:wen],yerr=err[wst:wen],linestyle='None',marker='s',markersize=5.,color='magenta',fillstyle='full',elinewidth=0.5,ecolor='magenta',capsize=2.0,capthick=1.,label='Kayo & Oguri 2012')  
      ax.errorbar(rpd[wst2:wen2],wpd[wst2:wen2],yerr=err[wst2:wen2],linestyle='None',marker='s',markersize=5.,color='magenta',fillstyle='full',elinewidth=0.5,ecolor='magenta',capsize=2.0,capthick=1.)  

      d2=ax.errorbar(rpd[wen2:len(cols)-1],wpd[wen2:len(cols)-1],yerr=err[wen2:len(cols)-1],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.,label='This work: KDE+eBOSS') 
      ax.errorbar(rpd[wKDE:wst2],wpd[wKDE:wst2],yerr=err[wKDE:wst2],linestyle='None',marker='o',markersize=3.,color='k',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.) 

      #plt.plot(rparr[0:29],wparr2[0:29],'k-.')

      #p2,= ax.plot(rparr[0:29],wp_10c_3pars[0:29],'g--',label=r"$10 x \bar c")
      
      #p2,= ax.plot(rparr,wparr4,'g-.',label=r"$f_{sat}=1.0$")
      p2,= ax.plot(rparr_10c,wp_10c,'--',color="orange",label=r'$\rm 10 \times \bar c$')
      #p3,= ax.plot(rparr_10c,wp_20c,'-.',color="cyan",label=r'$\rm 20 \times \bar c$')
 
      #p3,= ax.plot(rparr,wparr3,'r--',label=r"$f_{sat}=0.01$")


      lns = [ p2,p1]# [p1, p2, p3, d1, d2]
      ax.legend(handles=lns, loc='best')
      fig.savefig('Bestfit_wpmodel1KDE_2pars_10c.eps',bbox_inches='tight')


      plt.show()
      
      sys.exit()

     ###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      cols=["n","fsat", "Mm","dm","q","wp0","wp1","wp2","wp3","wp4","wp5","wp6","wp7","wp8","wp9","wp10","wp11","wp12","wp13","wp14","wp15","wp16","wp17","wp18","wp19","wp20","wp21","wp22"]#,"wp23","wp24"]
      df = pd.DataFrame(dat, columns=cols)
      ###g= sns.jointplot(x="Mm", y="fsat", data=df,space = 0.,stat_func=None,marker='.',joint_kws={'color':'green'})#,xlabel=r'$\rm M[h^{-1}M_{\odot}]$',ylabel=r'$f_{sat}$',r"$\Delta_m$"
      ###g.set_axis_labels(r'$\rm M_{m}[h^{-1}M_{\odot}]$', r"$f_{sat}$")
      ###g.savefig('scatter_Mm_fsat.eps', bbox_inches='tight')

      ###plt.show()
     ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ###gg=sns.jointplot(x="Mm", y="dm", data=df,space = 0.,stat_func=None,joint_kws=None, marginal_kws=dict(bins=15, rug=False),
                   ###annot_kws=dict(stat="r"),
                   ###s=40, edgecolor="w", linewidth=1)
      ###gg.set_axis_labels(r'$\rm M_{m}[h^{-1}M_{\odot}]$', r'$\Delta_m $', fontsize=20)
      ###gg.savefig('scatter_Mm_dm.eps', bbox_inches='tight')

      ###plt.show() 
     ###+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      #g = sns.JointGrid(x="dm", y="fsat", data=df, space=0)
      #g.plot_joint(sns.kdeplot, cmap="Purples_d",n_levels=6)
      #g.plot_marginals(sns.kdeplot, shade=True)
      #g.set_axis_labels(r'$\Delta_m $', r"$f_{sat}$", fontsize=12)
      #g.savefig('contour_dm_fsat.eps')#, bbox_inches='tight'
      #plt.show()
     ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      #g = (sns.jointplot(x="Mm", y="fsat",data=df,space = 0.,stat_func=None,color="grey",marker='.', s=2).plot_joint(sns.kdeplot, zorder=0, n_levels=7))
      #g.set_axis_labels(r'$\rm M_{m}[h^{-1}M_{\odot}]$', r"$f_{sat}$", fontsize=12)
      #g.savefig('contour_Mm_fsat.eps')#, bbox_inches='tight'
      #plt.show()
     ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      #g = sns.jointplot(x="Mm", y="fsat", data=df, kind="kde", color="m",space = 0.,stat_func=None)
      #g.plot_joint(plt.scatter, c="w", s=1, linewidth=1, marker=".")
      #g.ax_joint.collections[0].set_alpha(0)
      #g.set_axis_labels(r'$\rm M_{m}[h^{-1}M_{\odot}]$', r"$f_{sat}$", fontsize=12)
      #g.savefig('contour_scatter_Mm_fsat.eps')#, bbox_inches='tight'
      #plt.show()
     
     #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      #g = sns.JointGrid(x="dm", y="fsat", data=df,space = 0.)#,  xlim=(0, 50), ylim=(0, 8)
      #g = g.plot_joint(sns.kdeplot, cmap="Purples_d",space = 0.)
      #g = g.plot_marginals(sns.kdeplot, color="m", shade=True)
      #g.set_axis_labels(r'$\Delta_m $', r"$f_{sat}$", fontsize=12)
      #plt.show()
     ##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      #g = sns.jointplot(x="dm", y="fsat",data=df, color="grey", marker=".",s=1.,stat_func=None,space=0).plot_joint(sns.kdeplot, zorder=0, n_levels=6)
      #g.set_axis_labels(r'$\Delta_m $', r"$f_{sat}$", fontsize=12)
      #g.savefig('contour_dm_fsat.eps')#, bbox_inches='tight')
      #plt.show()
     
      #g = sns.jointplot(x="Mm", y="dm",data=df, color="grey", marker=".",s=1.,stat_func=None,space=0).plot_joint(sns.kdeplot, zorder=0, n_levels=6)
      #g.set_axis_labels(r"$M_{m}$", r'$\Delta_m $', fontsize=12)
      #g.savefig('contour_Mm_dm.eps')#, bbox_inches='tight')
      #plt.show()

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
      M = np.logspace(11,14,num=100,endpoint=True) 
      hmf= loadtxt('mVector_PLANCK-SMTz1.55.txt')
      Mr = hmf[:,0]
      dnm = hmf[:,5]
      f = interp1d(Mr, dnm)
    
      Ndist = []
      Nmorig = []
      Nother = []
      hostDMH = []
 
      if (pars == 3):
            #seed,nchain = 10, 1000
            dat2 = loadtxt('chain2_2pars_new_ONLY1KDE_ALL.dat')#'chain2_2pars_newdatapoints_ALL.dat')#'chain2_2pars_Jul18.txt')#'results_seed73_nchain2000_2pars.txt')#'results_seed'+str(seed)+'_nchain'+str(nchain)+'_2pars.txt')
            naccept = dat2[:,0]
            fsat = dat2[:,1]
            Mm = dat2[:,2]
            #chi2 = dat2[:,3]
            fsat = fsat[dropouts:len(dat)]
            Mm = Mm[dropouts:len(dat)]
            #chi2 = chi2[dropouts:len(dat)]
            #w=np.where(chi2 == min(chi2))
            Mm2par = np.mean(Mm)
            fsat2par = np.mean(fsat)
  
            for mass in M:
                Nother.append(fsat2par*meanNM(mass,0.75,Mm2par))#fsat2par*
            
      if (pars == 2) :
            
            dat2 = loadtxt('chain2_2pars_new_ONLY1KDE_ALL.dat')#'chain2_2pars_newdatapoints_ALL.dat')#'chain2_2pars_Jul18.txt')#'results_seed73_nchain3000_3pars.txt')#'results_seed'+str(seed)+'_nchain'+str(nchain)+'_3pars.txt')
            naccept = dat2[:,0]
            fsat = dat2[:,1]
            Mm = dat2[:,2]
            delm = dat[:,3]
            delm = delm[dropouts:len(dat)]
            #chi2 = dat2[:,4]
            fsat = fsat[dropouts:len(dat)]
            Mm = Mm[dropouts:len(dat)]
            chi2 = chi2[dropouts:len(dat)]
            #w=np.where(chi2 == min(chi2))
            Mm3par = np.mean(Mm)
            fsat3par = np.mean(fsat)
            delm3par = np.mean(delm)
            
            for mass in M:
                Nother.append(fsat3par*meanNM(mass,delm3par,Mm3par))#fsat3par*
            
      for mass in M:
             
             Ndist.append(meanNM(mass,bdelm,bMm))
             Nmorig.append(meanNM(mass,0.75,1.33e13))
             hostDMH.append(f(mass)*meanNM(mass,bdelm,bMm))# <N(M)>*dn/dm  ---> host halo mass distribution
             
             
      #plt.plot(M,Ndist,'b-')
      #plt.plot(M,Nmorig,'r--')
      #plt.plot(M,Nother,'g-.')
      plt.plot(M,hostDMH,'k--.')
      plt.axis([1.e11, 1.e15, 0., 0.035]) # for fN*<N(m)>
      plt.axis([1.e11, 1.e15, 0., 0.6])
      plt.xscale('log') 
      plt.ylabel(r'$\rm \langle N(M) \rangle$', fontsize=13)
      plt.xlabel(r'$\rm M[h^{-1}M_{\odot}]$', fontsize=13)
      plt.savefig('bestfit_Nm.eps')#, bbox_inches='tight'
      
      
     
def main():
    
       stats(73,3000,2,plots=True)

       plt.show()
       
if __name__ == '__main__':
  main()
      
      
