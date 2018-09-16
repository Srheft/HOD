# -*- coding: utf-8 -*-
import sys
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
from scipy.special import spherical_jn,j0,sici # scipy.special.sici(x[, out1, out2]) 
from scipy.integrate import romberg,quad
import scipy.integrate as integrate
from random import gauss
from numpy import loadtxt, zeros,inf
import numpy as np
from time import time, sleep
from astropy.table import Table, Column
from numpy.fft import fft
from scipy import fftpack
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import warnings
from scipy  import interpolate
from HODemceeFIT_v4 import *

#============================================================================
def meanNM(M,delm,Mm):
     
    return (1./(sqrt(2.*pi)*delm))*exp(-(ln(M/Mm)*ln(M/Mm))/(2.*delm*delm))

 #==================================== Virial radius from the dark matter halo mass =============================================
def rvir(om0,delta_vir,M):
    
    rho_bg = 2.78e11*om0
    
    return (3.*M/(4*pi*delta_vir*rho_bg))**(1/3.)
 
#============================   u(k,m): is the fourier transform of the halo desity profile   =======================
#
#NOTE: from Giocoli et al 2010
def ukmcal(k,M,z,delta_vir,om0,M_star):
          
      c = cbar(M,z,M_star)
      
      r_vir= rvir(om0,delta_vir,M) 
      
      r_s = rvir(om0,delta_vir,M) / c
      
      rho_s =  M/(4.*pi*r_s**3.0*(ln(c+1.)-(c/(1.+c))))
      
      sumu =0.
      dr = 0.001
      for r in np.arange(1.e-5,r_vir,0.001):
         
          rho_rM = rho_s/((r/r_s)*(1+r/r_s)**2.)
          
          sumu += (4*pi*r**2./M)*rho_rM*sin(k*r)*dr/(k*r)
          
      return sumu
#============================   u(k,m): is the fourier transform of the halo desity profile   =======================
#
#from Giocoli et al 2010

def ukmcal10c(k,M,z,delta_vir,om0,M_star):
          
      c =10.* cbar(M,z,M_star)
      
      r_vir= rvir(om0,delta_vir,M) 
      
      r_s = rvir(om0,delta_vir,M) / c
      
      rho_s =  M/(4.*pi*r_s**3.0*(ln(c+1.)-(c/(1.+c))))
      
      sumu =0.
      
      dr = 0.001
      
      for r in np.arange(1.e-5,r_vir,0.001):
         
          rho_rM = rho_s/((r/r_s)*(1+r/r_s)**2.)
          
          sumu += (4*pi*r**2./M)*rho_rM*sin(k*r)*dr/(k*r)
          
      return sumu
  
#============================   u(k,m): is the fourier transform of the halo desity profile   =======================
#
#from Giocoli et al 2010

def ukmcal20c(k,M,z,delta_vir,om0,M_star):
          
      c =20.* cbar(M,z,M_star)
      
      r_vir= rvir(om0,delta_vir,M) 
      
      r_s = rvir(om0,delta_vir,M) / c
      
      rho_s =  M/(4.*pi*r_s**3.0*(ln(c+1.)-(c/(1.+c))))
      
      sumu =0.
      
      dr = 0.001
      
      for r in np.arange(1.e-5,r_vir,0.001):
         
          rho_rM = rho_s/((r/r_s)*(1+r/r_s)**2.)
          
          sumu += (4*pi*r**2./M)*rho_rM*sin(k*r)*dr/(k*r)
          
      return sumu
  
  
#===============================================================================
#from Giocoli et al 2010

def ukmcalc0(k,M,z,delta_vir,om0,M_star,c0):
          
      c =cbarc0(c0,M,z,M_star)
      
      r_vir= rvir(om0,delta_vir,M) 
      
      r_s = rvir(om0,delta_vir,M) / c
      
      rho_s =  M/(4.*pi*r_s**3.0*(ln(c+1.)-(c/(1.+c))))
      
      sumu =0.
      dr = 0.001
      for r in np.arange(1.e-5,r_vir,0.001):
         
          rho_rM = rho_s/((r/r_s)*(1+r/r_s)**2.)
          
          sumu += (4*pi*r**2./M)*rho_rM*sin(k*r)*dr/(k*r)
          
      return sumu

#==================================================================
#SE: from Eq 6 and 7 Chan et al 2017 (https://arxiv.org/pdf/1712.01947.pdf)
 #   
def ukmcalc0(k,M,z,delta_vir,om0,M_star,c0):
          
      c =cbarc0(c0,M,z,M_star)
      
      r_vir= rvir(om0,delta_vir,M) 
      
      r_s = rvir(om0,delta_vir,M) / c
      
      rho_s =  M/(4.*pi*r_s**3.0*(ln(c+1.)-(c/(1.+c))))
      
      sumu =0.
      dr = 0.001
      for r in np.arange(1.e-5,r_vir,0.001):
         
          rho_rM = rho_s/((r/r_s)*(1+r/r_s)**2.)
          
          sumu += (4*pi*r**2./M)*rho_rM*sin(k*r)*dr/(k*r)
          
      return sumu


#======================================================================================================
def wpfunc(z,nqobs,rparr,Mm,fsat,delm,kmin,kmax,dx, B=0):
       
    bctrl = 2.7
    warnings.filterwarnings('ignore')
    
    e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)
    
    M_star,om0,delta_vir = 2.0e13,0.308,200.#,0.01,1000.

    b,fN,nqcalc = biasint2(z,nqobs,bctrl,delm,Mm)
    
    if B > 0:
        b = B
    
    wpcalc = []
  
    qmf = 'mVector_PLANCK-SMTz'+str(z)+'.txt' # quasar mass function
    mall = loadtxt(qmf)
    M = mall[:,0]
    dnm = mall[:,5]
    
    p2k = []
    p1k = []
    dkp = 0.5
    K   = np.logspace(log10(kmin),log10(kmax),300)
  
    for k in K:#np.logspace(log10(kmin),log10(kmax)-0.5,(kmax-kmin)/dkp):
            
       k0 = kmin
            
       p2k.append(b**2.*e.P(k,n=1.)) 
       p1calc = 0.
            
       for p in range(len(M)-1):
   
                  Nm = fN*meanNM(M[p],delm,Mm)
         
                  uk0 = ukmcal(k0,M[p],z,delta_vir,om0,M_star)

                  u = ukmcal(k,M[p],z,delta_vir,om0,M_star)/uk0
                  
                  dm =  M[p+1]/M[p]
                  
                  p1calc  += M[p]*ln(10)*dm*dnm[p]*(2.*fsat*(1.-fsat)*u+fsat**2.*abs(u)**2.)*Nm**2./nqobs**2.                   
                  
       p1k.append(p1calc)

    p1k = np.asarray(p1k)
    p2k = np.asarray(p2k)
        
  
    P1K = interpolate.interp1d(K, p1k)
    P2K = interpolate.interp1d(K, p2k) 
    
    
    for r in range(len(rparr)):
      
       #dx = 0.3   
       dkp = dx/rparr[r]
      
       wpcal = 0.
       for k in np.arange(K[1], np.max(K), dkp):
            
             wpcal += k*dkp*j0(k*rparr[r])*(P1K(k)+P2K(k))/(2.*pi) 
       wpcalc.append(wpcal)
    
    #wparr = np.asarray(wpcalc)
    
    return wpcalc
#================================================================================================    
def wpfunc10c(z,nqobs,rparr,Mm,fsat,delm,kmin,kmax,dx, B=0):
       
    bctrl = 2.7
    warnings.filterwarnings('ignore')
    
    e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)
    
    M_star,om0,delta_vir = 2.0e13,0.308,200.#,0.01,1000.

    b,fN,nqcalc = biasint2(z,nqobs,bctrl,delm,Mm)
     
    if B > 0:
        b = B
    
    wpcalc = []
  
    qmf = 'mVector_PLANCK-SMTz'+str(z)+'.txt' # quasar mass function
    mall = loadtxt(qmf)
    M = mall[:,0]
    dnm = mall[:,5]
    
    p2k = []
    p1k = []
    dkp = 0.5
    K   = np.logspace(log10(kmin),log10(kmax),300)
  
    for k in K:#np.logspace(log10(kmin),log10(kmax)-0.5,(kmax-kmin)/dkp):
            
       k0 = kmin
            
       p2k.append(b**2.*e.P(k,n=1.)) 
       p1calc = 0.
            
       for p in range(len(M)-1):
   
                  Nm = fN*meanNM(M[p],delm,Mm)
         
                  uk0 = ukmcal10c(k0,M[p],z,delta_vir,om0,M_star)

                  u = ukmcal10c(k,M[p],z,delta_vir,om0,M_star)/uk0
                  
                  dm =  M[p+1]/M[p]
                  
                  p1calc  += M[p]*ln(10)*dm*dnm[p]*(2.*fsat*(1.-fsat)*u+fsat**2.*abs(u)**2.)*Nm**2./nqobs**2.                   
                  
       p1k.append(p1calc)

    p1k = np.asarray(p1k)
    p2k = np.asarray(p2k)
  
    P1K = interpolate.interp1d(K, p1k)
    P2K = interpolate.interp1d(K, p2k) 
    

    for r in range(len(rparr)):
      
       #dx = 0.3   
       dkp = dx/rparr[r]
      
       wpcal = 0.
       for k in np.arange(K[1], np.max(K), dkp):
            
             wpcal += k*dkp*j0(k*rparr[r])*(P1K(k)+P2K(k))/(2.*pi) 
       wpcalc.append(wpcal)
    
    #wparr = np.asarray(wpcalc)
    
    return wpcalc


#================================================================================================    
def wpfunc20c(z,nqobs,rparr,Mm,fsat,delm,kmin,kmax,dx, B=0):
       
    bctrl = 2.7
    warnings.filterwarnings('ignore')
    
    e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)
    
    M_star,om0,delta_vir = 2.0e13,0.308,200.#,0.01,1000.

    b,fN,nqcalc = biasint2(z,nqobs,bctrl,delm,Mm)
     
    if B > 0:
        b = B
    
    wpcalc = []
  
    qmf = 'mVector_PLANCK-SMTz'+str(z)+'.txt' # quasar mass function
    mall = loadtxt(qmf)
    M = mall[:,0]
    dnm = mall[:,5]
    
    p2k = []
    p1k = []
    dkp = 0.5
    K   = np.logspace(log10(kmin),log10(kmax),300)
  
    for k in K:#np.logspace(log10(kmin),log10(kmax)-0.5,(kmax-kmin)/dkp):
            
       k0 = kmin
            
       p2k.append(b**2.*e.P(k,n=1.)) 
       p1calc = 0.
            
       for p in range(len(M)-1):
   
                  Nm = fN*meanNM(M[p],delm,Mm)
         
                  uk0 = ukmcal20c(k0,M[p],z,delta_vir,om0,M_star)

                  u = ukmcal20c(k,M[p],z,delta_vir,om0,M_star)/uk0
                  
                  dm =  M[p+1]/M[p]
                  
                  p1calc  += M[p]*ln(10)*dm*dnm[p]*(2.*fsat*(1.-fsat)*u+fsat**2.*abs(u)**2.)*Nm**2./nqobs**2.                   
                  
       p1k.append(p1calc)

    p1k = np.asarray(p1k)
    p2k = np.asarray(p2k)
  
    P1K = interpolate.interp1d(K, p1k)
    P2K = interpolate.interp1d(K, p2k) 
    

    for r in range(len(rparr)):
      
       #dx = 0.3   
       dkp = dx/rparr[r]
      
       wpcal = 0.
       for k in np.arange(K[1], np.max(K), dkp):
            
             wpcal += k*dkp*j0(k*rparr[r])*(P1K(k)+P2K(k))/(2.*pi) 
       wpcalc.append(wpcal)
    
    #wparr = np.asarray(wpcalc)
    
    return wpcalc

#===========================================================================================


def wpfuncc0(z,nqobs,rparr,Mm,fsat,delm,kmin,kmax,dx,c0, B=0):
       
    bctrl = 2.7
    warnings.filterwarnings('ignore')
    
    e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)
    
    M_star,om0,delta_vir = 2.0e13,0.308,200.#,0.01,1000.

    b,fN,nqcalc = biasint2(z,nqobs,bctrl,delm,Mm)
    
    if B > 0:
        b = B
    
    wpcalc = []
  
    qmf = 'mVector_PLANCK-SMTz'+str(z)+'.txt' # quasar mass function
    mall = loadtxt(qmf)
    M = mall[:,0]
    dnm = mall[:,5]
    
    p2k = []
    p1k = []
    dkp = 0.5
    K   = np.logspace(log10(kmin),log10(kmax),300)
  
    for k in K:#np.logspace(log10(kmin),log10(kmax)-0.5,(kmax-kmin)/dkp):
            
       k0 = kmin
            
       p2k.append(b**2.*e.P(k,n=1.)) 
       
       p1calc = 0.
            
       for p in range(len(M)-1):
   
                  Nm = fN*meanNM(M[p],delm,Mm)
         
                  uk0 = ukmcalc0(k0,M[p],z,delta_vir,om0,M_star,c0)

                  u = ukmcalc0(k,M[p],z,delta_vir,om0,M_star,c0)/uk0
                  
                  dm =  M[p+1]/M[p]
                  
                  p1calc  += M[p]*ln(10)*dm*dnm[p]*(2.*fsat*(1.-fsat)*u+fsat**2.*abs(u)**2.)*Nm**2./nqobs**2.                   
                  
       p1k.append(p1calc)

    p1k = np.asarray(p1k)
    p2k = np.asarray(p2k)
        
  
    P1K = interpolate.interp1d(K, p1k)
    P2K = interpolate.interp1d(K, p2k) 
    
    
    for r in range(len(rparr)):
      
       #dx = 0.3   
       dkp = dx/rparr[r]
      
       wpcal = 0.
       for k in np.arange(K[1], np.max(K), dkp):
            
             wpcal += k*dkp*j0(k*rparr[r])*(P1K(k)+P2K(k))/(2.*pi) 
       wpcalc.append(wpcal)
    
    #wparr = np.asarray(wpcalc)
    
    return wpcalc
 
#==============================================================================


def main(plot=True):
  
  #dat = loadtxt('WP_obs.txt')
  #rparr = dat[:,0]
  ##wpd = dat[:,1]
  ##err = dat[:,2]
  
  #rparr=rparr[0:14]
  #print len(rparr)
  
  z = 1.55
  
  rparr = np.logspace(-2, log10(100),num=30, endpoint=True )
  
  wparr= wpfunc(z,8.813713331743382e-06,rparr,5.24903365032e+12,0.0326323936368,0.823258600156, 0.01, 1000, 0.5, B=0) # fitting parameters are Mm,fsat,delm
  print 'wparr=',wparr
  
  #wparr2= wpfunc(z,8.35e-6,rparr,1.33e13,0.01,0.75, 0.01, 1000, 0.5, B=0) # fitting parameters are Mm,fsat,delm 
  #print 'wparr2=',wparr2
  
  #wparr3= wpfunc(z,8.35e-6,rparr,1.33e13,0.054,0.75, 0.01, 1000, 0.5) # fitting parameters are Mm,fsat,delm 
  #print 'wparr3=',wparr3
  
  if True:
   
    plt.plot(rparr,wparr, 'b-')
    #plt.plot(rparr,wparr2, 'g--')
    #plt.plot(rparr,wparr3, 'r--')
    
    dat = loadtxt('rp_wp_err.dat')
    rpd = dat[:,0]
    wpd = dat[:,1]
    err = dat[:,2]
    
    plt.errorbar(rpd[0:4],wpd[0:4],yerr=err[0:4],linestyle='None',marker='o',markersize=3.,color='b',fillstyle='full',elinewidth=0.5,ecolor='orange',capsize=2.0,capthick=1.)  
    
    #plt.errorbar(rpd[4:17],wpd[4:17],yerr=err[4:17],linestyle='None',marker='o',markersize=3.,color='r',fillstyle='full',elinewidth=0.5,ecolor='blue',capsize=2.0,capthick=1.)  

    plt.errorbar(rpd[4:18],wpd[4:18],yerr=err[4:18],linestyle='None',marker='o',markersize=3.,color='g',fillstyle='full',elinewidth=0.5,ecolor='black',capsize=2.0,capthick=1.) 
  
    plt.yscale('log')
    plt.xscale('log') 
    
    plt.savefig('Bestfit_wpmodel_z'+str(z)+'.eps')

    plt.show()

    
if __name__ == '__main__':
  main()
  
 