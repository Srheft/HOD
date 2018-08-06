# -*- coding: utf-8 -*-
# This code is written with the purpose of fitting the halo model description used by **Kayo & Oguri 2012** that 
# divides the power spectrum into 1-halo and 2-halo form
# P(k)=P1h(k)+P2h(k)
# wp(rp)=integral(kdk/2pi)P(k)J0(krp)
# P2h(k)=b^2 Plin(k) where b=(1/nq)integral(bh(M) <N(M)> dn/dM dM)

from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
from scipy.special import spherical_jn,j0,sici # scipy.special.sici(x[, out1, out2]) 
from scipy.integrate import romberg,quad
import scipy.integrate as integrate
from cosmo import dist
#from cosmo_gen import distance
from random import gauss
from numpy import loadtxt, zeros,inf
import numpy as np
from time import time, sleep
from astropy.table import Table, Column
from numpy.fft import fft
from scipy import fftpack
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import warnings

 # TO RUN:
 # >>> import HODemceeFIT_v3
 # >>> HODemceeFIT_v2.wpint(1.4,7.6*0.001,0.75,1.33e13,0.054,0.11,1.1054046050405841e11,3.135,40,80)  ; ~ 460
 # >>> HODemceeFIT_v2.P1hint_one(1.4,7.6*0.001,0.75,1.33e13,0.054,0.11,1.1054046050405841e11,3.135,40,80)  

class ellipsoidalcollapse:

    def __init__(self,omega=0.307,lamda=0.693,sig8=0.8137,gams=0.21,acc=1e-9,H0=67.8,omb=0.045):
        """Sheth, Mo & Tormen (2001) ellipsoidal collapse model

        INPUTS: Dark Matter halo mass, omega_matter_0, omega_lamda_0,
        sigma_8, power spectrum shape, integration accuracies, Hubble
        constant (written with dimensions to avoid confusion with
        little h as we only use it for Eisenstein & Hu models) and the
        fraction of the total universe in baryons.

        v1.0 Adam D. Myers: August 2006
        """
#        self.T = self.TEBW # choose a transfer function

# Eisenstein and Hu only T(k) parameters below
        self.brat = omb/float(omega)
        self.hspec = H0/100.
        self.mat = float(omega)*self.hspec*self.hspec
        self.bar = float(omb)*self.hspec*self.hspec
        
# Eisenstein and Hu only T(k) parameters above

        self.T = self.TEisHu
        #if self.T == self.TEisHu:
            #print "Using Eisenstein/Hu baryon formalism:::::"
            #print "Shape of power spectrum meaningless"
            #print "Using omega_baryons instead"
            #print "(omega_b/omega_m) =",self.brat
            #print "(omega_m*h*h) =",self.mat-----------------
            #print "(omega_b*h*h) =",self.bar
        self.om0 = float(omega)
        self.ol = float(lamda)
        self.acc = float(acc)
        self.D = dist(omega,lamda,h=1.).D # linear growth factor G(z)
        self.omz = dist(omega,lamda,h=1.).omz # omega_matter(z)
        self.shape = float(gams)
        self.h0 = 100.
        self.P0 = 1.
        self.r = 8. # For normalization purposes this is in
# units of h-1Mpc
        correc = (sig8/self.siglog())**2.
        self.P0 = correc

        regracc = str(100.*(self.siglog()-sig8)/sig8)+'%'
        #print "*******************************************"
        #print "REGRESSED sigma_8 =",self.siglog()
        #print "INTEGRATION ACCURATE TO", regracc
        #print "*******************************************"
            

    def MDMH(self,MDMH=3e12): # mass should be in units of h-1 M_sol
        """Initialize mass, return radius, variance for given mass"""
        
        self.rho0 = 2.78e11*self.om0 # local mean density of universe in
# units of h**2M_sol.
        self.r = (3.*MDMH/(4.*pi*self.rho0))**(1/3.)

        self.sig = self.siglog()

#        print "INITALIZED MASS TO",MDMH
#        print "EFFECTIVE RADIUS FOR THIS MASS",self.r
#        print "EFFECTIVE rms FLUCTUATIONS FOR THIS MASS", self.sig
        
        return self.sig, self.r

    def nucalc(self,M,z):
        """Return the value of nu the "peak-height" for a given mass and redshift
        see, e.g, Equation 3 of Tinker et al. (2010)
        delta_c is fixed to 1.686 as in Tinker et al.!"""

        self.MDMH(MDMH=M)
        delc = 1.686
        nu = delc/self.sig/self.D(z)
        
        return nu

    def Tinkerbq(self,M,z):
        """bias from Tinker (2010) fitting function"""
        
        delta = 200
        y = log10(delta)
        capA = 1.0 + 0.24*y*exp(-(4/y)**4)
        a = 0.44*y -0.88
        capB = 0.183
        b = 1.5
        capC = 0.019 + 0.107*y + 0.19*exp(-(4/y)**4)
        c = 2.4

        nu = self.nucalc(M,z)
#        print 'log (i.e. ln) nu is',log10(nu)
        
        delc = 1.686

        b = 1 - capA*(nu**a/(nu**a + delc**a)) + capB*(nu**b) + capC*(nu**c)

        return b, log10(nu)


    def bq(self,z):
        """bias from Sheth, Mo & Tormen (2001) collapse model"""

        try:
            self.sig
        except AttributeError,NameError:
            print "*******************************************"
            print "* RUN MDMH first to initialize mass scale *"
            print "*******************************************"
            raise IOError

        aa = 0.707
        sqrta = sqrt(0.707)
        bb = 0.5
        cc = 0.6
        delc = 0.15*((12*pi)**(2/3.))*(self.omz(z)**0.0055)
        anusq = aa*((delc/self.sig/self.D(z))**2.)

        return 1.+((1/sqrta/delc)*((anusq*sqrta)+(sqrta*bb*(anusq**(1-cc)))-((anusq**cc)/(anusq**cc)+(bb*(1-cc)*(1-(cc/2.))))))

    def siglog(self):
        """rms mass fluctuations (integrated in logspace)"""

        return sqrt((ln(10.)/2./pi/pi)*romberg(self.intefunclog,-5,5,tol=self.acc))

    def intefunclog(self,logk):
        """Function to integrate to get rms mass fluctuations (logspace)"""

        k = 10.**logk
        kr = k*self.r

        return k*k*k*self.P(k)*(w(kr)**2)

    def P(self,k,n=1.):
        """Linear power spectrum"""

        return self.P0*(self.T(k)**2)*(k**n)

# All TRANSFER FUNCTION BELOW ASSUME T=2.7K.  Changing to T=2.73K
# increases amplitudes by about 2% (2.73/2.7)**2.

    def TBard(self,k):
        """Bardeen (1986) transfer function for adiabatic CDM"""
# Note, only one h in denominator to maintain hMpc-1 units of k
        q = k/self.shape
        
        return (ln(1+2.34*q)/2.34/q)*(1+3.89*q+(16.1*q)**2.+(5.46*q)**3.+(6.71*q)**4.)**(-0.25)

    def TEBW(self,k):
        """Efstathiou, Bond, White (1992) transfer function"""
# Note the returned quantity is squared in the Bardeen formalism
# hence the 1/nu rather than 2/nu in the denominator
        aa = 6.4/self.shape
        bb = 3.0/self.shape
        cc = 1.7/self.shape
        nu = 1.13
        denom = (1.+((aa*k)+((bb*k)**1.5)+((cc*k)**2.))**nu)**(1./nu)
      
        return 1./denom

    def TEisHuPure(self,k):
        """Eisenstein, Hu (1998) *PURE (no baryons)* transfer"""
# See equations 28 and 29 of Eisenstein & Hu (1998)
        q = k/self.shape
        Cno = 14.2+(731./(1+(62.5*q)))
        Lno = ln((2*e)+(1.8*q))
        return Lno/(Lno+(Cno*q*q))

    def TEisHu(self,k):
        """Eisenstein, Hu (1998) baryons included transfer"""
# See equations 28, 29, 30 and 31 of Eisenstein & Hu (1998)

        sh = (44.5*ln(9.83/self.mat))/(sqrt(1+(10.*((self.bar)**0.75))))
        sh *= self.hspec
#        alph = 1.-(0.328*ln(431.*self.mat)*self.brat)-(0.38*ln(22.3*self.mat)*(self.brat**2.))

        alph = 1.-(0.328*ln(431.*self.mat)*self.brat)+(0.38*ln(22.3*self.mat)*(self.brat**2.))
        gameff = self.om0*self.hspec
        gameff*=(alph+((1.-alph)/(1.+((0.43*k*sh)**4.))))
        q = k/gameff
        Cno = 14.2+(731./(1+(62.5*q)))
        Lno = ln((2*e)+(1.8*q))
        return Lno/(Lno+(Cno*q*q))
    
    def Hz(self,z):
        
        ok=1.-self.om0-self.ol
        return sqrt(self.h0**2.*(self.om0*(1+z)**3.+ok*(1+z)**2.0+self.ol))

def w(kr):
    """FT of a spherical top hat"""
    return 3.*(sin(kr)-(kr*cos(kr)))/kr/kr/kr

def altw(kr):
    """Equivalent to w(kr) from SPHERICAL bessel function"""
    return 3.*sph_jn(1,kr)[0][1]/kr


#=================================================== OCCUPATION FUNCTION <N(M)> ========================

def meanNM(M,delm,Mm):
     
    return (1./(sqrt(2.*pi)*delm))*exp(-(ln(M/Mm)*ln(M/Mm))/(2.*delm*delm))

#===========================================================================================================
def density(M,Mmin,Mmax,dndm,fN,delm,Mm):
    
    
    n = integrate.quad(dndm*meanNM(M,fN,delm,Mm),Mmin,Mmax)
    
    return n

#=======================================================================================================
# SE: ADM suggested to follow the equation 4,5 of Press & Schechter 1974

def nst_PS(z):
   mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')

   M = mall[:,0]
   dndm = mall[:,5]
   nst = 0.
   
   for i in range(len(M)-1):
          if dndm[i] >= 0 :
            dlog10M = log10(M[i+1]/M[i])  
            nst += dndm[i]*M[i]*dlog10M*ln(10)
            
   return nst
#=======================================================================================================
# SE: ADM suggested to follow the equation 4,5 of Press & Schechter 1974
def mst_PS(z):


   mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')

   M = mall[:,0]

   dndm = mall[:,5]
   
   nst = nst_PS(z)
   summ = 0.            
   for i in range(len(M)-1):  
       if dndm[i] >= 0 :
            dlog10M = log10(M[i+1]/M[i])
            summ += M[i]*dndm[i]*M[i]*dlog10M*ln(10)
            
   return summ/nst
#=================================================================================================
# SE: ADM suggested to follow the equation 4,5 of Press & Schechter 1974

def rhobar_PS(z):

    mst= mst_PS(z)
    nst = nst_PS(z)
    
    return mst*nst
#================================================== LINEAR BIAS ========================================
def biasint(z,nqobs,delm,Mm):
    warnings.filterwarnings('ignore')
    """ 
    The dn/dM is read in from a file
    created by hmfcalc online http://hmf.icrar.org. The passed redshift is
    the equivalent redshift for the calculations. Make sure the dn/dM file
    matches the chosen set of cosmological parameters and redshift and is
    named 'mVector_PLANCK-SMTz'+z+'.txt'
    """
   
    mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')

    M = mall[:,0]
    dndm = mall[:,5]

    dndlog10m = mall[:,7]

    # the cosmologiocal parameters (h,om_m) are the one's from Planck's 2013 ArXiv-1303.5076 and the sigma_8=0.8137 and om_b*h*h=0.022032
    # omb=0.022032/(0.67*0.67)=0.04907997326798841
    e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)

    numsum = 0
    denomsum = 0

    #ADM the log interval in mass for dlog10M in the integrator
    
    for i in range(len(M)-1):
          dlog10M = log10(M[i+1]/M[i])
          if dndm[i] >= 0 :
            b = e.Tinkerbq(M[i],z)
            num = dndm[i]*b[0]*M[i]*dlog10M*ln(10)*meanNM(M[i],delm,Mm)
            numsum += num
            denom = dndm[i]*M[i]*dlog10M*ln(10)*meanNM(M[i],delm,Mm)
            denomsum += denom
    
    blin = numsum/denomsum
    
    denomsum = numsum/blin
    nqcalc=denomsum
         
    fN = nqobs/nqcalc
    #alternative integration using dndlog10dM

    #for i in range(len(M)):
        #if all([dndm[i] >= 0 ]):
            #b = e.Tinkerbq(M[i],z)
            #num = dndlog10m[i]*b[0]*dlog10M
            #numsum += num
            #denom = dndlog10m[i]*dlog10M
            #denomsum += denom

    #blog = 10**(numsum/denomsum)
    #nqcalclog = denomsum
    #fNlog = nqobs/nqcalclog
    
    return blin ,fN ,nqcalc#, blog, fNlog, nqcalclogs 

#==================================================================================
def biasint2(z,nqobs,bctrl,delm,Mm):
    warnings.filterwarnings('ignore')
    """ 
    The dn/dM is read in from a file
    created by hmfcalc online http://hmf.icrar.org. The passed redshift is
    the equivalent redshift for the calculations. Make sure the dn/dM file
    matches the chosen set of cosmological parameters and redshift and is
    named 'mVector_PLANCK-SMTz'+z+'.txt'
    """
   
    mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')

    M = mall[:,0]
    dndm = mall[:,5]

    dndlog10m = mall[:,7]

    # the cosmologiocal parameters (h,om_m) are the one's from Planck's 2013 ArXiv-1303.5076 and the sigma_8=0.8137 and om_b*h*h=0.022032
    # omb=0.022032/(0.67*0.67)=0.04907997326798841
    e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)

    numsum = 0
    denomsum = 0

    #ADM the log interval in mass for dlog10M in the integrator
    
    for i in range(len(M)-1):
          dlog10M = log10(M[i+1]/M[i])
          if dndm[i] >= 0 :
            b = e.Tinkerbq(M[i],z)
            num = dndm[i]*b[0]*M[i]*dlog10M*ln(10)*meanNM(M[i],delm,Mm)
            numsum += num
            denom = dndm[i]*M[i]*dlog10M*ln(10)*meanNM(M[i],delm,Mm)
            denomsum += denom
    
    blin = numsum/denomsum
    
    if (bctrl >= 1.1*blin):
       nqcalc=denomsum
    else:
       
       blin = bctrl
       denomsum = numsum/blin
       nqcalc=denomsum
         
    fN = nqobs/nqcalc
    #alternative integration using dndlog10dM
    #numsum = 0.
    #denomsum = 0.

    #for i in range(len(M)):
        #if all([dndm[i] >= 0 ]):
            #b = e.Tinkerbq(M[i],z)
            #num = dndlog10m[i]*b[0]*dlog10M
            #numsum += num
            #denom = dndlog10m[i]*dlog10M
            #denomsum += denom

    #blog = 10**(numsum/denomsum)
    #nqcalclog = denomsum
    #fNlog = nqobs/nqcalclog
    
    return blin ,fN ,nqcalc#, blog, fNlog, nqcalclogs 

#========================================TWO HALO TERM OF THE POWER SPECTRUM ==========================================

def P2hk(z,b,rp):
   # The two halo term according to eq.(8) of Kayo & Oguri 2012 is: P2h(k)=b**2 Plin(k) where "Plin" is the ``linear power spectrum''
   # from ADM's twopt.py code  Plin(k)= const. K**3 P(k), const= 4piV/(2pi)**3
   #V=dist(omega,lamda,h=1.).Vcfunc(z)
   #cosnt=4.*pi*/(2pi)**3
   
   #plin = loadtxt('CAMPoutput.txt')
   #k = plin[:,0]
   #e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.315*0.67,acc=1e-9,H0=67.8,omb=0.045)
   #plin=e.P(k,n=1.)
   #b=biasint(z,nqobs,delm,Mm)[0]
   #return plin*b**2.
   
   mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')

   M = mall[:,0]
   dndm = mall[:,5]
   dndlog10m = mall[:,7]

   kall = loadtxt('kVector_PLANCK-SMTz'+str(z)+'.txt')
   k = kall[:,0]

   plin = kall[:,1]
   numsum=0.
   for i in range(len(k)-1):
        #plin=e.P(k[i],n=1.)
        dlog10k = log10(k[i+1]/k[i])  
        num=plin[i]*((k[i])**2.0)*dlog10k*ln(10)*j0(k[i]*rp)/(2*pi)
        numsum += (num)*1.0
        
   return (b**2.)*numsum

#============================ mean concentration parameter of the density profile of the NFW dark matter halo ===============
#SE: the formula for the concentration parameter as a function of halo mass and redshift is taken from equation 78 of Cooray & Sheth 2002 
# NOTE that this is the *mean* concentration paramter that depends on mass unlike the width of its distribution. use the *c(r_s)* function  
# Richardson et al 2012 calls M_star the "non-lenear mass for collapse at z=0"
def cbar(M,z,M_star):
   
   return (9./(1.+z))*(M/M_star)**(-0.13)

#========================================================================
#SE: Richardson et al 2012 used the same formula(Eq. 3) but with c0=25 as opposed to 9 for z~1.4 modeling and changed that to c=45 fro z~3.2  
def cbarc0(c0,M,z,M_star):
   
   return (c0/(1.+z))*(M/M_star)**(-0.13)


#==================================== Virial radius from the dark matter halo mass =============================================
def rvir(om0,delta_vir,M):
    
    rho_bg = 2.78e11*om0
    
    return (3.*M/(4*pi*delta_vir*rho_bg))**(1/3.)

#============================  concentration parameter of the density profile of the NFW dark matter halo ===============
#def c(r_s,M):

   #'''
   #SE: According to the argument on page 19 of Cooray and Sheth 2002, virialized object is \Delta _vir= 18 pi**2. times the mean background density 
   #In the Spherical Collapse Model, an initially overdense region *turnaround* at \theta = pi and completely copllapse at \theta = 2*pi so equation 50 at theta = pi 
   #will become: R0/R(z_ta) = (3pi/4)/**(2/3.)
   #'''
   #c=r_vir/r_s , r_vir= R_ta/2.
   #cbar(M,z,M_star)
   #return r_vir/r_s
#=========================== Lagrangian size ========================================================
#def LagSize(z,theta):
    
    
    
#============================   Ci integral that is used in u(k,m) :Ci(x) =\int_{x}^{\inf} cos(t)/t dt ==============
# SE: the formula is taken from equation 82 of Cooray and Sheth (2002)
def ci_int(t):
    
    return (-1)*cos(t)/t
#============================   Si integral that is used in u(k,m) :Si(x) =\int_{0}^{x} sin(t)/t dt =================
# SE: the formula is taken from equation 82 of Cooray and Sheth (2002)
def si_int(t):
  
    return sin(t)/t
#============================   u(k,m): is the fourier transform of the halo desity profile   =======================
# SE: the formula is taken from equation 81 of Cooray and Sheth (2002)

#def u_km(k,M,rho_s,r_s):
      
      #coef= 4*pi*rho_s*r_s**3./M
      #term1 = ((integrate.quad(si_int,0,(c(M)+1)*k*r_s))[0]-(integrate.quad(si_int,0,k*r_s))[0])*sin(k*r_s)
      #term2 = sin(c(M)*k*r_s)/((1+c(M))*k*r_s)
      #term3 = ((integrate.quad(ci_int,(c(M)+1)*k*r_s,np.inf))[0]-(integrate.quad(ci_int,k*r_s,np.inf))[0])*cos(k*r_s)
      
      #return coef*(term1-term2+term3)
#============================   u(k,m): is the fourier transform of the halo desity profile   =======================
# SE: the formula is taken from equation 78 and 81 of Cooray and Sheth (2002)

def u_km(k,M,M_star,z,om0,delta_vir,kmax):
      
      r_s = rvir(om0,delta_vir,M) / cbar(M,z,M_star)
      myinf = np.inf#kmax*r_s
      
      rho_s = abs( M/(4.*pi*(r_s**3.0)*(log10(cbar(M,z,M_star)+1.)-(cbar(M,z,M_star)/(1.+cbar(M,z,M_star))))))
      coef = 4*pi*rho_s*r_s**3./M
      term1 = ((integrate.quad(si_int,0.,(cbar(M,z,M_star)+1.)*k*r_s))[0]-(integrate.quad(si_int,0.,k*r_s))[0])*sin(k*r_s)
      term2 = sin(cbar(M,z,M_star)*k*r_s)/((1+cbar(M,z,M_star))*k*r_s)
      term3 = ((integrate.quad(ci_int,(cbar(M,z,M_star)+1.)*k*r_s,myinf))[0]-(integrate.quad(ci_int,k*r_s,myinf))[0])*cos(k*r_s)
      
      return coef*(term1-term2+term3)
#====================================================================================================================
#SE: equation 53 from Cooray & Sheth 2002 
def delta_sc(z):
     
     return (3/5.)*(1+z)*(3*pi/2.)**(2/3.) 

#============this function doesnt give the correct mean density at z=0 while the formula and H(z) is correct !!!!!!!!!!===================
def rho_bar(z):
   #G=7.0136e10# in M_sun per (Mpc)^3
   #e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-6,H0=67.8,omb=0.05)
   #return 3.*e.Hz(z)**2./(8*pi*G)
   
   # SE: according to Cooray and Sheth 2002  mean density is :   \bar rho = \int m*n(m)*dm  and gives the same result as the rhobar_PS(z,lcut,hcut) that is introduced in this code [based on the ]
   mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')

   M = mall[:,0]

   dndm = mall[:,5]
   
   rhobar = 0.            
   for i in range(len(M)-1):  
       if dndm[i] >= 0 :
            dlog10M = log10(M[i+1]/M[i])
            rhobar += M[i]*dndm[i]*M[i]*dlog10M*ln(10)

   return rhobar
   
#====================================================================================================================
def Mstar(z,k,omega=0.308):
    
    #2.78e11*0.308
    return 4*pi*Rstar(z,k)**3.*rho_bar(z)/3.
#====================================================================================================================
def A(k):
    
     e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)
     plin=e.P(k,n=1.)
     
     return  plin/e.P0*k**(3./2.)
#====================================================================================================================
def Rstar(z,k):
    
    return ((A(k)*8.*sqrt(pi))/(delta_sc(z)**2.)*15.*pi**2.)**(2/3.)
#====================================================================================================================
    
def P1k(z,delm,Mm,fsat,nqobs):
    
    dk = 0.5
    kmax = 900.
    b,fN,nqcalc = biasint(z,nqobs,delm,Mm)
    nq2 = nqcalc**2.0 #nqobs**2. #nq2=(1.35e-6)**2.#biasint(z,nqobs,delm,Mm)[1]
    
    mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')
    M = mall[:,0]
    
    dlog10M = log10(M[30]/M[29])

    ukmfile = loadtxt('ukm_z'+str(z)+'_dk'+str(dk)+'_kmax'+str(kmax)+'.dat')
    mhalos = ukmfile[:,0]
    dndm =  ukmfile[:,1]
    sorted_k_accept =ukmfile[:,2]
    sorted_ukm_accept = ukmfile[:,3]

    
    numsum =0.
    for i in range(len(M)):
        if dndm[i] >= 0 :
               num = dndm[i]*mhalos[i]*dlog10M*ln(10)*((fN*meanNM(mhalos[i],delm,Mm))**2.)*( 2.*fsat*(1.-fsat)*(sorted_ukm_accept[i])+(fsat**2.)*(abs(sorted_ukm_accept[i])**2.) )/(nq2)
               numsum += (num)*1.0
               
    p1k = numsum      

    return p1k

#====================================================================================================================
def P2k(z,b,k,fN,delm,Mm):
    
    
    #kall = loadtxt('kVector_PLANCK-SMTz'+str(z)+'.txt')
    #kvec = kall[:,0]
    #kvec = kvec[0:1386] # k from 1.e-4 --> 100 h/Mpc
    #plin = kall[:,1]
    #plin = plin[0:1386]
    
    #p2 = interp1d(kvec,plin,kind='cubic')
    #p_intp= p2(k)
    #return (b**2.)*p_intp
    
    e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)
    plin = e.P(k,n=1.)
    
    return (b**2.)*plin

    
#====================================================================================================================
def onehaloterm(z,delm,Mm,fsat,rp,nqobs):
        
    dk = 0.5
    kmax = 900.
    b,fN,nqcalc = biasint(z,nqobs,delm,Mm)
    nq2 =(1.35e-6)**2.# nqcalc**2.0 #nqobs**2. #nq2=(1.35e-6)**2.#biasint(z,nqobs,delm,Mm)[1]
    
    mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')
    M = mall[:,0]
    dlog10M = log10(M[30]/M[29])

    ukmfile = loadtxt('ukm_z'+str(z)+'_dk'+str(dk)+'_kmax'+str(kmax)+'.dat')
    mhalos = ukmfile[:,0]
    dndm =  ukmfile[:,1]
    sorted_k_accept =ukmfile[:,2]
    sorted_ukm_accept = ukmfile[:,3]

    
    numsum =0.
    for i in range(len(mhalos)):
        
               num = dndm[i]*mhalos[i]*dlog10M*ln(10)*((fN*meanNM(mhalos[i],delm,Mm))**2.)*( (2.*fsat*(1.-fsat)*sorted_ukm_accept[i])+(fsat**2.*(sorted_ukm_accept[i])**2. ))*(sorted_k_accept[i])*dk*j0(sorted_k_accept[i]*rp)/(2*pi*nq2)
               numsum += (num)*1.0
               
    wprp= numsum        
            
    return wprp
#---------------------------------------------------------------------------------------------------------------------------------------


def P1hint_one(z,fN,delm,Mm,fsat,rp,M_star,kmin,kmax,dk,nq2):
  
#numbers are from Scoccimarro, Sheth, Hui and Jain 2001
#delm=340 for omega=0.3, delm=200[NFW] for omega=1 
#rho_star=rho_s=M_s/(4/3*pi*r_s**3.) = 1.1054046050405841e11   (M_sun/h)/(Mpc/h)^3
#M_star=M_s=1.07e13 for LCDM M_sun/h
#r_star=r_s=3.135 Mpc/h
    
    mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')
    M = mall[:,0]

    dlog10M = log10(M[2]/M[1])
    
    r_s = rvir(om0,delta_vir,M) / cbar(M,z,M_star) #6.0
    rho_s =1.e11
    ukmfile = loadtxt('ukm'+str(z)+'rs'+str(r_s)+'rhos'+str(rho_s)+'_dk'+str(dk)+'_kmax'+str(kmax)+'.dat')
    
    mhalos = ukmfile[:,0]
    dndm =  ukmfile[:,1]
    sorted_k_accept =ukmfile[:,2]
    sorted_ukm_accept = ukmfile[:,3]

    #nq2=(1.35e-6)**2.#biasint(z,nqobs,delm,Mm)[1]
    numsum =0.
    for i in range(len(mhalos)):
            
               
               num = dndm[i]*mhalos[i]*dlog10M*ln(10)*(meanNM(mhalos[i],fN,delm,Mm)**2.)*(sorted_ukm_accept[i])*(dk)*dk*j0(sorted_k_accept[i]*rp)/(2*pi)
               numsum += (num)*1.0
               sall = numsum
            
            
    return 2.*fsat*(1.-fsat)*sall/nq2
#====================================================================================================================

def P1hint_two(z,fN,delm,Mm,fsat,rp,M_star,kmin,kmax,dk,nq2):
  
    
    mall = loadtxt('mVector_PLANCK-SMTz'+str(z)+'.txt')
    M = mall[:,0]
    
    dlog10M = log10(M[2]/M[1])
    
    r_s = rvir(om0,delta_vir,M) / cbar(M,z,M_star)#6.0
    rho_s= 1.e11
    ukmfile = loadtxt('ukm'+str(z)+'rs'+str(r_s)+'rhos'+str(rho_s)+'_dk'+str(dk)+'_kmax'+str(kmax)+'.dat')
    
    mhalos = ukmfile[:,0]
    dndm =  ukmfile[:,1]
    sorted_k_accept =ukmfile[:,2]
    sorted_ukm_accept = ukmfile[:,3]


    #nq2=(1.35e-6)**2.#biasint(z,nqobs,delm,Mm)[1]  # SE: got the observational nq=1.35e-6 @ z=1.4 from Kayo & Oguri 2012 
    numsum = 0.
    for i in range(len(mhalos)):
       
          num = dndm[i]*mhalos[i]*dlog10M*ln(10)*(meanNM(mhalos[i],fN,delm,Mm)**2.)*(abs(sorted_ukm_accept[i])**2.)*((sorted_k_accept[i]))*dk*j0(sorted_k_accept[i]*rp)/(2*pi)
          numsum += (num)*1.0
          intg=numsum
      
        
    return (fsat**2.)*intg/nq2
#====================================================================================================================
#wpint(1.4,0.0076,0.75,1.33e13,0.054,0.11,1.1054046050405841e11,3.135,) #All values are from Kayo & Oguri 2012: delm=0.75, fN=0.0076, Mm=1.33e13, fsat=0.054(SIS profile), 0.048(NFW profile), c=10, rp in Mpc from Table2 
 #M_star = 1.07e13

def wpint(z,b,fN,delm,Mm,fsat,rp,M_star,kmin,kmax,dk):
        start = time()
     
        b = biasint(z,nqobs,delm,Mm)[0]
        
        part1=P1hint_one(z,fN,delm,Mm,fsat,rp,M_star,kmin,kmax,dk)+P1hint_two(z,fN,delm,Mm,fsat,rp,M_star,kmin,kmax,dk)
     
        part2=P2hk(z,b,fN,delm,Mm,rp)
     
        wp=part1+part2
     
        #print "Done!....took ",(time()-start)/60., "min"
        return wp#,P1hint_one(z,fN,delm,Mm,fsat,rp,rho_s,r_s)[2],P1hint_two(z,fN,delm,Mm,fsat,rp,rho_s,r_s), part2

#=============================================================================================================================
def fastwp(z,b,delm,Mm,fsat,rp,nqobs):
    
    p1= onehaloterm(z,delm,Mm,fsat,rp,nqobs) 
    
    p2=P2hk(z,b,rp)
        
    return p1+p2
      
#=============================================================================================================================
 
def wp_new(z,b,fN,delm,Mm,fsat,rp,rho_s,r_s,M_star,kmin,kmax,dk):
    
        k = np.arange(kmin,kmax,dk)

        w = 0.
        for i in range(len(k)-1):
            
            Pk = P1k(z,k[i],fN,delm,Mm,fsat,rho_s,r_s,M_star)+ P2k(z,b,k[i],fN,delm,Mm)
        
            w += Pk*k[i]*dk*j0(k[i]*rp)/(2.*pi)
    
        return w
#=============================================================================================================================

def main():
   
   start = time()
   
   z,delm,Mm,fsat,M_star,kmin,kmax,dk,om0,delta_vir,nqobs = 1.55,0.75,1.33e13,0.054,1.03e13,0.01,900,0.5,0.308,200.,8.813713331743382e-06
   
   warnings.filterwarnings('ignore')
   ##rpdis=[0.012,0.071,0.11,0.171,0.265,0.634,0.981,1.518,2.349,3.636,5.627]
   ##rpdis=np.logspace(-2, 2.3, num=30,endpoint=True)
   ##rpdis=np.logspace(-2, 2.3, num=30,endpoint=True)
   ##k = np.arange(0.01,10.0,0.01)
   
   ##e = ellipsoidalcollapse(omega=0.308,lamda=0.692,sig8=0.8137,gams=0.308*0.678,acc=1e-9,H0=67.8,omb=0.05)
   
   b,fN,nqcalc = biasint(z,nqobs,delm,Mm)
   
   rpdis=[0.012,0.071,0.11,0.171,0.265,0.634,3.636,5.627,10.,20.,30.]
   print ''
   print 'Working with ','ukm_z'+str(z)+'_dk'+str(dk)+'_kmax'+str(kmax)+'.dat'
   print ''
   
   for i in range(len(rpdis)):
     binstart = time()
     print rpdis[i],'Mpc/h ',fastwp(z,b,delm,Mm,fsat,rpdis[i],nqobs),(time()-binstart), "sec"
     #wpint(z,b,fN,delm,Mm,fsat,rpdis[i],M_star,kmin,kmax,dk,lcut,hcut)
     #wp_new(z,b,fN,delm,Mm,fsat,rpdis[i],rho_s,r_s,M_star,kmin,kmax,dk,lcut,hcut)#
if __name__ == '__main__':
  start = time()  
  main()
  
  
  print "Done!....took ",(time()-start)/60., "min"
