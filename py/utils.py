
from time import clock, sleep, time
import IO
from glob import glob
import pickle
from sys import argv
from numpy import *

if verbose:
         print 'Loaded random catalogue....t =',time()-self.s,'sec'
         print rantop-1, 'points in random catalogue'

         if flag == True:
            print ' '
            print 'WARNING: Fewer random points than data points.  This'
            print 'might lead to statistical imprecision in Landy-Szalay'
            print 'estimator.  INCREASE SIZE OF INPUT RANDOM CATALOGUE!!'
            print ' '
            while True:
               x = raw_input("Continue Anyway? (y/n)>: ")
               if x == 'n' or x == 'N':
                  return
               elif x == 'y' or x == 'yes':
                  break
               sleep(1)
         
         
print "took", clock()-s, "s"

#----------------------------------------------------------------
def splineDELsq(self,rk,z):

        aexp=1./(1.+float(z))    # expansion factor
        
# calculate matter density, vacuum density at desired redshift

        self.om_m=self.omega_m(aexp,self.om_m0,self.om_v0) 
        self.om_v=self.omega_v(aexp,self.om_m0,self.om_v0)

# calculate the amplitude of the power spectrum at desired redshift 
# using linear growth factors (gg given by
# Caroll, Press & Turner 1992, ARAA, 30, 499)

        grow=self.gg(self.om_m,self.om_v) 
        grow0=self.gg(self.om_m0,self.om_v0)
        self.amp=aexp*grow/grow0

#        print 'Effective spectral quantities:'

        self.rknl=interpolate.splev(z,self.rknlspline)
        self.rneff=interpolate.splev(z,self.rneffspline)
        self.rncur=interpolate.splev(z,self.rncurspline)

#        print 'rknl [h/Mpc] = '+str(self.rknl)+' rneff= '+str(self.rneff)+' rncur= '+str(self.rncur)+'\n'

        self.gam=0.86485+0.2989*self.rneff+0.1631*self.rncur
        a=1.4861+1.83693*self.rneff+1.67618*self.rneff*self.rneff+0.7940*self.rneff*self.rneff*self.rneff+0.1670756*self.rneff*self.rneff*self.rneff*self.rneff-0.620695*self.rncur
        self.a=10**a      
        self.b=10**(0.9463+0.9466*self.rneff+0.3084*self.rneff*self.rneff-0.940*self.rncur)
        self.c=10**(-0.2807+0.6669*self.rneff+0.3214*self.rneff*self.rneff-0.0793*self.rncur)
        self.xmu=10**(-3.54419+0.19086*self.rneff)
        self.xnu=10**(0.95897+1.2857*self.rneff)
        self.alpha=1.38848+0.3701*self.rneff-0.1452*self.rneff*self.rneff
        self.beta=0.8291+0.9854*self.rneff+0.3400*self.rneff**2

        if(abs(1-self.om_m) > 0.01):  # omega evolution 
            f1a=self.om_m**(-0.0732)
            f2a=self.om_m**(-0.1423)
            f3a=self.om_m**(0.0725)
            f1b=self.om_m**(-0.0307)
            f2b=self.om_m**(-0.0585)
            f3b=self.om_m**(0.0743)       
            frac=self.om_v/(1.-self.om_m) 
            self.f1=frac*f1b + (1-frac)*f1a
            self.f2=frac*f2b + (1-frac)*f2a
            self.f3=frac*f3b + (1-frac)*f3a
        else:
            self.f1=1.0
            self.f2=1.
            self.f3=1.

# linear power spectrum !! Remember => plin=k^3*P(k)*constant
#
# constant = 4*pi*V/(2*pi)^3 

        plin = self.amp*self.amp*self.p_cdm(rk)

# calculate nonlinear power according to halofit: pnl = pq + ph,
# where pq represents the quasi-linear (halo-halo) power and 
# where ph represents the self-correlation halo term. 

        return self.halofit(rk,plin) # halo fit formula
        
#----------------------------------------------------------------
def splinewkernelkz(self,k,z):

        dzdchi = self.HubPar(z)/self.speedc

        a = self.splineDELsq(k,z)/k/k
        b = j0(self.codist(z)*self.theta*k)
        c = interpolate.splev(z,self.dndz)**2.

        return a*b*c*dzdchi
#----------------------------------------------------------------
def wkernelkz(self,k,z):

        dzdchi = self.HubPar(z)/self.speedc

        a = self.DELsq(k,z)/k/k
        b = j0(self.codist(z)*self.theta*k)
        c = interpolate.splev(z,self.dndz)**2.
        
        return a*b*c*dzdchi
#----------------------------------------------------------------
