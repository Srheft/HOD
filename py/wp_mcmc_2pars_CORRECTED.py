import sys
import os
import random
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from math import *
from time import time
from astropy.table import Table, Column 
from astropy.io import ascii
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
from numpy import loadtxt, zeros, matmul
import numpy as np
from time import time, sleep
import HODemceeFIT_v4
from HODemceeFIT_v4 import biasint
from wpfunc4hod import wpfunc,wpfunc10c
import warnings
from numpy.linalg import inv

#
######################################################################
def median_error(X):
  
  size = len(X)
  X    = np.sort(X) 
  mean = np.median(X)
  
  u_err = X[int(round(0.84*size))] - mean
  l_err = mean - X[int(round(0.16*size))]
  
  return mean, u_err, l_err

######################################################################
def xi2(data, up_err, down_err, model):
  
  sigma   = 2*up_err*down_err/(up_err+down_err)
  sigma_p = abs(up_err-down_err)/(up_err+down_err) 
  
  Xi    = (model-data)/(sigma+sigma_p*(data-model))  
  
  return Xi**2,model
######################################################################
def xi2_cov(data,icov,model):
  
  Xi =0
  for i in range(len(data)):
       
      for j in range(len(data)):
          
         Xi2 += (data[i]-model[i])*icov[i,j]*(data[j]-model[j])
               
  return Xi2,model         
        

######################################################################
def Xi2(fsat, Mm):
  warnings.filterwarnings('ignore')
  #dat = loadtxt('rp_wp_err.dat')#wpobs_err_4mcmc.dat

  dat = loadtxt('rp_wp_err_ONLY1KDE.dat')
  rp = dat[:,0]
  wpobs = dat[:,1]
  uerr = dat[:,2]
  
  rp = np.asarray(rp)
  wpobs = np.asarray(wpobs)
  uerr= np.asarray(uerr)
  lerr = uerr
   
  z,nqobs,om0,delta_vir,kmin,kmax,M_star,dx,delm = 1.55,6.46178e-06,0.308,200.,0.01,1000,2.e13,0.5,0.75
  
  chisq,model = xi2(wpobs, uerr, lerr, wpfunc(z,nqobs,rp,Mm,fsat,delm,kmin,kmax,dx, B=0))
  
  return chisq,model
######################################################################
def Xi2_cov(fsat, Mm):
  warnings.filterwarnings('ignore')
  #dat = loadtxt('rp_wp_err.dat')#wpobs_err_4mcmc.dat

  dat = loadtxt('rp_wp_err_ONLY1KDE.dat')
  rp = dat[:,0]
  wpobs = dat[:,1]
  uerr = dat[:,2]
  cov = loadtxt('covmatrix_eboss.dat')
  icov = inv(cov)
  
  rp = np.asarray(rp)
  wpobs = np.asarray(wpobs)
  uerr= np.asarray(uerr)
  lerr = uerr
   
  z,nqobs,om0,delta_vir,kmin,kmax,M_star,dx,delm = 1.55,6.46178e-06,0.308,200.,0.01,1000,2.e13,0.5,0.75
  
  chisq,model = xi2_cov(wpobs, icov, wpfunc(z,nqobs,rp,Mm,fsat,delm,kmin,kmax,dx, B=0))
  
  return chisq,model
  ####################################################################
def main():

#### M A I N 
#### Initial Values - Arbitrary - Change them to see their effect
#### Convergence must happen with a good
#set of results .... 
 sd = sys.argv[1]
 
 nchain = sys.argv[2]

 z,nqobs,om0,delta_vir,kmin,kmax,M_star,dx,delm = 1.55,6.46178e-06,0.308,200.,0.01,1000,2.e13,0.5,0.75
 
 random.seed(int(sd))

 fsat = 0.01#started with Kayo&Oguro 2012 best fit parameters fsat = 0.054 
 Mm =1.0e+13 #from the 2000 batch# 5.86585430561e+12 from the 1000 batch #started with Mm=1.33e13 from Kayo & Oguri 2012 # a = 0.5
 
 sigma_Mm  = 0.3*Mm      # a guess value, arbitrarty
 sigma_fsat  = 0.3*fsat      # c guess value, arbitrary 

#### 
 #xhi2, model = Xi2(fsat, Mm)
 xhi2_cov, model = Xi2_cov(fsat, Mm)
 
 #chi2 = sum(xhi2)
 chi2 = sum(xhi2_cov)
 
 filename = 'results_seed'+str(sd)+'_nchain'+str(nchain)+'_2pars_fsat'+str(fsat)+'_Mm'+str(Mm)+'_delm'+str(delm)+'.txt'
  
 F = open(filename,'w')

 i = 0
 myTable = Table()
 empty = []
 myTable.add_column(Column(data=empty,name='n', dtype=np.dtype(int)))
 myTable.add_column(Column(data=empty,name='fsat', dtype=np.dtype(float)))
 myTable.add_column(Column(data=empty,name='Mm', format='%0.4f'))
 myTable.add_column(Column(data=empty,name='dm', dtype=np.dtype(float)))
 myTable.add_column(Column(data=empty,name='chi2', format='%1.2e'))
 
 for p in range(len(model)):
     myTable.add_column(Column(data=empty,name='wp'+str(p), format='%0.4f'))

 
 while i < int(nchain):

    fsat_new    = random.normalvariate(fsat, 0.3*sigma_fsat)
    Mm_new    = random.normalvariate(Mm, 0.3*sigma_Mm)
    
    if all([fsat_new > 0.005,fsat_new <0.09 ,Mm_new >1.e11 ,Mm_new <1.e15]):  
        
      xhi2_new, model_new = Xi2(fsat_new, Mm_new)
      chi2_new = sum(xhi2_new)
            
      delta_chi2 = chi2-chi2_new
      
      if delta_chi2>1: 
         ratio = 1
      
      else: 
      
         ratio = exp(delta_chi2)
    
      if random.uniform(0, 1.) < ratio:
           
         fsat = fsat_new
         Mm = Mm_new
         chi2 = chi2_new
         model = model_new
         
         #modr = np.logspace(log10(3.40e-2),log10(1.1),num=10, endpoint=True )
         
         #wpmodr= wpfunc(z,nqobs,modr,Mm,fsat,delm, 0.01, 1000, 0.5, B=0)

         #allwp = model[0:4]+wpmodr+model[4:19]
         allwp = model#[0:4]+model[4:model.shape[0]-1]
         row = [i, fsat, Mm, delm, chi2]
         for l in range(len(allwp)):
             row.append(allwp[l])
        
         myTable.add_row(row)
         
         if i%10==0:
             myTable.write(filename, format='ascii.fixed_width',delimiter=' ', bookend=False)
             print i, chi2
         i+=1   # loop index
 
 myTable.write(filename, format='ascii.fixed_width',delimiter=' ', bookend=False)

######################################################################
# Interpreting the results ....
## Calculating the media, uerr, and lerr values
#mytable = Table.read(filename, format='ascii', delimiter=',')
#Mm_lst = mytable['Mm']
#fsat_lst = mytable['fsat']

#print ''
#print 'Mm= ', median_error(Mm_lst)
#print 'fsat= ', median_error(fsat_lst)
#print ''

######################################################################

if __name__ == '__main__':
  main()
  

  
  
