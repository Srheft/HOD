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
from numpy import loadtxt, zeros
from time import time, sleep
import HODemceeFIT_v4
from HODemceeFIT_v4 import biasint
from wpfunc4hod import wpfunc,wpfunc10c
import warnings
from numpy.linalg import inv,norm
from sklearn.preprocessing import normalize

# SE: wpfile has to be givn in form of a string
def covmaker(wpfile):
    dat=loadtxt(wpfile)
    rp = dat[:,0]
    wp = dat[:,1]
    wperr = dat[:,2]
    cov=np.zeros((len(rp),len(rp)))

    if np.shape(dat)[1] < 4:
        
         print('No measured covariance, making it form wp measurement')
         print('Warning: CHECK THE INPUT FILE FORMAT! - The assumption for the format of the input file is: col1: rp[Mpc/h], col2:wp[Mpc/h], col3:wperr')
         
         for i in range(len(rp)):
             
            cov[i,i]=wperr[i]**2.
    else:
         jkarr = dat[0:,4:]
         for i in range(len(rp)):
             for j in range(len(rp)):     
                 for L in range(np.shape(jkarr)[1]):
          
                     cov[i,j] += (jkarr[i,L]-wp[i])*(jkarr[j,L]-wp[j])

    return cov
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
  
  Xi2 = 0.
  for i in range(len(data)):
       
      for j in range(len(data)):
        #SE: see eq 5 of Eftekharzadeh et al 2015 
         Xi2 += (data[i]-model[i])*icov[i,j]*(data[j]-model[j]) 
  return Xi2,model         
######################################################################

def Xi2(wpfile,fsat, Mm):
    
  warnings.filterwarnings('ignore')

  dat = loadtxt(wpfile)
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
def Xi2_cov(wpfile,fsat, Mm):
    
  warnings.filterwarnings('ignore')
  
  dat = loadtxt(wpfile)
  rp = dat[:,0]
  wpobs = dat[:,1]
  uerr = dat[:,2]
  #SE: if the NORMALIZED covariance matrix is stored in a separate ascii file, then we can go there   
  cov = covmaker(wpfile)#loadtxt('covar_4KDE+eBOSS.dat')  
  icov = inv(cov)
  rp = np.asarray(rp)
  wpobs = np.asarray(wpobs)
  uerr= np.asarray(uerr)
  lerr = uerr
   
  z,nqobs,om0,delta_vir,kmin,kmax,M_star,dx,delm = 1.55,6.46178e-06,0.308,200.,0.01,1000,2.e13,0.5,0.75
 
  chisq,model = xi2_cov(wpobs, icov, wpfunc(z,nqobs,rp,Mm,fsat,delm,kmin,kmax,dx, B=0))
  print(chisq,model)
  return chisq,model
######################################################################

def main():

#### M A I N 
#### Initial Values - Arbitrary - Change them to see their effect
#### Convergence must happen with a good
#set of results .... 
wpfile = sys.argv[1]
dat = loadtxt(wpfile)
if np.shape(dat)[1] < 4:      
         print('No measured covariance, making it form wp measurement')
         print('Warning: CHECK THE INPUT FILE FORMAT! - The assumption for the format of the input file is: col1: rp[Mpc/h], col2:wp[Mpc/h], col3:wperr')

 sd = sys.argv[2]
 
 nchain = sys.argv[3]

 z,nqobs,om0,delta_vir,kmin,kmax,M_star,dx,delm = 1.55,6.46178e-06,0.308,200.,0.01,1000,2.e13,0.5,0.75
 
 random.seed(int(sd))

 fsat = 0.01#started with Kayo&Oguro 2012 best fit parameters fsat = 0.054 
 Mm =1.0e+13 #from the 2000 batch# 5.86585430561e+12 from the 1000 batch #started with Mm=1.33e13 from Kayo & Oguri 2012 # a = 0.5
 
 sigma_Mm  = 0.3*Mm      # a guess value, arbitrarty
 sigma_fsat  = 0.3*fsat      # c guess value, arbitrary 

#### 
 print('Starting ...')
 start = time()
 xhi2_cov, model = Xi2_cov(wpfile,fsat, Mm)
 print('Done calculating the model =',(time()-start),'sec')

 chi2 = xhi2_cov 
 
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
    
    if all([fsat_new > 0.008,fsat_new <0.080 ,Mm_new >1.e12 ,Mm_new <1.e15]):  
      
      xhi2_new, model_new = Xi2_cov(wpfile,fsat_new, Mm_new)
      
      chi2_new = xhi2_new
            
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
         
         allwp = model
         row = [i, fsat, Mm, delm, chi2]
         for l in range(len(allwp)):
             row.append(allwp[l])
        
         myTable.add_row(row)
         
         if i%10==0:
             myTable.write(filename, format='ascii.fixed_width',delimiter=' ', bookend=False)
             print(i, chi2)
         i+=1   # loop index
 
 ### Do not remove the following line !!!
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
  

  
  
