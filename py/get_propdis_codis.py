# -*- coding: utf-8 -*-
from cosmo import dist
from sys import argv
from math import pi
import numpy as np

def propdis(sep,z,om,ol):
    
    d = dist(om,ol)  # h=1
    fac = pi/180./3.6
    da = dist(omega=0.315,lamda=0.685,h=1).da
    #Rprop=d.ptdflat(z,sep)
    Rprop=round(sep*fac*da(z),1)
    #print Rprop
    return Rprop   
if __name__ == '__main__':

   #if len(argv) == 1:    # Then no arguments were passed
        #print "Syntax ERROR: No redshift and angular separation is specified!"
        
   #else:
        #z = float(argv[1])
        #sep = float(argv[2])
        #om = float(argv[3])
        #ol = float(argv[4])
   #dummy = propdis(z,sep, om, ol)

   for z in np.linspace(0.1,3.0,num=50):
       print z,propdis(2.9,z,0.315,.685),propdis(2.9,z,0.315,.685)*(1+z)