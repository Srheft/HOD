import numpy as np
from math import *
from time import time,sleep
from astropy.table import Table, Column 
from astropy.io import ascii
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
from numpy import loadtxt, zeros, matmul
from numpy.linalg import inv,norm
from sklearn.preprocessing import normalize

dat = loadtxt('matrix_test.txt')


rp = dat[:,0]
wp = dat[:,1]
err = dat[:,2]
wpall = dat[:,3]

jkarr = dat[0:,4:]  #SE this represents all the coulmns including and after the 5th  

filename = 'covariance_matrix.txt'
myTable = Table()
empty =[]
F = open(filename,'w')

c= np.zeros((len(rp),len(rp)))

for i in range(len(rp)):
  #myTable.add_column(Column(data=empty,name='c'+str(i)+str(i), format='%0.4f'))
  for j in range(len(rp)): 
    
      for L in range(np.shape(jkarr)[1]):
          
          c[i,j] += (jkarr[i,L]-wp[i])*(jkarr[j,L]-wp[j])
  #myTable.add_row(list(c[i,:]))
  print(list(c[i,:]))


print(normalize(c))
#print(np.shape(c),c[1,2],c[2,1])
#myTable.write(filename, format='ascii.fixed_width',delimiter=' ', bookend=False)



