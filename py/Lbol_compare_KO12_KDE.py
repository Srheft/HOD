from numpy import *
import numpy as np
from time import time
from math import acos, pi, log10, sqrt
from cosmo import dist
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

start = time()

# SE: Cosmology params used in Eftekharzadeh++2017
dis = dist(0.307,0.693,0.677)

#SE: i-magnitudes and redshifts of 44 close pairs(88 values) that participate in the Wp measurement-> taken from Table 5 of Eftekharzadeh et al 2017  
kde_imags=[20.16,21.07,20.82,20.94,17.24,20.25,17.99,19.75,20.43,20.66,19.73,20.10,18.93,19.83,18.58,20.29,20.04,20.61,20.74,20.84,20.24,20.63,17.58,20.18,20.66,20.78,20.11,20.37,19.69,19.77,19.87,20.07,19.64,19.72,19.99,20.09,20.33,20.54,19.90,19.94,18.55,20.05,20.26,20.90,20.06,20.48,20.03,20.82,19.39,20.27,19.64,20.06,19.38,19.45,20.80,20.94,19.97,20.65,20.25,20.83,19.82,19.99,20.49,20.80,19.88,20.44,18.58,20.83,20.17,20.25,20.64,20.70,19.70,20.66,19.81,20.44,20.08,20.46,18.29,18.50,20.33,20.42,19.11,20.60,19.45,19.77,18.78,20.32]

kde_z=[1.838,1.545,0.778,1.961,2.195,1.911,1.677,1.744,1.956,1.666,2.173,1.164,2.062,1.795,1.445,1.596,1.529,2.180,1.686,1.597,1.721,1.799,1.534,1.324,1.245,0.912,1.506,1.151,1.153,1.584,1.376,1.905,0.863,1.809,1.249,1.535,1.494,0.870,2.018,0.77,1.587,1.961,2.080,1.597]

#SE: i-magnitude and redshifts of 26 pairs (52 values) taken from Table 1 of Kayo & Oguri 2012
ko12_imags=[18.41,19.66,18.61,19.04,18.96,20.17,18.52,19.60,18.99,19.70,18.65,19.14,19.10,20.28,18.81,20.01,19.03,20.11,18.47,19.55,18.34,19.55,19.06,18.63,18.91,19.27,18.92,19.93,18.34,19.27,19.01,19.68,18.82,19.19,18.94,19.63,18.86,19.88,18.67,19.73,18.20,18.62,18.86,19.43,18.31,18.38,19.03,20.07,17.63,17.77,18.87,19.02]

ko12_z=[0.978,0.980,0.626,0.627,1.712,1.712,1.218,1.223,1.833,1.883,1.212,1.215,1.745,1.740,1.678,1.678,1.216,1.218,1.495,1.495,1.200,1.195,1.246,1.241,2.051,2.041,1.555,1.543,1.877,1.867,1.246,1.261,1.506,1.506,0.799,0.799,1.249,1.256,1.644,1.648,1.567,1.567,0.750,0.750,0.769,0.769,1.775,1.775,1.890,1.879,1.897,1.897]

#SE:


w=open('KDEpairs_z_mi_Mi_Lbole-47.dat','w+') # make them without .dat 

kde_absimag=[]
kde_Lbol=[]
kde_2z=[]

for j in range(len(kde_imags)):
    
          #print(j,int(j/2))
          
          i = int(j/2)
          
          dmz = dis.dm(kde_z[i])
          
          Absimag = kde_imags[j]-dis.Kcorr(kde_z[i])-dmz
          
          kde_absimag.append(Absimag)

          kde_2z.append(kde_z[i])
          #SE: combining eqn 1 of Shen++2009 and eqn 1 of Richards++2006 => Mi(z=0)=Mi(z=2)+0.596, where Mi(z=2)=90-2.5 log10(Lbol/erg/s)
          Lbol =  10**((Absimag-90.596)/2.5)*1e47 # in ergs/s 

          kde_Lbol.append(Lbol)

          string=str(kde_z[i])+'    '+str(kde_imags[j])+'    '+str(Absimag)+'    '+str(Lbol)
          
          #print(string)
          
          w.write(string+'\n')
   
w=open('KO12pairs_z_mi_Mi_Lbole-47.dat','w+') # make them without .dat 
ko12_absimag=[]
ko12_Lbol=[]

for j in range(len(ko12_imags)):
    
          #print(j,int(j/2))
          
          i = j
          
          dmz = dis.dm(ko12_z[i])
          
          Absimag = ko12_imags[j]-dis.Kcorr(ko12_z[i])-dmz
          
          ko12_absimag.append(Absimag)
          
          #SE: combining eqn 1 of Shen++2009 and eqn 1 of Richards++2006 => Mi(z=0)=Mi(z=2)+0.596, where Mi(z=2)=90-2.5 log10(Lbol/erg/s)
          Lbol =  10**((Absimag-90.596)/2.5)*1e47  # in ergs/s 
          
          ko12_Lbol.append(Lbol)
          
          string=str(ko12_z[i])+'    '+str(ko12_imags[j])+'    '+str(Absimag)+'    '+str(Lbol)
          #print(ko12_z[i],ko12_imags[j],Absimag,Lbol)
          
          w.write(string+'\n')
          
 
#SE: binning in delta_z=0.2 from 0.6 to 2.2

zm,zmax,delz=0.6,2.2,0.2
kde_2z = np.asarray(kde_2z)
ko12_z = np.asarray(ko12_z)
ko12_Lbol = np.asarray(ko12_Lbol)
kde_Lbol = np.asarray(kde_Lbol)

Lbin_kde = []
Lbin_ko12 = []
zmids = []
for i in range(8):
    
    zx = zm+delz
    print(i,zm,zx)
    zmids.append((zm+zx)/2.)
    w = np.where(kde_2z <=  zx) and np.where(kde_2z > zm)
    print(w)
    Lbin_kde.append(mean(kde_Lbol[w]))
    w = np.where(ko12_z <=  zx) and np.where(ko12_z > zm)
    Lbin_ko12.append(mean(ko12_Lbol[w]))
    zm = zx
          
counts_kde, bin_edges_kde = np.histogram(kde_Lbol, bins=5)
counts_ko12, bin_edges_ko12 = np.histogram(ko12_Lbol, bins=5)
print('HISTcounts_kde',counts_kde,'HISTcounts_ko12',counts_ko12) 
print('Average bolometric luminosity of KDE quasars:',np.mean(kde_Lbol),'e-47 ergs/s','Average bolometric luminosity of KO12 quasars:',np.mean(ko12_Lbol),'e-47 ergs/s','mean(kde_Lbol)/mean(ko12_Lbol):',np.mean(kde_Lbol)/np.mean(ko12_Lbol)) 





fig=plt.figure()
ax= fig.add_subplot(111)
# axes limits
x1=0.6; x2=2.5
y1=0.4; y2=50.0
   
# Genral axes
ax.set_xlim(x1, x2)
ax.set_ylim(y1, y2)
ax.minorticks_on()
  
# additional Y-axis (on the right)
y_ax = ax.twinx()
y_ax.set_ylim(y1, y2)
y_ax.set_yticklabels([])
y_ax.minorticks_on()

# additional X-axis (on the top)
x_ax = ax.twiny()
x_ax.set_xlim(x1, x2)
x_ax.set_xticklabels([])
x_ax.minorticks_on()
   
ax.set_ylabel(r'$\rm L_{bol}\, [\times 10^{-47} \, s^{-1}erg]$', fontsize=12)
ax.set_xlabel(r'$\rm z$', fontsize=12)
      
# Grid
ax.grid(False)

p1,= ax.plot(kde_2z,kde_Lbol,marker='o',markersize=4.,color='orange',linestyle='None', label=r'$\rm KDE \, QSOs$')   
p2,= ax.plot(ko12_z,ko12_Lbol,marker='s',markersize=4.,color='magenta',linestyle='None', label=r'$\rm KO12 \, QSOs$')   

lns = [p1,p2]

ax.legend(handles=lns, loc='best')
fig.savefig('Lum_compare.eps',bbox_inches='tight')




fig=plt.figure()
ax= fig.add_subplot(111)
# axes limits
x1=0.6; x2=2.3
y1=0.4; y2=12.0
   
# Genral axes
ax.set_xlim(x1, x2)
ax.set_ylim(y1, y2)
ax.minorticks_on()
  
# additional Y-axis (on the right)
y_ax = ax.twinx()
y_ax.set_ylim(y1, y2)
y_ax.set_yticklabels([])
y_ax.minorticks_on()

# additional X-axis (on the top)
x_ax = ax.twiny()
x_ax.set_xlim(x1, x2)
x_ax.set_xticklabels([])
x_ax.minorticks_on()


ax.set_ylabel(r'$\rm L_{bol}\, [\times 10^{-47} \, s^{-1}erg]$', fontsize=12)
ax.set_xlabel(r'$\rm z$', fontsize=12)
      
# Grid
ax.grid(False)
#SE: custom marker
h=delz
v=0.15


p1,= ax.plot(zmids,Lbin_kde,marker='o',markersize=2.,color='orange',linestyle='None', label=r'$\rm KDE \, QSOs$')  
a=(list(zmids),list(Lbin_kde))
a_zipped = zip(*a)
for a_x, a_y in a_zipped:

    ax.add_patch(Rectangle(xy=(a_x-h/2, a_y-v/2) ,width=h, height=v, linewidth=1, color='orange', fill=True))
    
    

p2,= ax.plot(zmids,Lbin_ko12,marker='s',markersize=2.,color='magenta',linestyle='None', label=r'$\rm KO12 \, QSOs$') 
a=(list(zmids),list(Lbin_ko12))
a_zipped = zip(*a)
for a_x, a_y in a_zipped:

    ax.add_patch(Rectangle(xy=(a_x-h/2, a_y-v/2) ,width=h, height=v, linewidth=1, color='magenta', fill=True))


lns = [p1,p2]

ax.legend(handles=lns, loc='best')
fig.savefig('binnedLum_compare.eps',bbox_inches='tight')

#print('kde_Lbol',kde_Lbol)
#print('ko12_Lbol',ko12_Lbol)
#kwargs = dict(histtype='stepfilled', alpha=0.3, normed=True, bins=5)
#plt.hist(kde_Lbol, **kwargs,label=r'$KDE-L_{bol}$')
#plt.hist(ko12_Lbol, **kwargs,label=r'$KO12-L_{bol}$')
plt.show()

  
print('Done....t =',(time()-start),'sec')

