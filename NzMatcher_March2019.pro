pro NzMatcher_March2019,hardcopy=hardcopy

binsize = 0.2
zcp = [1.838,1.545,0.778,1.961,2.195,1.911,1.677,1.744,1.956,1.666,2.173,1.164,2.062,1.795,1.445,1.596,1.529,2.180,1.686,1.597,1.721,1.799,1.534,1.324,1.245,0.912,1.506,1.151,1.153,1.584,1.376,1.905,0.863,1.809,1.249,1.535,1.494,0.870,2.018,0.770,1.587,1.961,2.080,1.597]

zmin = min(zcp)
zmax = max(zcp)

Nbins= floor((zmax-zmin)/binsize)
edges= indgen(Nbins+1)*binsize+zmin
midz = edges[0:Nbins]+binsize/2.

datfits = 'eboss_v1.9f-QSO-N-eboss_v1.9f.dat.fits'
ranfits = 'eboss_v1.9f-QSO-N-eboss_v1.9f.ran.fits'  ; 5135183

a = mrdfits(datfits,1)
r = mrdfits(ranfits,1)
wz = where(a.z le zmax and a.z ge zmin)
sp = a[wz]
wzr = where(r.z le zmax and r.z ge zmin)
ran = r[wzr]  ; 3978978


hzc = histogram(zcp,bin=binsize,locations=loc)
heb = histogram(sp.z,bin=binsize,locations=loc)
hr = histogram(ran.z,bin=binsize,locations=loc)

eb_z = []
ran_z = []
zcp_m = []
eb_remain_ind = []
ran_remain_ind = []

for i=0L, Nbins-1 do begin

    print, 'from ', edges[i],' to ', edges[i+1]
    
    wb = where(sp.z lt edges[i+1] and sp.z gt edges[i],cb)
    
    wcp = where(zcp lt edges[i+1] and zcp gt edges[i],c)

    N2 = float(n_elements(sp))
    
    q2 = heb[i]
   
    N1 = float(n_elements(zcp))
    
    q1 = hzc[i]
    
    n = floor(N2/N1)
    
    m = ceil(q2/n - q1)
    
    if (m gt 0.) then begin
    
   	 ind = Round( RandomU(666,m) * n_elements(wb))
    	 remove, ind, wb
         PRINT, '******************* Removed QSOs ******************'
         PRINT, n_elements(ind), zcp[ind]
         PRINT,'********************************************************'

    	 
    endif 
    
    eb_z = [eb_z, sp[wb].z]
    eb_remain_ind = [eb_remain_ind,wb]
    zcp_m = [zcp_m, zcp[wcp]]

endfor


mwrfits,sp[eb_remain_ind],'matched_eBOSS_'+strtrim(n_elements(eb_z),1)+'_zmin'+strtrim(zmin,1)+'_zmax'+strtrim(zmax,1)+'.fits',/create


spm=sp[eb_remain_ind]

spmg=spm[where(spm.MODELMAG[1] ge  17.68 and spm.MODELMAG[1] le 20.85)]
mwrfits,spmg,'matched_eBOSS_'+strtrim(n_elements(eb_z),1)+'_zmin'+strtrim(zmin,1)+'_zmax'+strtrim(zmax,1)+'_g20.85.fits',/create

wpmjd = where(sp[eb_remain_ind].MJD gt 0.)
obs_dates = sp[eb_remain_ind].MJD

 PRINT,'############### the remaining targets are taken over dates: ',min(obs_dates[wpmjd]),'  to',max(obs_dates[wpmjd])

;print,histogram(eb_z,bin=binsize),histogram(eb_z,bin=binsize)/total(eb_z)
;print,histogram(zcp_m,bin=binsize),histogram(zcp_m,bin=binsize)/total(zcp_m)

histZeb = histogram(eb_z,bin=binsize)
histZcp = histogram(zcp_m,bin=binsize)


for i=0L, Nbins-1 do begin

   print, 'from ', edges[i],' to ', edges[i+1]
    wb = where(sp.z lt edges[i+1] and sp.z gt edges[i],cb)
    wr = where(ran.z lt edges[i+1] and ran.z gt edges[i],c)

    N2 = float(n_elements(ran))
    q2 = hr[i]
    
   
    N1 = float(n_elements(sp[eb_remain_ind]))
    q1 = histZeb[i]
    
    nn = floor(N2/N1)
    
    
    mr = ceil(q2/nn - q1)
    
    
    if (mr gt 0.) then begin
    	 indr = Round( RandomU(666,mr) * n_elements(wr))
    	 remove, indr, wr
    endif 
    
    ran_remain_ind = [ran_remain_ind,wr]
    ran_z = [ran_z, ran[wr].z]

endfor


mwrfits,ran[eb_remain_ind],'matched_random_'+strtrim(n_elements(ran_z),1)+'_zmin'+strtrim(zmin,1)+'_zmax'+strtrim(zmax,1)+'.fits',/create

histZran = histogram(ran_z,bin=binsize)

print, '------------------------------------------------------'
print, 'The remaining KDE pairs ', total(histZcp),'*2'
print, 'The remaining eBOSS quasars ', total(histZeb)
print, 'The remaining eBOSS randoms ', total(histZran)
print, '------------------------------------------------------'


if (keyword_set(hardcopy)) then PS_start,filename='NEWmatched_Nz.eps',/encapsulated, bits_per_pixel=24, /color, /helv,xsize=10.5, ysize=7

  plot,[midz],histZeb/total(histZeb), psym=10, xtitle='z',ytitle='(1/N)(dN/dz)',xthick=7,ythick=7,charthick=2.5,xstyle=1,ystyle=1,charsize=2.0,yrange=[0,0.35],xrange=[0.7,2.5]

  oplot,[midz],histZeb/total(histZeb), psym=10, color=cgcolor('blue')
  oploterror,midz,histZeb/total(histZeb),sqrt(histZeb)/total(histZeb), errcolor=cgcolor('blue'), psym=3, color=cgcolor('blue')

  oplot,[midz],histogram(zcp_m,bin=binsize)/total(histogram(zcp_m,bin=binsize)), psym=10, color=cgcolor('red'),linestyle=5
  oploterror,midz,histZcp/total(histZcp),sqrt(histZcp)/total(histZcp), errcolor=cgcolor('red'), psym=3, color=cgcolor('red')
  
  ; ********************  creating edges for the histogram's first and last bins ****************
  
  oplot,intarr(6)+edges[0],0.001+findgen(6,increment=(((histZeb/total(histZeb))[0])+0.02)*(1/6.)),color=cgcolor('blue')
  oplot,intarr(6)+edges[n_elements(edges)-1],0.001+findgen(6,increment=(((histZeb/total(histZeb))[n_elements(histZeb)-1])+0.017)*(1/6.)),color=cgcolor('blue')
 
  oplot,intarr(6)+edges[0],0.001+findgen(6,increment=(((histZcp/total(histZcp))[0])+0.01)*(1/6.)),color=cgcolor('red'),linestyle=5
  oplot,intarr(6)+edges[n_elements(edges)-1],0.001+findgen(6,increment=(((histZcp/total(histZcp))[n_elements(histZcp)-1])-0.001)*(1/6.)),color=cgcolor('red'),linestyle=5

  oplot,[edges[0],midz[0]],intarr(2)+(histZcp/total(histZcp))[0],color=cgcolor('red'),linestyle=5
  oplot,[edges[0],midz[0]],intarr(2)+(histZeb/total(histZeb))[0],color=cgcolor('blue')
  
  oplot,reverse([edges[n_elements(edges)-1],midz[n_elements(midz)-2]]),intarr(2)+(histZcp/total(histZcp))[n_elements(histZcp)-1],color=cgcolor('red'),linestyle=5
  oplot,reverse([edges[n_elements(edges)-1],midz[n_elements(midz)-2]]),intarr(2)+(histZeb/total(histZeb))[n_elements(histZeb)-1],color=cgcolor('blue')
  
  ; ********************  ***************************************************** ****************
 
   ;legend,['eBOSS','KDE'],color=[cgcolor('blue'),cgcolor('red')],lines=[0,5],charsize=1.2,box=0,thick=5,pos=[2.3,0.30]
   xyouts,2.1,0.30,'___  eBOSS',color= cgcolor('blue'),charsize=1.2
   xyouts,2.1,0.28,'_ _ _  KDE ',color= cgcolor('red'),charsize=1.2
   
   ;legend,['eBOSS','KDE'],color=[cgcolor('blue'),cgcolor('red')],pos=[2.3,0.30]
  
  if (keyword_set(hardcopy)) then PS_end,/png


if (keyword_set(hardcopy)) then PS_start,filename='NEWmatched_RAN_eBOOS_NZ_v2.eps',/encapsulated, bits_per_pixel=24, /color, /helv,xsize=10.5, ysize=7

  plot,[midz],histZeb/total(histZeb), psym=10, xtitle='z',ytitle='(1/N)(dN/dz)',xthick=7,ythick=7,charthick=2.5,xstyle=1,ystyle=1,charsize=2.0,yrange=[0,0.35],xrange=[0.7,2.5]

  oplot,[midz],histZeb/total(histZeb), psym=10, color=cgcolor('blue')
  oplot,[midz],histZeb/total(histZeb), psym=10, color=cgcolor('forest green'),linestyle=3
   
  oploterror,midz,histZeb/total(histZeb),sqrt(histZeb)/total(histZeb), errcolor=cgcolor('blue'), psym=3, color=cgcolor('blue')
  oploterror,midz,histZeb/total(histZeb),sqrt(histZeb)/total(histZeb), errcolor=cgcolor('Forest green'), psym=3, color=cgcolor('blue')

  ;oplot,[midz],histZran/total(histZran), psym=10, color=cgcolor('Forest green'),linestyle=3
  ;oploterror,midz,histZran/total(histZran),sqrt(histZran)/total(histZran), errcolor=cgcolor('Forest green'), psym=3, color=cgcolor('Forest green')
  
  ; ********************  creating edges for the histogram's first and last bins ****************
    
  oplot,intarr(6)+edges[0],0.001+findgen(6,increment=(((histZeb/total(histZeb))[0])+0.02)*(1/6.)),color=cgcolor('blue')
  oplot,intarr(6)+edges[n_elements(edges)-1],0.001+findgen(6,increment=(((histZeb/total(histZeb))[n_elements(histZeb)-1])+0.017)*(1/6.)),color=cgcolor('blue')
  oplot,intarr(6)+edges[n_elements(edges)-1],0.001+findgen(6,increment=(((histZeb/total(histZeb))[n_elements(histZeb)-1])+0.017)*(1/6.)),color=cgcolor('Forest green'), linestyle=3
 
  oplot,intarr(6)+edges[0],0.001+findgen(6,increment=(((histZran/total(histZran))[0])+0.01)*(1/6.)),color=cgcolor('Forest green'),linestyle=3
  
  oplot,[edges[0],midz[0]],intarr(2)+(histZeb/total(histZeb))[0],color=cgcolor('blue')
  oplot,[edges[0],midz[0]],intarr(2)+(histZeb/total(histZeb))[0],color=cgcolor('forest green'),linestyle=3
  
  oplot,reverse([edges[n_elements(edges)-1],midz[n_elements(midz)-2]]),intarr(2)+(histZran/total(histZran))[n_elements(histZran)-1],color=cgcolor('Forest green'),linestyle=3
  
  oplot,reverse([edges[n_elements(edges)-1],midz[n_elements(midz)-2]]),intarr(2)+(histZeb/total(histZeb))[n_elements(histZeb)-1],color=cgcolor('blue')
 

  ; ********************  ***************************************************** ****************
   xyouts,2.1,0.30,'___  eBOSS',color= cgcolor('blue'),charsize=1.2
   xyouts,2.1,0.28,'_ . _ . _  Random ',color= cgcolor('Forest green'),charsize=1.2
   
  
  if (keyword_set(hardcopy)) then PS_end,/png
  
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  K-S test - CDF plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
;  KSTWO Returns the Kolmogorov-Smirnov statistic and associated probability that two arrays of data values are drawn from the same distribution :KSTWO, data1, data2, d, prob
;  D is the maximum difference between the two distribution functions
;  prob is the significance of statistic calculated by internal function PROB_KS


print, ''
KSTWO, eb_z, zcp_m, D_matched, prob_matched & print, 'The probablity of KDE-eBOSS are drawn from the same distribution is:', prob_matched, 'Max dist. b/w the 2 distributions :', D_matched
print, ''

print, ''
KSTWO, sp.z, zcp, D_orig, prob_orig & print, 'The probablity of original KDE-eBOSS to be drawn from the same distribution is:', prob_orig, 'Max dist. b/w the 2 distributions :', D_orig
print, ''

;parentKDE = mrdfits('NGC_parentsample_aftercuts_369556_WITHreplacements_ALLzassigned_da.fits',1)
;KSTWO,parentKDE.DR6Z_PSREPLACED,sp.z, D_parent, prob_parent & print, 'The probablity of parent KDE-eBOSS to be drawn from the same distribution is:', prob_parent, 'Max dist. ;b/w the 2 distributions :', D_parent
;print, ''


pdf_Zeb = histZeb/total(histZeb)
cdf_Zeb = TOTAL(pdf_Zeb, /CUMULATIVE) 

pdf_Zcp = histZcp/total(histZcp)
cdf_Zcp = TOTAL(pdf_Zcp, /CUMULATIVE) 


if (keyword_set(hardcopy)) then PS_start,filename='KDE_eBOSS_CDFs.eps',/encapsulated, bits_per_pixel=24, /color, /helv,xsize=10.5, ysize=7 

PLOT,[midz], cdf_Zcp, XTITLE='z',YTITLE='Cumulative redshift distribution', /nodata, charsize=1.8,yrange=[0,1.05],YSTYLE=1, XTHICk=7,YTHICk=7

;*********************************************   MARCH2019 *******************************************
cdf_Zcp = [0.10000000, 0.17500001, 0.27500001, 0.44999999, 0.69999999, 0.92499995, 0.99999994]
cdf_Zeb = [0.096277773, 0.18696139, 0.26145970, 0.46514622, 0.70130591, 0.89429583, 0.99999994]



OPLOT,[midz], cdf_Zcp, COLOR=cgcolor('red'), linestyle=2

OPLOT,[midz], cdf_Zeb, COLOR=cgcolor('blue')
   
;legend,['eBOSS','KDE'],color=[cgcolor('blue'),cgcolor('red')],pos=[0.9,0.012]
xyouts,0.8,0.9,'___  eBOSS',color= cgcolor('blue'),charsize=1.5
xyouts,0.8,0.85,'_ _ _  KDE ',color= cgcolor('red'),charsize=1.5

if (keyword_set(hardcopy)) then PS_end,/png
 
;NOTE: this was used for past HOD paper
  
stop
end









