pro NzMatcher_leftPANEL_FIG1,hardcopy=hardcopy

binsize = 0.2
zcp = [1.838,1.545,0.778,1.961,2.195,1.911,1.677,1.744,1.956,1.666,2.173,1.164,2.062,1.795,1.445,1.596,1.529,2.180,1.686,1.597,1.721,1.799,1.534,1.324,1.245,0.912,1.506,1.151,1.153,1.584,1.376,1.905,0.863,1.809,1.249,1.535,1.494,0.870,2.018,0.770,1.587,1.961,2.080,1.597]

zmin = min(zcp)
zmax = max(zcp)

Nbins= floor((zmax-zmin)/binsize)
edges= indgen(Nbins+1)*binsize+zmin
midz = edges[0:Nbins]+binsize/2.

datfits = 'eboss_v1.9f-QSO-N-eboss_v1.9f.dat.fits'

a = mrdfits(datfits,1)

wz = where(a.z le zmax and a.z ge zmin)
sp = a[wz]
wg= where(sp.MODELMAG[1] lt 20.85)
spg=sp[wg]

hzc = histogram(zcp,bin=binsize,locations=loc)
heb = histogram(spg.z,bin=binsize,locations=loc)

eb_z = []
ran_z = []
zcp_m=[]
eb_remain_ind=[]

for i=0L, Nbins-1 do begin

    print, 'from ', edges[i],' to ', edges[i+1]
    wb = where(spg.z lt edges[i+1] and spg.z gt edges[i],cb)
    wcp = where(zcp lt edges[i+1] and zcp gt edges[i],c)

    N2 = float(n_elements(spg))
    q2 = heb[i]
    p2 = q2/N2
   
    N1 = float(n_elements(zcp))
    q1 = hzc[i]
    p1 = q1/N1
    
    print, n_elements(wb),'  p1= ',p1,'  p2=',p2,'  p1*N2= ',p1*N2, '  p2*N2= ',p2*N2
    
    m = ceil(((p1*N2*1.)-(p2*N2*1.0))/(p1-1.))
    
    
    if (m gt 0.) then begin
         print,'bin #'+strtrim(i,1)+'m ='+strtrim(m,1)+' p2 ='+strtrim((q2-m)/(N2-m),1)+' p1= '+strtrim(p1,1)
   	 ind = Round( RandomU(666,m) * n_elements(wb))
    	 remove, ind, wb
    endif else begin
         m = floor(((p2*N1*1.)-(p1*N1*1.0))/(p2-1.))
         
         if (m gt 0.) then begin
            print,'bin #'+strtrim(i,1)+'m =',strtrim(m,1)+'matched p2 ='+strtrim((q2-m)/(N2-m),1)+'p1= '+strtrim(p1,1)
            ind = Round( RandomU(666,m) * n_elements(wcp))
            
            remove, ind, wcp
            PRINT, '******************* Removed KDE pairs******************'
            PRINT, ind, zcp[ind]
            PRINT,'********************************************************'
         endif
    endelse 
    
    eb_z = [eb_z, spg[wb].z]
    eb_remain_ind = [eb_remain_ind,wb]
    zcp_m = [zcp_m, zcp[wcp]]

endfor


mwrfits,spg[eb_remain_ind],'matched_eBOSS_'+strtrim(n_elements(eb_z),1)+'_zmin'+strtrim(zmin,1)+'_zmax'+strtrim(zmax,1)+'_gmag_LessThen20.85.fits',/create


wpmjd = where(spg[eb_remain_ind].MJD gt 0.)
obs_dates = spg[eb_remain_ind].MJD

PRINT,'############### the remaining targets are taken over dates: ',min(obs_dates[wpmjd]),'  to',max(obs_dates[wpmjd])

histZeb = histogram(eb_z,bin=binsize)
histZcp = histogram(zcp_m,bin=binsize)


print, '------------------------------------------------------'
print, 'The remaining KDE pairs ', total(histZcp),'*2'
print, 'The remaining eBOSS quasars ', total(histZeb)


print, '------------------------------------------------------'

nh_cp = (1/binsize)*histogram(zcp_m,bin=binsize)/total(histogram(zcp_m,bin=binsize))
alt_nh_cp  = [0.675,0.453,0.558,0.969,0.99,0.99,0.405]

nhistZcp= alt_nh_cp*binsize
if (keyword_set(hardcopy)) then PS_start,filename='FIG1_dNdz.eps',/encapsulated, bits_per_pixel=24, /color, /helv,xsize=10.5, ysize=7

  plot,[midz],(1/binsize)*histZeb/total(histZeb), psym=10, xtitle='z',ytitle='(1/N)(dN/dz)',xthick=7,ythick=7,charthick=2.5,xstyle=1,ystyle=1,charsize=2.0,yrange=[0,1.8],xrange=[0.7,2.5]

  oplot,[midz],(1/binsize)*histZeb/total(histZeb), psym=10, color=cgcolor('blue')
  oploterror,midz,(1/binsize)*histZeb/total(histZeb),(1/binsize)*sqrt(histZeb)/total(histZeb), errcolor=cgcolor('blue'), psym=3, color=cgcolor('blue')

  oplot,[midz],alt_nh_cp, psym=10, color=cgcolor('red'),linestyle=5
  oploterror,midz,alt_nh_cp,(1/binsize)*sqrt(histZcp)/total(histZcp), errcolor=cgcolor('red'), psym=3, color=cgcolor('red')
  
  ; ********************  creating edges for the histogram's first and last bins ****************
  ;vertical line for the starting eboss bin
  oplot,intarr(6)+edges[0],0.001+findgen(6,increment=(((1/binsize)*(histZeb/total(histZeb))[0])+0.11)*(1/6.)),color=cgcolor('blue')
  
  ;vertical line for the ending eboss bin
  oplot,intarr(6)+edges[n_elements(edges)-1],0.02+findgen(6,increment=(((1/binsize)*(histZeb/total(histZeb))[n_elements(histZeb)-1])+0.075)*(1/6.)),color=cgcolor('blue')
  
  ;vertical line for the starting KDE  bin 
  oplot,intarr(6)+edges[0],0.02+findgen(6,increment=(((1/binsize)*(histZcp/total(histZcp))[0])+0.11)*(1/6.)),color=cgcolor('red'),linestyle=5
  
  ;vertical line for the ending KDE bin
  oplot,intarr(6)+edges[n_elements(edges)-1],0.02+findgen(6,increment=(((1/binsize)*(histZcp/total(histZcp))[n_elements(histZcp)-1])+0.01)*(1/6.)),color=cgcolor('red'),linestyle=5

  ;horizontal line for starting bin 
  oplot,[edges[0],midz[0]],intarr(2)+(alt_nh_cp)[0],color=cgcolor('red'),linestyle=5
  oplot,[edges[0],midz[0]],intarr(2)+(1/binsize)*(histZeb/total(histZeb))[0],color=cgcolor('blue')
  
  ;Horizontal line for the ending bin
  oplot,reverse([edges[n_elements(edges)-1],midz[n_elements(midz)-2]]),intarr(2)+(alt_nh_cp)[n_elements(histZcp)-1],color=cgcolor('red'),linestyle=5
  oplot,reverse([edges[n_elements(edges)-1],midz[n_elements(midz)-2]]),intarr(2)+(1/binsize)*(histZeb/total(histZeb))[n_elements(histZeb)-1],color=cgcolor('blue')
  
  ; ********************  ***************************************************** ****************
 
   xyouts,2.1,1.60,'___  eBOSS',color= cgcolor('blue'),charsize=1.4
   xyouts,2.1,1.48,'_ _ _  KDE ',color= cgcolor('red'),charsize=1.4
   
  
  if (keyword_set(hardcopy)) then PS_end,/png

 

stop
end









