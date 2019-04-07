pro kdebox_OCT2018,hardcopy=hardcopy


dir='/home/ehsan/HaloOccupation/KDE-complete-box/'
readcol,dir+'z_R_rp_2.9arcsec.dat',z,Rpm,mcodis ; run get_propdis_codis.py to get this list 
readcol,dir+'z_R_rp_7.7arcsec.dat',z,Rpx,xcodis ; run get_propdis_codis.py to get this list 

;---------------------------------------------------------------------------------------------------------------------------
readcol,'allbinaries_lensREMOVED_25_49_23_21_6_2_1_145.cosmo',N,sep,id1,zem1,id2,zem2,x2,v,propsep,Format='F,F,A,F,A,F,F,F,F'

good = where(abs(v) lt 2000. ) 
;-------------------------------------------------------------------------------------------------------------------------------
readcol,'allbinaries_OCT2018.cosmo',N_supply,sep_supply,id1_supply,zem1_supply,id2_supply,zem2_supply,x2_supply,v_supply,propsep_supply,cosep_supply,codis,Format='F,F,A,F,A,F,F,F,F,F'
;-------------------------------------------------------------------------------------------------------------------------------
readcol,'allbinaries_AUG16_4plot.cosmo',Ng8,sepg8,id1g8,zem1g8,id2g8,zem2g8,x2g8,vg8,propsepg8,Format='F,F,A,F,A,F,F,F,F'
wg8=where(zem1g8 gt 0. and zem2g8 gt 0.)
wwg8=where(abs(vg8[wg8]) lt 2000.)
;-------------------------------------------------------------------------------------------------------------------------------
zmin=0.75
zmax=2.2

wm=where(z ge zmin)
print, wm[0]
wx=where(z ge zmax)
print, wx[0]
;-------------------------------------------------------------------------------------------------------------------------------


if (keyword_set(hardcopy)) then ps_start, /encapsulated, filename = 'Rp_z_OCT2018.eps', bits_per_pixel=24, /color, /helv,xsize=12, ysize=7
;v=['10','20','30','40','50','60','70','80']

plot,[0],[0],/xlog,xrange=[10.,200.],yrange=[0.2,3.0],xtitle='rp (kpc/h)',ytitle='z (redshift)',charsize=2.0,charthick=5.0,xstyle=1,ystyle=1,xthick=7,ythick=7;,XTICKV=['10','20','30','40','50','60','70','80'],xticks=7

;***************************************** PLOTTING LINES **************************************************
hline,zmin,linestyle=1
hline,zmax,linestyle=1

vline,mcodis[wm[0]],linestyle=2
xyouts,mcodis[wm[0]]+0.3,3.05,strnsignif(mcodis[wm[0]],3)+' kpc/h',charsize=1.1,charthick=5

vline,xcodis[wx[0]],linestyle=2
xyouts,xcodis[wx[0]]+0.5,3.05,strnsignif(xcodis[wx[0]],3)+' kpc/h',charsize=1.1,charthick=5

vline,43.35,linestyle=0,color=cgcolor('red')
vline,95,linestyle=0,color=cgcolor('red')

oplot,mcodis,z,linestyle=3,color=cgcolor('blue'),thick=7
xyouts,16.2*2.5,2.7,'theta=2.9',charsize=1.3,charthick=7,color=cgcolor('blue')

oplot,xcodis,z,linestyle=3,color=cgcolor('forest green'),thick=7
xyouts,43.5*2.5,2.7,'theta=7.7',charsize=1.3,charthick=7,color=cgcolor('forest green')

;****************************************PLOTTING OBJECTS **************************************************
plotsym,0,fill=0,thick=5
oplot,propsep*(1+zem1),zem1,psym=8,color=cgcolor('hot pink'),symsize=1.3,thick=5

oplot,propsepg8*(1+zem1g8),zem1g8,psym=8,color=cgcolor('hot pink'),symsize=1.3,thick=5

plotsym,0,fill=1,thick=5
oplot,propsep[good]*(1+zem1[good]),zem1[good],psym=8,color=cgcolor('orange'),symsize=1.1,thick=5

plotsym,0,fill=1,thick=5
oplot,propsepg8[wwg8]*(1+zem1g8[wwg8]),zem1g8[wwg8],psym=8,color=cgcolor('orange'),symsize=1.1,thick=5

if (keyword_set(hardcopy)) then PS_end,/png

print,'i+1, N[good],sep[good],id1[good],zem1[good],id2[good],zem2[good],x2[good],v[good],propsep[good],codis[good]'
for i=0, n_elements(zem1[good])-1 do begin
   print,i+1, N[good[i]],sep[good[i]],id1[good[i]],zem1[good[i]],id2[good[i]],zem2[good[i]],x2[good[i]],v[good[i]],propsep[good[i]],propsep[good[i]]*(1+zem1[good[i]])
endfor

stop
end
