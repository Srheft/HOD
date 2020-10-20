pro normalization_DD

;;;;;;  NOTE :
;         -  "_fib" and "_par" subscripts in the code blow are 
:           explained right below equation 9 here: http://export.arxiv.org/pdf/2007.09005
;         - "lseps" and "hseps" are the lower and higher edges of the logarithmic angular bins [in degrees] 
;        
;;;;;


datfilename='eBOSS_LRG_NGC_pip_v7_2_new.dat.fits'


lseps = [0. ,        0.01,       0.01217357, 0.01481958, 0.01804072];, 0.021962,0.0267356,  0.03254677, 0.03962104, 0.04823295, 0.05871672, 0.07147921,
0.08701572];, 0.1059292 , 0.12895366, 0.15698264, 0.19110393, 0.23264171,0.28320803, 0.34476529, 0.41970246, 0.51092774, 0.62198149, 0.75717355, 0.9217
5056];, 1.12209954, 1.36599578 ,1.66290459,2.02434862,2.46435506]

hseps = [0.01,    0.01217357, 0.01481958, 0.01804072, 0.021962];,0.0267356,  0.03254677, 0.03962104, 0.04823295, 0.05871672, 0.07147921,0.08701572, 0.1
059292]; , 0.12895366, 0.15698264, 0.19110393, 0.23264171, 0.28320803, 0.34476529, 0.41970246, 0.51092774, 0.62198149, 0.75717355, 0.92175056, 1.122099
54];, 1.36599578 ,1.66290459,2.02434862,2.46435506,3.0]

name = strsplit(datfilename,'pip',/EXTRACT)
tgt = name[0]

catalog = mrdfits(datfilename,1)

;compute normalization factors
;subsample cat
cat=catalog
c=0
for i=0L,(n_elements(cat)-1)/100 do begin
   cat[i]=cat[i*100]
   c+=1
endfor
cat=cat[0:c-1]

norm_par=long64(n_elements(cat))*long64(n_elements(cat)-1)/2.
a = where(cat.FIBER OR cat.CLUSTERING)
print,'unweighted counts all/clustering', n_elements(cat), n_elements(a)
cat = cat[a]

norm_fib=0.
q=intarr(n_elements(cat))
for i=0L,n_elements(cat)-2 do begin
   for j=i+1L,n_elements(cat)-1 do begin
      PIP_fib=0
      for k=0,59 do begin
         a=cat[i].WEIGHT_BW[k] AND cat[j].WEIGHT_BW[k]
         PIP_fib += bit_population(a)
      endfor
      norm_fib += 1860./PIP_fib 
   endfor
endfor


; "_fib" and "_par" subscripts are explained right below equation 9 here: http://export.arxiv.org/pdf/2007.09005

print,'normalization params (norm_par, norm_fib) are ', norm_par, norm_fib

