; quick code to prepare datacubes for BUDDI.
; The code will:
;   - convert flux back into counts in a really dirty way (/counts)
;   - log-rebin the datacube in the wave_log direction if
;     it's linearly binned (/log_rebin)
;

;===============================================================


;===============================================================

pro BUDDI_manga_prep,directory,file,MUSE=muse,JY=jy,BADPIX=badpix
  ;buddi_prep,'BUDDI_input.txt',/JY,/BADPIX

fits_read,directory+file+'.fits',input,h

x=sxpar(h,'NAXIS1')
y=sxpar(h,'NAXIS2')
z=sxpar(h,'NAXIS3')


newfile = directory+file+'.fits'
fits_read,newfile,wave_lin,h_wave,extname='WAVE'
wave_log=alog10(wave_lin)

sxaddpar,h,'CRVAL3',wave_log[0]
sxaddpar,h,'CDELT3',wave_log[1]-wave_log[0]
sxaddpar,h,'CD3_3',wave_log[1]-wave_log[0]



input_counts=fltarr(x,y,z)
if keyword_set(Jy) then begin
    print,'*** Now converting units to Jy'

    for zz=0,z-1,1 do begin
      ;convert from 10^-17 erg/s/cm2/AA to Jy
      input_counts[*,*,zz]=(input[*,*,zz]*1e-17)*wave_lin[zz]*wave_lin[zz]*3.34e4
    endfor
    sxaddpar,h,'BUNIT','Jy', 'Specific intensity (per spaxel)'

endif else input_counts=input







;output=fltarr(x,y,z)
output=input_counts


;print,'New wave_log solution for BUDDI input file:'
;print,'D00) wave_log start= '+string(wave_log[0])
;print,'D01) wave_log step=  '+string(wave_log[1]-wave_log[0])


fits_write,directory+file+'_FLUX.fits',output,h
temp=mrdfits(directory+file+'.fits',0,h_temp)
modfits,directory+file+'_FLUX.fits',0,h_temp

if keyword_set(Jy) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
sxaddpar,h_temp,'CD1_1',sxpar(h,'CD1_1')
sxaddpar,h_temp,'CD2_2',sxpar(h,'CD2_2')



;extract IVAR datacube and convert to sigma images
fits_read,directory+file+'.fits',IVAR,hdr2,extname='IVAR'
s=size(IVAR)
sigma=fltarr(s[1],s[2],s[3])
sigma[*,*,*]=1/sqrt(IVAR[*,*,*])

if keyword_set(Jy) then begin
  for zz=0,z-1,1 do sigma[*,*,zz]=(sigma[*,*,zz]*1e-17)*wave_lin[zz]*wave_lin[zz]*3.34e4
endif
    
;fits_write,directory+file+'_SIGMA.fits',sigma,extname='SIGMA'
;sxaddpar,h2,'EXTNAME','SIGMA'

if keyword_set(Jy) then sxaddpar,hdr2,'BUNIT','Jy', 'Specific intensity (per spaxel)'
if keyword_set(Jy) then sxaddpar,h_temp,'BUNIT','Jy', 'Specific intensity (per spaxel)'
;modfits,directory+file+'_SIGMA.fits',0,h_temp
;modfits,directory+file+'_SIGMA.fits',1,h2,extname='SIGMA'
  
fits_write,directory+file+'_SIGMA.fits',sigma,hdr2,extname='SIGMA'
;temp=mrdfits(directory+file+'.fits',0,h_temp)
modfits,directory+file+'_SIGMA.fits',0,h_temp
;modfits,directory+file+'_SIGMA.fits',1,hdr2,extname='SIGMA'

  ;==================================
  ;bad pixel datacube

if keyword_set(BADPIX) AND NOT file_test(directory+file+'_BADPIX.fits') then begin
    print,'*** Now creating bad pixel datacube'
    print,'creating bad pixel mask'
    fits_read,newfile,badpix,hdr3,extname='MASK'
    sxaddpar,hdr3,'EXTNAME','BADPIX'
    fits_write,directory+file+'_BADPIX.fits',badpix,hdr3,extname='BADPIX'
;    temp=mrdfits(directory+file+'.fits',0,h_temp)
    modfits,directory+file+'_BADPIX.fits',0,h_temp
;    modfits,directory+file+'_BADPIX.fits',1,hdr3,extname='BADPIX'
endif


print,'#############################'
print,'## CRVAL3='+string(wave_log[0])+' ##'
print,'## CDELT3='+string(wave_log[1]-wave_log[0])+' ##'
print,'#############################'

  
end
