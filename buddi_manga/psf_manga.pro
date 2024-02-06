; quick code to create PSF datacubes for MaNGA datacubes
; in preparation of fitting the datacube with BUDDI

pro PSF_manga,dir,cube,dir_static

  ;Best way to do this for SDSS MaNGA
  ;data is to read in the FWHM for each waveband, and then interpolate
  ;between their central wavelengths.
  ;
  ;read in PSF extensions for each filter

  fits_read,dir+cube+'.fits',input_g,header_g,extname='GPSF'
  fits_read,dir+cube+'.fits',input_r,header_r,extname='RPSF'
  fits_read,dir+cube+'.fits',input_i,header_i,extname='IPSF'
  fits_read,dir+cube+'.fits',input_z,header_z,extname='ZPSF'

  h = headfits(dir+cube+'.fits')
  h2 = headfits(dir+cube+'.fits',exten=1)
  h_wave = headfits(dir+cube+'.fits',exten=6)
  

  sxaddpar,h,'EXPTIME',1
  sxaddpar,h,'GAIN',1
  x_side=sxpar(header_g,'NAXIS1')
  y_side=sxpar(header_g,'NAXIS2')
  z_side=sxpar(h2,'NAXIS3')
  ugriz_wave=[3543,4770,6231,7635,9134]



  ;for each image slice, use the wavelength and the filter transmission curves
  ;to determine the fraction of light from each filter and coadd them in the
  ;correct proportions.
  ;
  readcol,dir_static+'g.dat',format='F,X,F,X,X',wave_g,transmission_g,/SILENT
  readcol,dir_static+'r.dat',format='F,X,F,X,X',wave_r,transmission_r,/SILENT
  readcol,dir_static+'i.dat',format='F,X,F,X,X',wave_i,transmission_i,/SILENT
  readcol,dir_static+'z.dat',format='F,X,F,X,X',wave_z,transmission_z,/SILENT

;  wavelength0_lin=sxpar(h2,'CRVAL3')
;  wavelength0=alog10(wavelength0_lin)
;  step_lin=sxpar(h2,'CD3_3')
;  step=alog10(wavelength0_lin+step_lin)-wavelength0
;  wave_log=fltarr(sxpar(h2,'NAXIS3'))
;  for m=0,sxpar(h2,'NAXIS3')-1,1 do wave_log[m]=wavelength0+(m*step)
  fits_read,dir+cube+'.fits',wave_lin,h_wave,extname='WAVE'
  wave_log=alog10(wave_lin)

  psf_out=fltarr(x_side,y_side,z_side)
  psf_temp=fltarr(x_side,y_side)

  for j=0,z_side-1,1 do begin
    wave=wave_lin[j]
    total_light=0
    if wave ge wave_g[0] and wave le wave_g[-1] then $
      total_g=linear_interpolate(wave,wave_g,transmission_g)$
    else total_g=0
    if wave ge wave_r[0] and wave le wave_r[-1] then $
      total_r=linear_interpolate(wave,wave_r,transmission_r)$
    else total_r=0
    if wave ge wave_i[0] and wave le wave_i[-1] then $
      total_i=linear_interpolate(wave,wave_i,transmission_i)$
    else total_i=0
    if wave ge wave_z[0] and wave le wave_z[-1] then $
      total_z=linear_interpolate(wave,wave_z,transmission_z)$
    else total_z=0
    total_light=total_g+total_r+total_i+total_z

    for x=0,x_side-1,1 do begin
      for y=0,y_side-1,1 do psf_temp[x,y]=((total_g/total_light)*input_g[x,y])+$
        ((total_r/total_light)*input_r[x,y])+((total_i/total_light)*input_i[x,y])+$
        ((total_z/total_light)*input_z[x,y])
    endfor
    psf_out[*,*,j]=psf_temp

  endfor


  mkhdr,h,psf_out
  sxaddpar,h,'CRPIX1',sxpar(header_g,'CRPIX1')
  sxaddpar,h,'CRPIX2',sxpar(header_g,'CRPIX2')
  sxaddpar,h,'CRVAL1',sxpar(header_g,'CRVAL1')
  sxaddpar,h,'CRVAL2',sxpar(header_g,'CRVAL2')
  sxaddpar,h,'CD1_1',sxpar(header_g,'CD1_1')
  sxaddpar,h,'CD2_2',sxpar(header_g,'CD2_2')
  sxaddpar,h,'CTYPE1',sxpar(header_g,'CTYPE1')
  sxaddpar,h,'CTYPE2',sxpar(header_g,'CTYPE2')
  sxaddpar,h,'CUNIT1',sxpar(header_g,'CUNIT1')
  sxaddpar,h,'CUNIT2',sxpar(header_g,'CUNIT2')
  sxaddpar,h,'IFURA',sxpar(header_g,'IFURA')
  sxaddpar,h,'IFUDEC',sxpar(header_g,'IFUDEC')
  sxaddpar,h,'OBJRA',sxpar(header_g,'OBJRA')
  sxaddpar,h,'OBJDEC',sxpar(header_g,'OBJDEC')
  sxaddpar,h,'BUNIT',sxpar(header_g,'BUNIT')
  sxaddpar,h,'GFWHM',sxpar(header_g,'GFWHM')
  sxaddpar,h,'RFWHM',sxpar(header_r,'RFWHM')
  sxaddpar,h,'IFWHM',sxpar(header_i,'IFWHM')
  sxaddpar,h,'ZFWHM',sxpar(header_z,'ZFWHM')
  sxaddpar,h,'EXPTIME',900



  sxaddpar,h,'CRVAL3',wave_log[0]
  sxaddpar,h,'CDELT3',wave_log[1]-wave_log[0]
  sxaddpar,h,'CD3_3',wave_log[1]-wave_log[0]

  for n=0,sxpar(h,'NAXIS3')-1 do begin
    if total(total(psf_out[*,*,n])) gt 0 then psf_out[*,*,n]=psf_out[*,*,n]/total(psf_out[*,*,n])
  endfor
  fits_write,dir+cube+'_PSF.fits', psf_out, h


end