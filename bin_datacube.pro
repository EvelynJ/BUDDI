;
; 
; This code slices the wavelength and kinematics corrected data cube into image , 
; slices, and coadds them to produce higher S/N images covering the full spectral 
; range.
; 
; Modified to bin the data using the resistant mean. i.e. for each pixel position, 
; the mean value is measured over the range of image slices used by excluding those
; values >3sigma from the mean value, and then multiplied by the number of image 
; slices to get the total value. This technique reduces the effects of intermittent 
; effects, such as cosmic rays or dead pixels.
; 
; Note- this version of the code makes a PSF file using the FWHM information in 
; the header
; 

pro bin_datacube, datacube, no_bins, root,dir,decomp,binned_dir,slices_dir,median_dir,$
                  galaxy_ref,file,start_wavelength, end_wavelength, wavelength_arr, binned_wavelengths,$
                  x_centre,y_centre,GALAXY=galaxy,PSF=psf,MANGA=manga,CALIFA=califa

directory=root+dir
no_bins-=3   
; set number of binned images to 50 as default=> 48 binned images 
; (including one bin going past the upper wavelength limit since 
; it's unlikely the number of images will be exactly divisible by 
; the number of bins) plus the first and last slices. The -3 is to 
; allow for the last bin to go beyond the final wavelength by dividing 
; by one bin less

;read in corrected data cube
if keyword_set(galaxy) then begin
  fits_read,directory+decomp+galaxy_ref+'_smoothed_kinematics.fits', spec_in, header1
  side1=sxpar(header1,'NAXIS1')
  side2=sxpar(header1,'NAXIS2')
  images=sxpar(header1,'NAXIS3')

  name='image_'

  ;slice up data cube into individual images

  result = FILE_TEST(directory+decomp+slices_dir, /DIRECTORY) 
  if result eq 0 then file_mkdir,directory+decomp+slices_dir

;  wavelength=fltarr(images)
;  wavelength_selected=fltarr(images)
;  wavelength0=sxpar(header1,'CRVAL3')
  step=wavelength_arr[1]-wavelength_arr[0]
  
  


  ;create wavelength array and identify elements within wavelength range
  final_image=images-1
  first_image=images-1
  
  ;read in header again with headfits. this produces a header without 
  ;extensions, which galfit and galfitm can read with no problems
  h = headfits(directory+file+'.fits')
  for n=0,images-1,1 do begin
    next_wavelength=wavelength_arr[n]+step
    previous_wavelength=wavelength_arr[n]-step
    ;Identify which image slices are to be 
    ;used based on the wavelength range desired.
    ;*** Note: IDL82+ array[-1]= final element of the array, no longer gives an error message
    if keyword_set(manga) then begin
      if n gt 0 and wavelength_arr[n] ge alog10(start_wavelength) and previous_wavelength lt alog10(start_wavelength) then first_image=n
      if n lt images-5 and wavelength_arr[n] le alog10(end_wavelength) and next_wavelength gt alog10(end_wavelength) then final_image=n
      
      sxaddpar,h,'Wavelength',10^(wavelength_arr[n])
      sxaddpar,h,'CRVAL3',wavelength_arr[n]
    endif else if keyword_set(califa) then begin
      if n gt 0 and wavelength_arr[n] ge alog(start_wavelength) and previous_wavelength lt alog(start_wavelength) then first_image=n
      if n lt images-5 and wavelength_arr[n] le alog(end_wavelength) and next_wavelength gt alog(end_wavelength) then final_image=n
      sxaddpar,h,'Wavelength',exp(wavelength_arr[n])
      sxaddpar,h,'CRVAL3',wavelength_arr[n]
      sxaddpar,h,'CD3_3',step
      sxaddpar,h,'CDELT3',step
    endif
    if n ge first_image and n le final_image then fits_write,directory+decomp+slices_dir+name+string(n,format='(i4.4)')+'.fits',spec_in[*,*,n], h
;    wavelength_selected
  endfor


  ;co-add slices to make binned images with better S/N. 
  ;First need to identify whether the number of images 
  ;required is exactly divisible by the number of bins
  total_no_images=final_image-first_image+1         ;total number of images to be used in binning
  divisible=total_no_images mod no_bins         ;check to see if the number of images will divide cleanly into the bins
  bins_no_images=long(total_no_images/no_bins)  ;identify the number of images to be included in each bin

  result = FILE_TEST(directory+decomp+binned_dir, /DIRECTORY) 
  if result eq 0 then file_mkdir,directory+decomp+binned_dir
  
  binned_image=fltarr(side1,side2)
  binned_wavelengths=fltarr(no_bins+3)
  
  
  ;print first image slice
  binned_image=spec_in[*,*,first_image-1]
  h = headfits(directory+decomp+slices_dir+name+string(first_image,format='(i4.4)')+'.fits')
  ;sxaddpar,h,'Wavelength',wavelength[first_image-1]
  wavelength=sxpar(h,'WAVELENG')
  sxaddpar,h,'Wavelength',wavelength
  fits_write,directory+decomp+binned_dir+name+string(0,format='(i4.4)')+'.fits', binned_image, h
  wavelength1=wavelength
  binned_wavelengths[0]=wavelength
  
  odd_even=bins_no_images mod 2
  
  
  ;bin remaining image slices
  for run=0,no_bins,1 do begin
      a=first_image+(run*bins_no_images)
      if run ne no_bins then b=first_image+((run+1)*bins_no_images)-1 $
        else b=final_image-1
        
      image=spec_in[*,*,a:b]
      ;if run eq no_bins then final_image=first_image+((run+1)*bins_no_images)
      for x=0,side1-1,1 do begin
        for y=0,side2-1,1 do begin
            ;binned_image[x,y]=total(image[x,y,*])
            RESISTANT_Mean,image[x,y,*],3,mean_temp     ;don't include pixels with values >3sigma from mean
            binned_image[x,y]=mean_temp;*bins_no_images
        endfor
      endfor
      if odd_even eq 0 then element=first_image+run*(bins_no_images/2)+(step/2)
      if odd_even ne 0 then element=first_image+run*((bins_no_images+1)/2)
      ;h = headfits(directory+galaxy_ref+'-LOGCUBE.fits')
      ;sxaddpar,h,'Wavelength',wavelength[element]
      h = headfits(directory+decomp+slices_dir+name+string(0.5*(a+b),format='(i4.4)')+'.fits')
      wavelength=sxpar(h,'WAVELENG')
      sxaddpar,h,'Wavelength',wavelength
      binned_wavelengths[run+1]=wavelength
      ;print,run,first_image+(run*bins_no_images),first_image+((run+1)*bins_no_images)-1,wavelength
      fits_write,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits', binned_image, h
  endfor  

  ;finally print last image slice
  binned_image=spec_in[*,*,final_image]
  h = headfits(directory+decomp+slices_dir+name+string(final_image,format='(i4.4)')+'.fits')
  wavelength=sxpar(h,'WAVELENG')
  sxaddpar,h,'Wavelength',wavelength
  binned_wavelengths[run+1]=wavelength
  fits_write,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits', binned_image, h
  
  wavelength2=sxpar(h,'Waveleng')
  
  
  
  openw,45,directory+decomp+binned_dir+'info.txt'
  printf,45,'First_image_slice    ',first_image
  printf,45,'Last_image_slice  ',final_image
  printf,45,'No_of_bins          ',no_bins+3
  printf,45,'No_of_images_per_bin ',bins_no_images
  printf,45,'Start_wavelength     ',wavelength1
  printf,45,'End_wavelength     ',wavelength2
  close,45
  
  openw,45,directory+decomp+slices_dir+'info.txt'
  printf,45,'First_image_slice    ',first_image
  printf,45,'Last_image_slice  ',final_image
  printf,45,'No_of_bins          ',no_bins+3
  printf,45,'No_of_images_per_bin ',bins_no_images
  printf,45,'Start_wavelength     ',wavelength1
  printf,45,'End_wavelength     ',wavelength2
  close,45
  
  
endif



;Now create PSF for the decomposition. best way to do this for SDSS MaNGA
;data is to read in the FWHM for each waveband, and then interpolate 
;between their central wavelengths.
;
;read in PSF extensions for each filter

;fits_read,directory+file+'.fits',input_u,header_u,extname='UPSF'
fits_read,directory+file+'.fits',input_g,header_g,extname='GPSF'
fits_read,directory+file+'.fits',input_r,header_r,extname='RPSF'
fits_read,directory+file+'.fits',input_i,header_i,extname='IPSF'
fits_read,directory+file+'.fits',input_z,header_z,extname='ZPSF'

h = headfits(directory+file+'.fits')
sxaddpar,h,'EXPTIME',1
sxaddpar,h,'GAIN',1
x_side=sxpar(header_g,'NAXIS1')
y_side=sxpar(header_g,'NAXIS2')
ugriz_wave=[3543,4770,6231,7635,9134]

;for each image slice, use the wavelength anf the filter transmission curves 
;to determine the fraction of light from each filter and coadd them in the 
;correct proportions.
;
;readcol,root+'u.dat',format='F,X,F,X,X',wave_u,transmission_u
readcol,root+'g.dat',format='F,X,F,X,X',wave_g,transmission_g
readcol,root+'r.dat',format='F,X,F,X,X',wave_r,transmission_r
readcol,root+'i.dat',format='F,X,F,X,X',wave_i,transmission_i
readcol,root+'z.dat',format='F,X,F,X,X',wave_z,transmission_z

for j=0,images-1,1 do begin
  wave=wavelength_arr[j]
  total_light=0
;  if wave ge wave_u[0] and wave le wave_u[-1] then $
;    total_u=linear_interpolate(wave,wave_u,transmission_u)
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
  
  psf_out=fltarr(x_side,y_side)
  for x=0,x_side-1,1 do begin
    for y=0,y_side-1,1 do psf_out[x,y]=((total_g/total_light)*input_g[x,y])+$
      ((total_r/total_light)*input_r[x,y])+((total_i/total_light)*input_i[x,y])+$
      ((total_z/total_light)*input_z[x,y])
  endfor
  
  sxaddpar,header_g,'Wavelength',wave
  result = FILE_TEST(directory+decomp+slices_dir+'PSF/', /DIRECTORY) 
  if result eq 0 then file_mkdir,directory+decomp+slices_dir+'PSF/'
  fits_write,directory+decomp+slices_dir+'PSF/'+string(j,format='(i4.4)')+'.fits', psf_out, header_u
endfor


;PSF for binned images
for j=0,no_bins+2,1 do begin
  h = headfits(directory+decomp+binned_dir+name+string(j,format='(i4.4)')+'.fits')
  wave=sxpar(h,'WAVELENG')
  
  total_light=0
;  if wave ge wave_u[0] and wave le wave_u[-1] then $
;    total_u=linear_interpolate(wave,wave_u,transmission_u)
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
  
  psf_out=fltarr(x_side,y_side)
  for x=0,x_side-1,1 do begin
    for y=0,y_side-1,1 do psf_out[x,y]=((total_g/total_light)*input_g[x,y])+$
      ((total_r/total_light)*input_r[x,y])+((total_i/total_light)*input_i[x,y])+$
      ((total_z/total_light)*input_z[x,y])
  endfor
  
  sxaddpar,header_g,'Wavelength',wave
  result = FILE_TEST(directory+decomp+binned_dir+'PSF/', /DIRECTORY) 
  if result eq 0 then file_mkdir,directory+decomp+binned_dir+'PSF/'
  fits_write,directory+decomp+binned_dir+'PSF/'+string(j,format='(i4.4)')+'.fits', psf_out, header_u
endfor


;
;fits_write,directory+decomp+slices_dir+'Upsf.fits',input_u,h;header_IFU
;fits_write,directory+decomp+slices_dir+'Gpsf.fits',input_g,h;header_IFU
;fits_write,directory+decomp+slices_dir+'Rpsf.fits',input_r,h;header_IFU
;fits_write,directory+decomp+slices_dir+'Ipsf.fits',input_i,h;header_IFU
;fits_write,directory+decomp+slices_dir+'Zpsf.fits',input_z,h;header_IFU

;fits_write,directory+decomp+binned_dir+'Upsf.fits',input_u,h;header_IFU
fits_write,directory+decomp+binned_dir+'Gpsf.fits',input_g,h;header_IFU
fits_write,directory+decomp+binned_dir+'Rpsf.fits',input_r,h;header_IFU
fits_write,directory+decomp+binned_dir+'Ipsf.fits',input_i,h;header_IFU
fits_write,directory+decomp+binned_dir+'Zpsf.fits',input_z,h;header_IFU

combined_psf=fltarr(x_side,y_side)
for column=0,x_side-1,1 do begin
  for row=0,y_side-1,1 do begin
    combined_psf[column,row]=mean([input_g[column,row],input_r[column,row],input_i[column,row],input_z[column,row]])
  endfor
endfor

result = FILE_TEST(directory+decomp+median_dir, /DIRECTORY) 
if result eq 0 then file_mkdir,directory+decomp+median_dir

fits_write,directory+decomp+median_dir+'psf.fits',combined_psf,h;header_IFU


end