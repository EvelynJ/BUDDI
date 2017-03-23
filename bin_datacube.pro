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

pro bin_datacube, datacube, setup,$
                  file,start_wavelength, end_wavelength, wavelength_arr, binned_wavelengths,$
                  x_centre,y_centre,GALAXY=galaxy,PSF=psf,MANGA=manga,CALIFA=califa

no_bins=setup.no_bins
root=setup.root
decomp=setup.decomp
binned_dir=setup.binned_dir
slices_dir=setup.slices_dir
median_dir=setup.median_dir
galaxy_ref=setup.galaxy_ref
psf_cube=setup.psf_file
x_centre=fix(setup.x_centre-1)             ;x position of centre of galaxy, -1 to convert to position in array
y_centre=fix(setup.y_centre-1)             ;y position of centre of galaxy, -1 to convert to position in array
sigma_cube=setup.sigma_cube
badpix_cube=setup.badpix_cube
directory=root
no_bins-=3   

;identify if sigma and badpix cubes are provided
result = FILE_TEST(root+sigma_cube+'.fits')
if result eq 1 then sigma_TF='T' else sigma_TF='F'
result = FILE_TEST(root+badpix_cube+'.fits')
if result eq 1 then badpix_TF='T' else badpix_TF='F'

; set number of binned images to 50 as default=> 48 binned images 
; (including one bin going past the upper wavelength limit since 
; it's unlikely the number of images will be exactly divisible by 
; the number of bins) plus the first and last slices. The -3 is to 
; allow for the last bin to go beyond the final wavelength by dividing 
; by one bin less

;read in corrected data cube
if keyword_set(galaxy) then begin
  fits_read,directory+decomp+galaxy_ref+'_smoothed_FLUX.fits', spec_in, header1
  if sigma_TF eq 'T' then fits_read,directory+decomp+galaxy_ref+'_smoothed_SIGMA.fits', sigma_in, h_sig
  if badpix_TF eq 'T' then fits_read,directory+decomp+galaxy_ref+'_smoothed_BADPIX.fits', badpix_in, h_bp
  side1=sxpar(header1,'NAXIS1')
  side2=sxpar(header1,'NAXIS2')
  images=sxpar(header1,'NAXIS3')

  name='image_'
  if sigma_TF eq 'T' then sigma_name='sigma_'
  if badpix_TF eq 'T' then badpix_name='badpix_'

  ;slice up data cube into individual images

  result = FILE_TEST(directory+decomp+slices_dir, /DIRECTORY) 
  if result eq 0 then file_mkdir,directory+decomp+slices_dir
  if result eq 0 then file_mkdir,directory+decomp+slices_dir+'badpix/'
  if result eq 0 then file_mkdir,directory+decomp+slices_dir+'sigma/'

  step=wavelength_arr[1]-wavelength_arr[0]
  
  


  ;create wavelength array and identify elements within wavelength range
  final_image=images-1
  first_image=0
  
  ;read in header again with headfits. this produces a header without 
  ;extensions, which galfit and galfitm can read with no problems
  fits_read,directory+file+'.fits',temp_input,h
;  h = headfits(directory+file+'.fits')
  

  ;print out only the image slices required for final wavelength range to save disc space
  for n=0,images-1,1 do begin
    next_wavelength=wavelength_arr[n]+step
    previous_wavelength=wavelength_arr[n]-step
    ;Identify which image slices are to be 
    ;used based on the wavelength range desired.
    ;*** Note: IDL82+ array[-1]= final element of the array, no longer gives an error message
    if n gt 0 and wavelength_arr[n] ge alog10(start_wavelength) and previous_wavelength lt alog10(start_wavelength) then first_image=n
    if n lt images-5 and wavelength_arr[n] le alog10(end_wavelength) and next_wavelength gt alog10(end_wavelength) then final_image=n
  endfor
 
  ;temp=mrdfits(directory+file+'_FLUX.fits',0,h_temp)
  h_temp=headfits(directory+file+'.fits')
  tempf=mrdfits(directory+file+'.fits',1,h_flux)
  if sigma_TF eq 'T' then temps=mrdfits(directory+sigma_cube+'.fits',1,h_sig)
  if badpix_TF eq 'T' then tempb=mrdfits(directory+badpix_cube+'.fits',1,h_bp)
  
  sxaddpar,h_flux,'NAXIS',2
  sxdelpar,h_flux,'NAXIS3'
  if sigma_TF eq 'T' then sxaddpar,h_sigma,'NAXIS',2
  if sigma_TF eq 'T' then sxdelpar,h_sigma,'NAXIS3'
  if badpix_TF eq 'T' then sxaddpar,h_bp,'NAXIS',2
  if badpix_TF eq 'T' then sxdelpar,h_bp,'NAXIS3'


  sxaddpar,h_temp,'EXPTIME',sxpar(h,'EXPTIME')/sxpar(h,'NEXP')
  sxaddpar,h_temp,'NCOMBINE',sxpar(h,'NEXP')

  sxaddpar,h_flux,'EXPTIME',sxpar(h,'EXPTIME')/sxpar(h,'NEXP')
  sxaddpar,h_flux,'NCOMBINE',sxpar(h,'NEXP')

  if sigma_TF eq 'T' then begin
    sxaddpar,h_sig,'EXPTIME',sxpar(h,'EXPTIME')/sxpar(h,'NEXP')
    sxaddpar,h_sig,'NCOMBINE',sxpar(h,'NEXP')
  endif

  if badpix_TF eq 'T' then begin
    sxaddpar,h_bp,'EXPTIME',sxpar(h,'EXPTIME')/sxpar(h,'NEXP')
    sxaddpar,h_bp,'NCOMBINE',sxpar(h,'NEXP')
  endif
  
  
  ;print out slices for wavelength range
  for n=first_image,final_image,1 do begin
    mkhdr,hdr0,spec_in[*,*,n]
  
    sxaddpar,h_temp,'Wavelength',10^(wavelength_arr[n])
    sxaddpar,h_temp,'CRVAL3',wavelength_arr[n]
    sxaddpar,hdr0,'Wavelength',10^(wavelength_arr[n])
    sxaddpar,hdr0,'CRVAL3',wavelength_arr[n]
    sxaddpar,h_flux,'Wavelength',10^(wavelength_arr[n])
    sxaddpar,h_flux,'CRVAL3',wavelength_arr[n]
    fits_write,directory+decomp+slices_dir+name+string(n,format='(i4.4)')+'.fits',spec_in[*,*,n],hdr0
;    fits_write,directory+decomp+slices_dir+name+string(n,format='(i4.4)')+'.fits',spec_in[*,*,n],extname='FLUX'
;    modfits,directory+decomp+slices_dir+name+string(n,format='(i4.4)')+'.fits',0,h_temp
;    modfits,directory+decomp+slices_dir+name+string(n,format='(i4.4)')+'.fits',1,h_flux,extname='FLUX'
    
    if sigma_TF eq 'T' then begin
      sxaddpar,h_sig,'Wavelength',10^(wavelength_arr[n])
      sxaddpar,h_sig,'CRVAL3',wavelength_arr[n]
      fits_write,directory+decomp+slices_dir+'sigma/'+sigma_name+string(n,format='(i4.4)')+'.fits',sigma_in[*,*,n],hdr0
;      fits_write,directory+decomp+slices_dir+sigma_name+string(n,format='(i4.4)')+'.fits',sigma_in[*,*,n],extname='SIGMA'
;      modfits,directory+decomp+slices_dir+sigma_name+string(n,format='(i4.4)')+'.fits',0,h_temp
;      modfits,directory+decomp+slices_dir+sigma_name+string(n,format='(i4.4)')+'.fits',1,h_sig,extname='SIGMA'
     endif
    
    if badpix_TF eq 'T' then begin
      sxaddpar,h_bp,'Wavelength',10^(wavelength_arr[n])
      sxaddpar,h_bp,'CRVAL3',wavelength_arr[n]
      fits_write,directory+decomp+slices_dir+'badpix/'+badpix_name+string(n,format='(i4.4)')+'.fits',badpix_in[*,*,n],hdr0
;      fits_write,directory+decomp+slices_dir+badpix_name+string(n,format='(i4.4)')+'.fits',badpix_in[*,*,n],extname='BADPIX'
;      modfits,directory+decomp+slices_dir+badpix_name+string(n,format='(i4.4)')+'.fits',0,h_temp
;      modfits,directory+decomp+slices_dir+badpix_name+string(n,format='(i4.4)')+'.fits',1,h_bp,extname='BADPIX'
    endif
   
    
   endfor

  openw,45,directory+decomp+slices_dir+'info.txt'
  printf,45,'First_image_slice    ',first_image
  printf,45,'Last_image_slice  ',final_image
  printf,45,'No_of_bins          ',no_bins+3
  printf,45,'No_of_images_per_bin ',long((final_image-first_image-1)/no_bins)
  printf,45,'Start_wavelength     ',10^(wavelength_arr[first_image])
  printf,45,'End_wavelength     ',10^(wavelength_arr[final_image])
  close,45







  ;co-add slices to make binned images with better S/N. 
  ;First need to identify whether the number of images 
  ;required is exactly divisible by the number of bins
  final_image=images-1
  first_image=0
  total_no_images=final_image-first_image-1         ;total number of images to be used in binning
  divisible=total_no_images mod no_bins         ;check to see if the number of images will divide cleanly into the bins
  bins_no_images=long(total_no_images/no_bins)  ;identify the number of images to be included in each bin

  result = FILE_TEST(directory+decomp+binned_dir, /DIRECTORY) 
  if result eq 0 then file_mkdir,directory+decomp+binned_dir
  
  binned_image=fltarr(side1,side2)
  binned_sigma=fltarr(side1,side2)
  binned_badpix=fltarr(side1,side2)
  binned_wavelengths=fltarr(no_bins+3)
  
  
  ;print first image slice
  binned_image=spec_in[*,*,first_image-1]
  if badpix_TF eq 'T' then binned_badpix=badpix_in[*,*,first_image-1]
  if sigma_TF eq 'T' then binned_sigma=sigma_in[*,*,first_image-1]

    
  mkhdr,hdr0,spec_in[*,*,0]
  wavelength=10^(wavelength_arr[0]);sxpar(header1,'WAVELENG')
  sxaddpar,hdr0,'Wavelength',wavelength
  sxaddpar,h_temp,'Wavelength',wavelength
  sxaddpar,h_flux,'Wavelength',wavelength
  fits_write,directory+decomp+binned_dir+name+string(0,format='(i4.4)')+'.fits', spec_in[*,*,0],hdr0
;  fits_write,directory+decomp+binned_dir+name+string(0,format='(i4.4)')+'.fits', spec_in[*,*,0],extname='FLUX'
;  modfits,directory+decomp+binned_dir+name+string(0,format='(i4.4)')+'.fits',0,h_temp
;  modfits,directory+decomp+binned_dir+name+string(0,format='(i4.4)')+'.fits',1,h_flux,extname='FLUX'


  
  if sigma_TF eq 'T' then begin
    sxaddpar,h_temp,'Wavelength',wavelength
    sxaddpar,h_sig,'Wavelength',wavelength
    fits_write,directory+decomp+binned_dir+sigma_name+string(0,format='(i4.4)')+'.fits', sigma_in[*,*,0],hdr0
;    fits_write,directory+decomp+binned_dir+sigma_name+string(0,format='(i4.4)')+'.fits', sigma_in[*,*,0],extname='SIGMA'
;    modfits,directory+decomp+binned_dir+sigma_name+string(0,format='(i4.4)')+'.fits',0,h_temp
;    modfits,directory+decomp+binned_dir+sigma_name+string(0,format='(i4.4)')+'.fits',1,h_sig,extname='SIGMA'
  endif
  if badpix_TF eq 'T' then begin
    sxaddpar,h_temp,'Wavelength',wavelength
    sxaddpar,h_bp,'Wavelength',wavelength
    fits_write,directory+decomp+binned_dir+badpix_name+string(0,format='(i4.4)')+'.fits', badpix_in[*,*,0],hdr0
;    fits_write,directory+decomp+binned_dir+badpix_name+string(0,format='(i4.4)')+'.fits', badpix_in[*,*,0],extname='BADPIX'
;    modfits,directory+decomp+binned_dir+badpix_name+string(0,format='(i4.4)')+'.fits',0,h_temp
;    modfits,directory+decomp+binned_dir+badpix_name+string(0,format='(i4.4)')+'.fits',1,h_bp,extname='BADPIX'
  endif
  
  ;modfits,directory+decomp+binned_dir+name+string(0,format='(i4.4)')+'.fits',0,hdr0
  wavelength1=wavelength
  binned_wavelengths[0]=wavelength
  
  odd_even=bins_no_images mod 2
  
  
  ;bin remaining image slices
  for run=0,no_bins,1 do begin
      a=first_image+(run*bins_no_images)
      if run ne no_bins then b=first_image+((run+1)*bins_no_images)-1 $
        else b=final_image-1
        
      image=spec_in[*,*,a:b]
      if sigma_TF eq 'T' then sig=sigma_in[*,*,a:b]
      if badpix_TF eq 'T' then bp=badpix_in[*,*,a:b]
     ;if run eq no_bins then final_image=first_image+((run+1)*bins_no_images)
      for x=0,side1-1,1 do begin
        for y=0,side2-1,1 do begin
            ;binned_image[x,y]=total(image[x,y,*])
            
            RESISTANT_Mean,image[x,y,*],3,mean_temp     ;don't include pixels with values >3sigma from mean
            binned_image[x,y]=mean_temp
            
            if sigma_TF eq 'T' then begin
              temp_val=where(finite(sig[x,y,*]))
              if n_elements(temp_val) eq 1 then temp_val=[0,1,2]
              if n_elements(temp_val) gt 1 then RESISTANT_Mean,sig[x,y,temp_val],3,mean_temp     ;don't include pixels with values >3sigma from mean
              binned_sigma[x,y]=mean_temp
            endif
            
            if badpix_TF eq 'T' then begin
              RESISTANT_Mean,bp[x,y,*],3,mean_temp     ;don't include pixels with values >3sigma from mean
              binned_badpix[x,y]=mean_temp
            endif
        endfor
      endfor
      if odd_even eq 0 then element=first_image+run*(bins_no_images/2)+(step/2)
      if odd_even ne 0 then element=first_image+run*((bins_no_images+1)/2)
 
      wavelength=10^(wavelength_arr[fix(0.5*(a+b))])
      sxaddpar,hdr0,'Wavelength',wavelength
      sxaddpar,h_temp,'Wavelength',wavelength
      sxaddpar,h_flux,'Wavelength',wavelength
      binned_wavelengths[run+1]=wavelength
      fits_write,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits', binned_image,hdr0
;      fits_write,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits', binned_image,extname='FLUX'
;      modfits,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits',0,h_temp
;      modfits,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits',1,h_flux,extname='FLUX'
 
      if sigma_TF eq 'T' then begin
        
        sxaddpar,h_sigma,'Wavelength',wavelength
        fits_write,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits', binned_sigma,hdr0
;        fits_write,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits', binned_sigma,extname='SIGMA'
;        modfits,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits',0,h_temp
;        modfits,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits',1,h_sig,extname='SIGMA'
      endif
      
      if badpix_TF eq 'T' then begin
        sxaddpar,h_bp,'Wavelength',wavelength
        fits_write,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits', binned_badpix,hdr0
;        fits_write,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits', binned_badpix,extname='BADPIX'
;        modfits,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits',0,h_temp
;        modfits,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits',1,h_bp,extname='BADPIX'
      endif
      
 endfor  

  ;finally print last image slice
  binned_image=spec_in[*,*,-1]
  if badpix_TF eq 'T' then binned_badpix=badpix_in[*,*,-1]
  if sigma_TF eq 'T' then binned_sigma=sigma_in[*,*,-1]
  
  wavelength=10^(wavelength_arr[-1])
  sxaddpar,hdr0,'Wavelength',wavelength
  sxaddpar,h_temp,'Wavelength',wavelength
  sxaddpar,h_flux,'Wavelength',wavelength
  binned_wavelengths[run+1]=wavelength
  
  fits_write,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits', binned_image,hdr0
;  fits_write,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits', binned_image,extname='FLUX'
;  modfits,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits',0,h_temp
;  modfits,directory+decomp+binned_dir+name+string(run+1,format='(i4.4)')+'.fits',1,h_flux,extname='FLUX'

  if sigma_TF eq 'T' then begin
    sxaddpar,h_sigma,'Wavelength',wavelength
    fits_write,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits', binned_sigma,hdr0
;    fits_write,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits', binned_sigma,extname='SIGMA'
;    modfits,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits',0,h_temp
;    modfits,directory+decomp+binned_dir+sigma_name+string(run+1,format='(i4.4)')+'.fits',1,h_sig,extname='SIGMA'
  endif

  if badpix_TF eq 'T' then begin
    sxaddpar,h_bp,'Wavelength',wavelength
    fits_write,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits', binned_badpix,hdr0
;    fits_write,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits', binned_badpix,extname='BADPIX';, h
;    modfits,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits',0,h_temp
;    modfits,directory+decomp+binned_dir+badpix_name+string(run+1,format='(i4.4)')+'.fits',1,h_bp,extname='BADPIX'
  endif
  
  wavelength2=sxpar(h,'Waveleng')
  
  
  
  openw,45,directory+decomp+binned_dir+'info.txt'
  printf,45,'First_image_slice    ',first_image
  printf,45,'Last_image_slice  ',final_image
  printf,45,'No_of_bins          ',no_bins+3
  printf,45,'No_of_images_per_bin ',bins_no_images
  printf,45,'Start_wavelength     ',wavelength1
  printf,45,'End_wavelength     ',wavelength_arr[-1]
  close,45
  
  
  
endif






;repeat for median image
result = FILE_TEST(directory+decomp+median_dir, /DIRECTORY) 
if result eq 0 then file_mkdir,directory+decomp+median_dir

wavelength=10^(median(wavelength_arr[first_image:final_image]))
sxaddpar,h_temp,'Wavelength',wavelength
sxaddpar,h_flux,'Wavelength',wavelength
sxaddpar,h_sigma,'Wavelength',wavelength
sxaddpar,h_bp,'Wavelength',wavelength


for x=0,side1-1,1 do begin
  for y=0,side2-1,1 do begin
    RESISTANT_Mean,spec_in[x,y,first_image:final_image],3,mean_temp     ;don't include pixels with values >3sigma from mean
    binned_image[x,y]=mean_temp
  endfor
endfor
mkhdr,hdr0,binned_image
sxaddpar,hdr0,'Wavelength',wavelength
fits_write,directory+decomp+median_dir+'image.fits', binned_image,hdr0;,extname='FLUX'
;fits_write,directory+decomp+median_dir+'image.fits', binned_image,extname='FLUX'
;modfits,directory+decomp+median_dir+'image.fits',0,h_temp
;modfits,directory+decomp+median_dir+'image.fits',1,h_flux,extname='FLUX'


if badpix_TF eq 'T' then begin
  for x=0,side1-1,1 do begin
    for y=0,side2-1,1 do begin
      RESISTANT_Mean,badpix_in[x,y,first_image:final_image],3,mean_temp     ;don't include pixels with values >3sigma from mean
      binned_badpix[x,y]=mean_temp
    endfor
  endfor
  fits_write,directory+decomp+median_dir+'badpix.fits', binned_badpix,hdr0;,extname='BADPIX'
  ;fits_write,directory+decomp+median_dir+'badpix.fits', binned_badpix,extname='BADPIX'
  ;modfits,directory+decomp+median_dir+'badpix.fits',0,h_temp
  ;modfits,directory+decomp+median_dir+'badpix.fits',1,h_bp,extname='BADPIX'
endif
if sigma_TF eq 'T' then begin
  for x=0,side1-1,1 do begin
    for y=0,side2-1,1 do begin
      RESISTANT_Mean,sigma_in[x,y,first_image:final_image],3,mean_temp     ;don't include pixels with values >3sigma from mean
      binned_sigma[x,y]=mean_temp
    endfor
  endfor
  fits_write,directory+decomp+median_dir+'sigma.fits', binned_sigma,hdr0;,extname='SIGMA'
  ;fits_write,directory+decomp+median_dir+'sigma.fits', binned_sigma,extname='SIGMA'
  ;modfits,directory+decomp+median_dir+'sigma.fits',0,h_temp
  ;modfits,directory+decomp+median_dir+'sigma.fits',1,h_sig,extname='SIGMA'
endif







;read in PSF datacube and split into image slices and binned PSF images
;repeat for bad pixel mask if included as datacube

result = FILE_TEST(directory+psf_cube) 
if result eq 1 then begin
  fits_read,directory+psf_cube,psf_cube_input,h_psf
  sxdelpar,h_psf,'CRVAL3'
  sxdelpar,h_psf,'CDELT3'
  sxdelpar,h_psf,'CD3_3'
  
  ;***median PSF image***
  fits_read,directory+file+'.fits',temp_input,h
  combined_psf=fltarr(sxpar(h_psf,'NAXIS1'),sxpar(h_psf,'NAXIS2'))
  for column=0,sxpar(h_psf,'NAXIS1')-1,1 do begin
    for row=0,sxpar(h_psf,'NAXIS2')-1,1 do begin
      combined_psf[column,row]=mean(psf_cube_input[column,row,*])
    endfor
  endfor
  
  fits_write,directory+decomp+median_dir+'psf.fits',combined_psf,h_psf



  ;***slices PSF images***
  file_mkdir,directory+decomp+slices_dir+'PSF/'
    
  readcol,directory+decomp+slices_dir+'info.txt',format='x,f',param,comment='#',/SILENT
  for j=param[0],param[1],1 do begin
    if total(psf_cube_input[*,*,j]) ne 0 then $
      fits_write,directory+decomp+slices_dir+'PSF/'+string(j,format='(i4.4)')+'.fits', psf_cube_input[*,*,j], h_psf $
    else fits_write,directory+decomp+slices_dir+'PSF/'+string(j,format='(i4.4)')+'.fits', combined_psf, h_psf
  endfor
  
  ;***binned PSF images***
  for j=0,no_bins+2,1 do begin
    ;note- sometimes the first and last PSF images will contain no information.
    ;this will lead to a segmentation faiult in GALFIT.
    ;solution-replace with the valid PSF image before or after
    h = headfits(directory+decomp+binned_dir+name+string(j,format='(i4.4)')+'.fits') 
    wave=sxpar(h,'WAVELENG')
    
    ;identify slice within psf_datacube with the wavelength closest to 
    ;central wavelength of that bin, and use that as PSF
    element=where(10^(wavelength_arr) ge wave)
    sxaddpar,h_psf,'Wavelength',wave
    result = FILE_TEST(directory+decomp+binned_dir+'PSF/', /DIRECTORY) 
    if result eq 0 then file_mkdir,directory+decomp+binned_dir+'PSF/'
    
    if total(psf_cube_input[*,*,j]) ne 0 then $
      fits_write,directory+decomp+binned_dir+'PSF/'+string(j,format='(i4.4)')+'.fits', psf_cube_input[*,*,element[0]], h_psf $
    else fits_write,directory+decomp+binned_dir+'PSF/'+string(j,format='(i4.4)')+'.fits', combined_psf, h_psf
  endfor  
endif






delvarx,psf_cube,combined_psf,psf_cube_input,spec_in,binned_image,binned_badpix,binned_sigma,tempf,tempb,temps





end