; this code will identify bad pixels in the MaNGA data cube, 
; which for now are those with a value of 0 outside of the 
; MaNGA field of view
; 
; Modified to create separate masks for the first and last 
; images where all pixels are masked. This will reduce 
; mis-fits due to the different binning techniques for these 
; images (Nov 2014)
; 

pro badpixelmask, setup, badpix
  root=setup.root
  decomp=setup.decomp
  binned_dir=setup.binned_dir
  slices_dir=setup.slices_dir
  median_dir=setup.median_dir
  galaxy_ref=setup.galaxy_ref
  badpix_cube=setup.badpix_cube  ;badpixel cube
  badpix_file=setup.badpix_file  ;bad pixel file for foreground stars

;if no badpixel cube preset, run as before by masking 0-value pixels
;if badpixel cube present, ignore manual mask and simply add foreground stars to array

result = FILE_TEST(root+badpix_cube+'.fits')
if result eq 1 then badpix_TF='T' else badpix_TF='F'

if badpix_TF eq 'F' then begin
  directory=root+decomp
  fits_read,directory+galaxy_ref+'_smoothed_flux.fits', spec_in, header1
  side1=sxpar(header1,'NAXIS1')
  side2=sxpar(header1,'NAXIS2')
  images=sxpar(header1,'NAXIS3')
  
  badpix=intarr(side1,side2)
  badpix_end=intarr(side1,side2)
  
  image_slice=spec_in[*,*,fix(images/2)]

  ;Pick a random image slics, and identify bad pixels 
  ;as those with a pixel value of 0.0
  for x=0,side1-1,1 do begin
    for y=0,side2-1,1 do begin
      badpix_end[x,y]=1
      if image_slice[x,y] eq 0.0 then badpix[x,y]=1
      if x eq 0 or y eq 0 then badpix[x,y]=1
      endif 
    endfor
  endfor
endif else begin
  fits_read,directory+galaxy_ref+'_smoothed_badpix.fits', badpix, header1
  badpix_end=intarr(side1,side2)
  badpix_end[*,*]=1
endelse


if file_test(root+badpix_file+'.fits') eq 1 then begin
  ;print,'***YESYESYESYESYES***'
 
  readcol,root+badpix_file,format='f,f',x_bad,y_bad,comment='#',/SILENT
  if badpix_TF eq 'F' then begin
    for j=0,n_elements(x_bad)-1,1 do badpix[x_bad,y_bad]=1
  endif
  if badpix_TF eq 'T' then begin
    for j=0,n_elements(x_bad)-1,1 do badpix[x_bad,y_bad,*]=1
  endif
endif
fits_write,directory+median_dir+'badpix.fits',byte(badpix),extname='BADPIX'
;temp=mrdfits(directory+file+'_BADPIX.fits',0,h_0)
h_temp=headfits(directory+file+'_BADPIX.fits')
modfits,directory+binned_dir+'badpix_end.fits',0,h_0
modfits,directory+median_dir+'badpix.fits',0,h_0
if badpix_TF eq 'T' then begin 
  tempb=mrdfits(directory+file+'_BADPIX.fits',1,h_bp)
  modfits,directory+binned_dir+'badpix_end.fits',1,h_bp,extname='BADPIX'
  modfits,directory+median_dir+'badpix.fits',1,h_bp,extname='BADPIX'
endif

;print out bad pixel masks if no input datacube used
if badpix_TF eq 'F' then begin
  fits_write,directory+binned_dir+'badpix.fits',byte(badpix)
  fits_write,directory+slices_dir+'badpix.fits',byte(badpix)
  modfits,directory+binned_dir+'badpix_end.fits',0,h_0
  modfits,directory+median_dir+'badpix.fits',0,h_0
endif

;print out bad pixel masks when input datacube provided
if badpix_TF eq 'T' then begin
  s=size(badpix)
  for n=0,s[3]-1,1 do begin
    fits_write,directory+binned_dir+'badpix_'+string(n,format='(i4.4)')+'.fits',byte(badpix),extname='BADPIX'
    fits_write,directory+slices_dir+'badpix_'+string(n,format='(i4.4)')+'.fits',byte(badpix),extname='BADPIX'
    modfits,directory+binned_dir+'badpix_'+string(n,format='(i4.4)')+'.fits',0,h_0
    modfits,directory+binned_dir+'badpix_'+string(n,format='(i4.4)')+'.fits',1,h_bp,extname='BADPIX'
    modfits,directory+slices_dir+'badpix_'+string(n,format='(i4.4)')+'.fits',0,h_0
    modfits,directory+slices_dir+'badpix_'+string(n,format='(i4.4)')+'.fits',1,h_bp,extname='BADPIX'
  endfor
endif

delvarx,badpix,badpix_end
end