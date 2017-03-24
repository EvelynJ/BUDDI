; this code will identify bad pixels in the MaNGA data cube, 
; which for now are those with a value of 0 outside of the 
; MaNGA field of view
; 
; Modified to create separate masks for the first and last 
; images where all pixels are masked. This will reduce 
; mis-fits due to the different binning techniques for these 
; images (Nov 2014)
; 
; If a badpixel cube is provided, bin_datacube creates correct masks.
; This code will only create bad pixel masks if no input datacube is 
; provided, or if a text file with additional masks is provided
; 

pro badpixelmask, setup
  root=setup.root
  decomp=setup.decomp
  binned_dir=setup.binned_dir
  slices_dir=setup.slices_dir
  median_dir=setup.median_dir
  galaxy_ref=setup.galaxy_ref
  badpix_cube=setup.badpix_cube  ;badpixel cube
  badpix_file=setup.badpix_file  ;bad pixel file for foreground stars
  file=setup.file
  
  
;if no badpixel cube preset, run as before by masking 0-value pixels
;if badpixel cube present, ignore manual mask and simply add foreground stars to array

result = FILE_TEST(root+badpix_cube+'.fits')
if result eq 1 then badpix_TF='T' else badpix_TF='F'
directory=root+decomp

if badpix_TF eq 'F' then begin
  fits_read,directory+galaxy_ref+'_smoothed_FLUX.fits', spec_in, header1
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
    endfor
  endfor
  
  ;add positions of stars to mask if provided
  if file_test(root+badpix_file) eq 1 then begin
    readcol,root+badpix_file,format='f,f',x_bad,y_bad,comment='#',/SILENT
    for j=0,n_elements(x_bad)-1,1 do badpix[x_bad,y_bad]=1
  endif
  mkhdr,hdr0,badpix
  fits_write,directory+median_dir+'badpix.fits',byte(badpix),hdr0;extname='BADPIX'
  fits_write,directory+binned_dir+'badpix.fits',byte(badpix),hdr0;extname='BADPIX'
  fits_write,directory+binned_dir+'badpix_end.fits',byte(badpix_end),hdr0;extname='BADPIX'
  fits_write,directory+slices_dir+'badpix.fits',byte(badpix),hdr0;extname='BADPIX'
;  h_temp=headfits(root+badpix_cube+'.fits')
;  h_bp=headfits(root+badpix_cube+'.fits',exten=1)
;  modfits,directory+binned_dir+'badpix_end.fits',0,h_temp
;  modfits,directory+binned_dir+'badpix.fits',0,h_temp
;  modfits,directory+median_dir+'badpix.fits',0,h_temp
;  modfits,directory+slices_dir+'badpix.fits',0,h_temp
;  modfits,directory+binned_dir+'badpix_end.fits',1,h_bp,extname='BADPIX'
;  modfits,directory+binned_dir+'badpix.fits',1,h_bp,extname='BADPIX'
;  modfits,directory+median_dir+'badpix.fits',1,h_bp,extname='BADPIX'
;  modfits,directory+slices_dir+'badpix.fits',1,h_bp,extname='BADPIX'
  
endif


;if badpixel cube has been included, we need to instead red in the 
;masks in each directory and add the new mask 
if badpix_TF eq 'T' and file_test(root+badpix_file) eq 1 then begin
  readcol,root+badpix_file,format='f,f',x_bad,y_bad,comment='#',/SILENT
  
  fits_read,directory+galaxy_ref+'_smoothed_FLUX.fits', spec_in, header1
  side1=sxpar(header1,'NAXIS1')
  side2=sxpar(header1,'NAXIS2')
  badpix=intarr(side1,side2)
  badpix_end=intarr(side1,side2)
  badpix_end[*,*]=1
  
  mkhdr,hdr0,badpix_end
  fits_write,directory+binned_dir+'badpix_end.fits',byte(badpix_end),hdr0
  fits_write,directory+slices_dir+'badpix_end.fits',byte(badpix_end),hdr0
  
  for j=0,n_elements(x_bad)-1,1 do badpix[x_bad,y_bad]=1
  
  
  ;median image
  fits_read,directory+median_dir+'badpix.fits',input
  output=input+badpix
  index = WHERE(output gt 0)
  s = SIZE(output)
  ncol = s(1)
  col = index MOD ncol
  row = index / ncol
  output[col,row]=1
  

  
;  temp=mrdfits(directory+median_dir+'badpix.fits',1,h_bp,/SILENT);
;  modfits,directory+median_dir+'badpix.fits',output,h_bp,exten_no=1
  ;temp=mrdfits(directory+median_dir+'badpix.fits',0,/SILENT);
  modfits,directory+median_dir+'badpix.fits',output;,h_bp,exten_no=1
  
  ;binned images
  bp_files= file_search(directory+binned_dir+'badpix*.fits',COUNT=nfiles)
  for n=0,nfiles-1,1 do begin
    fits_read,bp_files[n],input,h_bp
    output=input+badpix
    index = WHERE(output gt 0)
    s = SIZE(output)
    ncol = s(1)
    col = index MOD ncol
    row = index / ncol
    output[col,row]=1
    
;    temp=mrdfits(bp_files[n],1,h_bp,/SILENT);
    modfits,bp_files[n],output;,h_bp,exten_no=1    
  endfor
  
  
  ;image slices
  bp_files= file_search(directory+slices_dir+'badpix*.fits',COUNT=nfiles)
  for n=0,nfiles-1,1 do begin
    fits_read,bp_files[n],input
    output=input+badpix
    index = WHERE(output gt 0)
    s = SIZE(output)
    ncol = s(1)
    col = index MOD ncol
    row = index / ncol
    output[col,row]=1
    ;temp=mrdfits(bp_files[n],1,h_bp,/SILENT);
    modfits,bp_files[n],output;,h_bp,exten_no=1
  endfor
  
 
  
  
  
endif


delvarx,badpix,badpix_end,input,output,index,temp
end