; this code will identify bad pixels in the MaNGA data cube, 
; which for now are those with a value of 0 outside of the 
; MaNGA field of view
; 
; Modified to create separate masks for the first and last 
; images where all pixels are masked. This will reduce 
; mis-fits due to the different binning techniques for these 
; images (Nov 2014)
; 

pro badpixelmask, root,decomp, galaxy_ref, badpix, binned_dir, median_dir, slices_dir,badpix_file
directory=root+decomp
fits_read,directory+galaxy_ref+'_smoothed_kinematics.fits', spec_in, header1
;fits_read,directory+decomp+median_dir+galaxy_ref+'_median_image.fits', spec_in, header1
side1=sxpar(header1,'NAXIS1')
side2=sxpar(header1,'NAXIS2')
images=sxpar(header1,'NAXIS3')
openw,55,directory+binned_dir+'badpix.pl'
openw,56,directory+slices_dir+'badpix.pl'
openw,57,directory+median_dir+'badpix.pl'
openw,58,directory+binned_dir+'badpix_end.pl'

badpix=intarr(side1,side2)
badpix_end=intarr(side1,side2)

image_slice=spec_in[*,*,fix(images/2)]

;Pick a random image slics, and identify bad pixels 
;as those with a pixel value of 0.0
for x=0,side1-1,1 do begin
  for y=0,side2-1,1 do begin
;    if image_slice[x,y] eq 0 OR y le 6 or y ge 34 then begin
    printf,58,x+1,y+1
    badpix_end[x,y]=1
    if image_slice[x,y] eq 0.0 then begin
        printf,55,x+1,y+1
        printf,56,x+1,y+1
        printf,57,x+1,y+1
        badpix[x,y]=1
        
    endif 
    if x eq 0 or y eq 0 then begin
        printf,55,x+1,y+1
        printf,56,x+1,y+1
        printf,57,x+1,y+1
        badpix[x,y]=1
        
    endif 
  endfor
endfor

if file_test(root+badpix_file) eq 1 then begin
  ;print,'***YESYESYESYESYES***'
 
  readcol,root+badpix_file,format='f,f',x_bad,y_bad,comment='#',/SILENT
  for j=0,n_elements(x_bad)-1,1 do begin
    printf,55,x_bad+1,y_bad+1
    printf,56,x_bad+1,y_bad+1
    printf,57,x_bad+1,y_bad+1
    badpix[x_bad,y_bad]=1
    
  endfor
endif

fits_write,directory+binned_dir+'badpix.fits',byte(badpix)
fits_write,directory+slices_dir+'badpix.fits',byte(badpix)
fits_write,directory+median_dir+'badpix.fits',byte(badpix)
fits_write,directory+binned_dir+'badpix_end.fits',byte(badpix_end)
;printf,57,'21 12'
;printf,57,'21 11'
;printf,57,'20 12'
;printf,57,'20 11'

close,55
close,56
close,57

delvarx,badpix,badpix_end
end