; this code will identify bad pixels in the MaNGA data cube, 
; which for now are those with a value of 0 outside of the 
; MaNGA field of view
; 
; Modified to create separate masks for the first and last 
; images where all pixels are masked. This will reduce 
; mis-fits due to the different binning techniques for these 
; images (Nov 2014)
; 

pro badpixelmask, directory, galaxy_ref, badpix, binned_dir, median_dir, slices_dir
fits_read,directory+galaxy_ref+'_smoothed_kinematics.fits', spec_in, header1
side1=sxpar(header1,'NAXIS1')
side2=sxpar(header1,'NAXIS2')
images=sxpar(header1,'NAXIS3')
openw,55,directory+binned_dir+'badpix.pl'
openw,56,directory+slices_dir+'badpix.pl'
openw,57,directory+median_dir+'badpix.pl'

badpix=fltarr(side1,side2)

image_slice=spec_in[*,*,fix(images/2)]
;Pick a random image slics, and identify bad pixels 
;as those with a pixel value of 0.0
for x=0,side1-1,1 do begin
  for y=0,side2-1,1 do begin
;    if image_slice[x,y] eq 0 OR y le 6 or y ge 34 then begin
    if image_slice[x,y] eq 0.0 then begin
        printf,55,x+1,y+1
        printf,56,x+1,y+1
        printf,57,x+1,y+1
        badpix[x,y]=1
    endif 
  endfor
endfor


fits_write,directory+binned_dir+'badpix.fits',badpix
fits_write,directory+slices_dir+'badpix.fits',badpix
fits_write,directory+median_dir+'badpix.fits',badpix
;printf,57,'21 12'
;printf,57,'21 11'
;printf,57,'20 12'
;printf,57,'20 11'

close,55
close,56
close,57
end