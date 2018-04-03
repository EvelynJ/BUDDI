; code to take the residual data cube from BUDDI, 
; and use the ds9 region file created to extract 
; the spectra in a circular aperture around each 
; GC. 
; 

pro GC_residual_extraction,input_file

read_input, input_file, setup

;read in residual data cube
root=setup.root
decomp=setup.decomp
galaxy_ref=setup.galaxy_ref
slices_dir=setup.slices_dir
binned_dir=setup.binned_dir
median_dir=setup.median_dir
decomp_dir=setup.decomp_dir
x_centre=fix(setup.x_centre-1)             ;x position of centre of galaxy, -1 to convert to position in array
y_centre=fix(setup.y_centre-1)             ;y position of centre of galaxy, -1 to convert to position in array
Redshift=setup.Redshift
n_comp=setup.n_comp
comp3_type=setup.comp3_type
comp4_type=setup.comp4_type
no_slices=setup.no_slices

;read in region file with x and y coordinates
fits_read, root+decomp+decomp_dir+'residuals_cube.fits',resid,h
readcol,root+decomp+binned_dir+'detected_objects_psf.reg',format='f,f',x_in,y_in,comment='#',/silent

;extract 1D spectra from residual datacube
radius=2  ;set radius for circular region to extract
;pos_arr=fltarr(2*radius+3,2*radius+3)
radii_arr=fltarr(2*radius+3,2*radius+3)
cen=radius+1
for x=0,2*radius+2,1 do begin
  for y=0,2*radius+2,1 do begin
    ;pos_arr[x,y]=[x-cen,y-cen]
    distance=sqrt((x-cen)^2+(y-cen)^2)
    radii_arr[x,y]=distance
    ;print,x,y,distance
  endfor
endfor

;identify those pixels in array that are ourside of the predefined rasius
index = WHERE(radii_arr gt radius)
s = SIZE(radii_arr)
ncol = s(1)
col = index MOD ncol
row = index / ncol


s1=size(resid)
resid_spectra=fltarr(s1[3],n_elements(x_in))
for n=0,n_elements(x_in)-1,1 do begin
  small=resid[x_in[n]-radius-1:x_in[n]+radius+1,y_in[n]-radius-1:y_in[n]+radius+1,*]
  small[col+x_in,row+y_in,*]=0
  for m=0,s1[3]-1,1 do resid_spectra[m,n]=total(small[*,*,m])
endfor
stop
;plot spectra for quick look

;print spectra for all GCs to a single fits file 



end