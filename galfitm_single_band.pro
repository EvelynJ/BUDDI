; 
; This code will print out the GalfitM feedme files for single 
; band fits for IFU data. Note that GALFITM is being used fo all 
; fits to maintain consistency in the results and in the codes.
; 
;  MEDIAN keyword- to be used for the median combined image from the datacube
;  SLICES keyword- to be used with individual image slices
;  SINGLE keyword- to be used for single Sersic fit
;  DOUBLE keyword- to be used for double Sersic fit
;
pro Galfitm_single_band,output,median_dir,slices_dir,galaxy_ref,info,x,y,x_centre,$
  y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,estimates_comp4,$
  n_comp,disc_n_poly,bulge_n_poly,MEDIAN=median, SLICES=slices, SINGLE=single,$
  DOUBLE=double
  
first_image=info[0]
final_image=info[1]
no_bins=info[2]
images_per_bin=info[3]
start_wavelength=info[4]
end_wavelength=info[5]
x_centre+=1
y_centre+=1

;identify whether disc n is allowed to vary
if disc_n_poly ne 0 then disc_n_poly=1
if bulge_n_poly ne 0 then bulge_n_poly=1


if keyword_set(single) then s_d='single'
if keyword_set(double) then s_d='double'

if estimates_bulge[0] eq 0 then bulge_type='sersic' else bulge_type='psf'
if estimates_disk[0] eq 0 then disk_type='sersic' else disk_type='psf'
if estimates_comp3[0] eq 0 then comp3_type='sersic' else comp3_type='psf'
if estimates_comp4[0] eq 0 then comp4_type='sersic' else comp4_type='psf'


if keyword_set(median) then begin
  close,60
  openw,60,output+median_dir+'galfitm_'+s_d+'.feedme'
  
  printf,60,'#==============================================================================='
  printf,60,'# IMAGE and GALFIT CONTROL PARAMETERS';+$
  printf, 60, 'A) '+galaxy_ref+'_median_image.fits             # Input data image (FITS file)'
  printf, 60, 'A1) MaNGA             # Band labels (can be omitted if fitting a single band)'
  printf, 60, 'A2) '+string(0.5*(start_wavelength+end_wavelength))+'             # Band wavelengths'
  printf, 60, 'B) imgblock_'+s_d+'.fits       # Output data image block
  printf, 60, 'C) none                # Sigma image name (made from data if blank or "none") 
;  printf, 60, 'D) '+galaxy_ref+'_median_psf.fits           # Input PSF image and (optional) diffusion kernel
  printf, 60, 'D) psf.fits           # Input PSF image and (optional) diffusion kernel
  printf, 60, 'E) 1                   # PSF fine sampling factor relative to data 
  printf, 60, 'F) badpix.fits                # Bad pixel mask (FITS image or ASCII coord list)
;if keyword_set(single) then  printf, 60, 'G) galfitm.constraints_single                # File with parameter constraints (ASCII file)' 
if keyword_set(single) then  printf, 60, 'G) none                # File with parameter constraints (ASCII file)' 
if keyword_set(double) then  printf, 60, 'G) galfitm.constraints                # File with parameter constraints (ASCII file)' 
  printf, 60, 'H) 1    '+string(x,format='(I3.3)')+'   1  '+string(y,format='(I3.3)')+'    # Image region to fit (xmin xmax ymin ymax)'
  printf, 60, 'I) '+string(x,format='(I3.3)')+'    '+string(x,format='(I3.3)')+'          # Size of the convolution box (x y)'
  printf, 60, 'J) '+string(magzpt,format='(F4.1)')+'            # Magnitude photometric zeropoint '
  printf, 60, 'K) '+string(scale[0],format='(F4.2)')+'    '+string(scale[1],format='(F4.2)')+'      # Plate scale (dx dy)    [arcsec per pixel]'
  printf, 60, 'O) regular             # Display type (regular, curses, both)'
  printf, 60, 'P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps'
 
 printf, 60, ' '
 printf, 60, ' '
 ; 
  
  printf, 60, '  # INITIAL FITTING PARAMETERS'
  printf, 60, '#'
  printf, 60, '#   For object type, the allowed functions are: '
  printf, 60, '#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, '
  printf, 60, '#       ferrer, powsersic, sky, and isophote. '
  printf, 60, '#  '
  printf, 60, '#   Hidden parameters will only appear when theyre specified:'
  printf, 60, '#       C0 (diskyness/boxyness), '
  printf, 60, '#       Fn (n=integer, Azimuthal Fourier Modes),'
  printf, 60, '#       R0-R10 (PA rotation, for creating spiral structures).'
  printf, 60, '# '
  printf, 60, '# -----------------------------------------------------------------------------'
  printf, 60, '#   par)    par value(s)    fit toggle(s)    # parameter description '
  printf, 60, '# -----------------------------------------------------------------------------'
  
  printf, 60, ' '
  printf, 60, ' '
  printf, 60, ' '
  printf, 60, ' '

  printf, 60, '# Object number: 1'
  printf, 60, ' 0) sky                    #  object type'
  printf, 60, '  1) 0   1 band #  sky background at center of fitting region [ADUs]'
  printf, 60, '  2) 0      0 band  #  dsky/dx (sky gradient in x)'
  printf, 60, '  3) 0      0 band  #  dsky/dy (sky gradient in y)'
  printf, 60, '  Z) 0                      #  output option (0 = resid., 1 = Dont subtract) '
   
  if keyword_set(SINGLE) then begin 
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
       
      printf, 60, '# Object number: 2'    ;disc
      printf, 60, ' 0) '+disk_type+'                 #  object type'
      printf, 60, ' 1) '+string(x_centre)+'   1 band  #  position x'
      printf, 60, ' 2) '+string(y_centre)+'   1 band  #  position y'
      printf, 60, ' 3) '+string(estimates_disk[1])+'   1 band  #  Integrated magnitude' 
      printf, 60, ' 4) '+string(estimates_disk[2])+'   1 band  #  R_e (half-light radius)   [pix]'
      printf, 60, ' 5) '+string(estimates_disk[3])+'   '+string(disc_n_poly,format='(I1.1)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
      printf, 60, ' 9) '+string(estimates_disk[4])+'   1 band  #  axis ratio (b/a)  '
      printf, 60, '10) '+string(estimates_disk[5])+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
      ;printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'

      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
 
       
      if n_comp eq 1001 or n_comp eq 1011 or n_comp eq 1111 then begin
      printf, 60, '0) '+comp4_type+'                # object type'
      printf, 60, ' 1) '+string(estimates_comp4[1],format='(f6.2)')+'  1  band #  position x, y'
      printf, 60, ' 2) '+string(estimates_comp4[2],format='(f6.2)')+'   1 band  #  position x, y'
      printf, 60, '3) '+string(estimates_comp4[3],format='(f5.2)')+'      1  band     # Integrated magnitude   '  
      if comp4_type eq 'sersic' then begin
        printf, 60, '4) '+string(estimates_comp4[4],format='(f5.2)')+'      1  band     # R_e (half-light radius)   [pix] '  
        printf, 60, '5) '+string(estimates_comp4[5],format='(f5.2)')+'      1  band     # Sersic index n (de Vaucouleurs n=4)  '  
        printf, 60, '9) '+string(estimates_comp4[6],format='(f5.2)')+'      1  band     # axis ratio (b/a)    '  
        printf, 60, '10) '+string(estimates_comp4[7],format='(f5.2)')+'      1  band     # position angle (PA) [deg: Up=0, Left=90]  '  
        printf, 60, 'Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
      endif
    endif  
  endif
  if keyword_set(DOUBLE) then begin 
;      res=read_sersic_results_2comp(output+median_dir+'imgblock_single.fits', 1, bd=0)
      
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
       
  
  
      printf, 60, '# Object number: 2'    ;disc
      printf, 60, ' 0) '+disk_type+'                 #  object type'
      printf, 60, ' 1) '+string(x_centre)+'   1 band  #  position x'
      printf, 60, ' 2) '+string(y_centre)+'   1 band  #  position y'
      printf, 60, ' 3) '+string(estimates_disk[1])+'   1 band  #  Integrated magnitude' 
      printf, 60, ' 4) '+string(estimates_disk[2])+'   1 band  #  R_e (half-light radius)   [pix]'
      printf, 60, ' 5) '+string(estimates_disk[3])+'   '+string(disc_n_poly,format='(I1.1)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
      printf, 60, ' 9) '+string(estimates_disk[4])+'   1 band  #  axis ratio (b/a)  '
      printf, 60, '10) '+string(estimates_disk[5])+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 

      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
       
    if n_comp ge 1100 then begin
      printf, 60, ' # Object number: 3   '    ;bulge
      printf, 60, '   0) '+bulge_type+'                 #  object type'
      printf, 60, ' 1) '+string(x_centre)+'   1 band  #  position x'
      printf, 60, ' 2) '+string(y_centre)+'   1 band  #  position y'
      printf, 60, ' 3) '+string(estimates_bulge[1])+'   1 band  #  Integrated magnitude' 
      printf, 60, ' 4) '+string(estimates_bulge[2])+'   1 band  #  R_e (half-light radius)   [pix]'
      printf, 60, ' 5) '+string(estimates_bulge[3])+'   '+string(bulge_n_poly,format='(I1.1)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
      printf, 60, ' 9) '+string(estimates_bulge[4])+'   1 band  #  axis ratio (b/a)  '
      printf, 60, '10) '+string(estimates_bulge[5])+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
      printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'
  
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
    endif
    if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
      printf, 60, ' # Object number: 3   '    ;bulge
      printf, 60, '   0) '+comp3_type+'                 #  object type'
      printf, 60, ' 1) '+string(x_centre)+'   1 band  #  position x'
      printf, 60, ' 2) '+string(y_centre)+'   1 band  #  position y'
      printf, 60, ' 3) '+string(estimates_comp3[1])+'   1 band  #  Integrated magnitude' 
      if comp3_type eq 'sersic' then begin
        printf, 60, ' 4) '+string(estimates_comp3[2])+'   1 band  #  R_e (half-light radius)   [pix]'
        printf, 60, ' 5) '+string(estimates_comp3[3])+'   1 band  #  Sersic index n (de Vaucouleurs n=4) '
        printf, 60, ' 9) '+string(estimates_comp3[4])+'   1 band  #  axis ratio (b/a)  '
        printf, 60, '10) '+string(estimates_comp3[5])+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
      endif 
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
    endif


    if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
      printf, 60, '0) '+comp4_type+'                # object type'
      printf, 60, ' 1) '+string(estimates_comp4[1],format='(f6.2)')+'  1  band #  position x, y'
      printf, 60, ' 2) '+string(estimates_comp4[2],format='(f6.2)')+'   1 band  #  position x, y'
      printf, 60, '3) '+string(estimates_comp4[3],format='(f5.2)')+'      1  band     # Integrated magnitude   '  
      if comp4_type eq 'sersic' then begin
        printf, 60, '4) '+string(estimates_comp4[4],format='(f5.2)')+'      1  band     # R_e (half-light radius)   [pix] '  
        printf, 60, '5) '+string(estimates_comp4[5],format='(f5.2)')+'      1  band     # Sersic index n (de Vaucouleurs n=4)  '  
        printf, 60, '9) '+string(estimates_comp4[6],format='(f5.2)')+'      1  band     # axis ratio (b/a)    '  
        printf, 60, '10) '+string(estimates_comp4[7],format='(f5.2)')+'      1  band     # position angle (PA) [deg: Up=0, Left=90]  '  
      endif
      printf, 60, 'Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
    endif  
    
  endif
  close,60
endif  

if keyword_set(slices) then begin
  close,60
  openw,60,output+slices_dir+'galfitm.feedme'
  
  close,60
endif

end