; 
; Note- Galfitm -o3 galfitm****.feedme doesn't work to
; produce subcomps files. therefore, we should use the 
; imgblocks to produce Galfit feedme files for each 
; slice, and run those to create subcomps images.
; 
; 
pro subcomps,setup,info,wavelength,CALIFA=califa,MANGA=manga

  root=setup.root
  decomp=setup.decomp
  slices_dir=setup.slices_dir
  decomp_dir=setup.decomp_dir
  no_slices=setup.no_slices
  n_comp=setup.n_comp
  comp3_type=setup.comp3_type
  comp4_type=setup.comp4_type
  galfitm=setup.galfitm

  dir=root+decomp

first_image=info[0]
final_image=info[1]
no_bins=info[2]
images_per_bin=info[3]
start_wavelength=info[4]
end_wavelength=info[5]

galfit_or_galfitm='galfitm'

no_images=no_slices
total_images=final_image-first_image+1
no_loops=long(total_images/no_images)       ;number of feedme files containing no_images
no_images_final=total_images mod no_images   ;number of images in the final feedme file



if galfit_or_galfitm eq 'galfitm' then begin
  
;  result = FILE_TEST(dir+slices_dir+'models/', /DIRECTORY) 
;  if result eq 0 then file_mkdir,dir+slices_dir+'models/'

  ;spawn,'rm '+decomp+slices_dir+'models/*'
  
;  spawn,'cp '+decomp+slices_dir+'imgblock_*.band '+decomp+slices_dir+'models/'
;  spawn,'cp '+decomp+slices_dir+'image* '+decomp+slices_dir+'models/'
;  spawn,'cp '+decomp+slices_dir+'psf.fits '+decomp+slices_dir+'models/'
  
  mag_bulge=fltarr(total_images)
  mag_disk=fltarr(total_images)
  Re_bulge=fltarr(total_images)
  Re_disk=fltarr(total_images)
  n_bulge=fltarr(total_images)
;  wavelength=fltarr(total_images)
  comp3=fltarr(total_images)
  comp4=fltarr(total_images)
  sky=fltarr(total_images)
  h=headfits(dir+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits')
;  wavelength[0]=sxpar(h,'WAVELENG')
;  step=sxpar(h,'CD3_3')
  exptime=sxpar(h,'EXPTIME')
  
;  if keyword_set(manga) then begin
;    for x=1,total_images-1,1 do wavelength[x]=10^(alog10(wavelength[x-1])+(step))
;  endif else if keyword_set(califa) then begin
;    for x=1,total_images-1,1 do wavelength[x]=exp(alog(wavelength[x-1])+(step))
;  endif
  
  close,70
  openw,70,dir+slices_dir+'run_subcomps.sh'
  
  imgblock = file_search(dir+slices_dir+'galfitm_*.feedme',COUNT=nfiles_img)
  ;print,nfiles_img
  
  
  for n=0,nfiles_img-1,1 do begin
    result=file_test(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'.galfit.01.band')
    if result eq 1 then begin
      printf,70,galfitm+' -o3 imgblock_'+string(n,format='(I4.4)')+'.galfit.01.band'
      printf,70,'mv imgblock_'+string(n,format='(I4.4)')+'.fits subcomps_'+string(n,format='(I4.4)')+'.fits'
      if n eq nfiles_img-1 then nbands=no_images_final else nbands=no_slices
      
      nband=nbands
;      if n_comp eq 100 then res=read_sersic_results_2comp(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nbands, bd=0) $
;      else if n_comp eq 110 then res=read_sersic_results_2comp(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nbands, bd=1) $
;      else if n_comp eq 111 and comp3_type eq 'psf' then res=read_sersic_results_3psf(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nbands, bd=1) $
;      else if n_comp eq 111 and comp3_type eq 'sersic' then res=read_sersic_results_3sersic(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nbands, bd=1) $ 
;      else if n_comp eq 101 and comp3_type eq 'psf' then res=read_sersic_results_3psf(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nbands, bd=0) $
;      else if n_comp eq 101 and comp3_type eq 'sersic' then res=read_sersic_results_3sersic(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nbands, bd=0) 
  ;    res=read_sersic_results(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nbands, bd=1)

  if n_comp eq 1000 then res=read_sersic_results_2comp(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1100 then res=read_sersic_results_2comp(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'psf' then res=read_sersic_results_3psf(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1001 and comp4_type eq 'psf' then res=read_sersic_results_3psf(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1001 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  
  else if n_comp eq 1010  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1010  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1110  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1110  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) 


      a=n*no_slices
      b=a+nbands-1
      
      if n_comp ge 1100 then begin
        Re_bulge[a:b]=res.RE_GALFIT_BAND_B[0:nbands-1]
        n_bulge[a:b]=res.n_GALFIT_BAND_B[0:nbands-1]
      endif
      Re_disk[a:b]=res.RE_GALFIT_BAND_D[0:nbands-1]
      ;if n_comp ge 1100 then print,wavelength[a],n,a,b,Re_bulge[a],Re_disk[a],n_bulge[a]
      if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then comp3[a:b]=10^((res.mag_galfit_band_comp3[0:nbands-1]-15)/(-2.5))*exptime
      if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011  then comp4[a:b]=10^((res.mag_galfit_band_comp4[0:nbands-1]-15)/(-2.5))*exptime
      sky[a:b]=res.sky_galfit_band[0:nbands-1]
    endif
  endfor   
  close,70


result = FILE_TEST(dir+'decomposed_data/') 
if result eq 0 then spawn,'mkdir '+dir+'decomposed_data/'
  

wavelength=10^wavelength

set_plot,'ps'
device,file=dir+'decomposed_data/decomp_parameters.eps',/landscape;xoffset=0,yoffset=0,xsize=11,ysize=8,/inches,/color;,/landscape
!P.thick=3
!p.charthick=3
!p.charsize=1.3
!p.multi=[0,1,3] 

plot,wavelength,Re_disk,yrange=[0.8*min(Re_disk),1.2*max(Re_disk)],$
    xrange=[wavelength[0]-100,wavelength[total_images-1]+100],$
    /xstyle,/ystyle,xthick=3,ythick=3,$
    ytitle='Disk Re',$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')'
if n_comp ge 1100 then plot,wavelength,Re_bulge,yrange=[0.8*min(Re_bulge),1.2*max(Re_bulge)],$
    xrange=[wavelength[0]-100,wavelength[total_images-1]+100],$
    /xstyle,/ystyle,xthick=3,ythick=3,$
    ytitle='Bulge Re',$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')'
if n_comp ge 1100 then plot,wavelength,n_bulge,yrange=[0.8*min(n_bulge),1.2*max(n_bulge)],$
    xrange=[wavelength[0]-100,wavelength[total_images-1]+100],$
    /xstyle,/ystyle,xthick=3,ythick=3,$
    ytitle='Bulge n',$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')'

resistant_mean,sky,3,y_mean,y_sigma
plot,wavelength,sky,yrange=[0.8*min(sky),1.2*max(sky)],$
    xrange=[wavelength[0]-100,wavelength[total_images-1]+100],$
    /xstyle,/ystyle,xthick=3,ythick=3,$
    ytitle='Sky (ADU)',$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')'
if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111  then begin
  resistant_mean,comp3,3,y_mean,y_sigma
  plot,wavelength,comp3,yrange=[0.8*min(comp3),1.2*max(comp3)],$
      xrange=[wavelength[0]-100,wavelength[total_images-1]+100],$
      /xstyle,/ystyle,xthick=3,ythick=3,$
      ytitle='comp3 flux',$
      xtitle='Wavelength ('+cgSymbol("angstrom")+')'
endif
if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011  then begin
  resistant_mean,comp4,3,y_mean,y_sigma
  plot,wavelength,comp4,yrange=[0.8*min(comp4),1.2*max(comp4)],$
      xrange=[wavelength[0]-100,wavelength[total_images-1]+100],$
      /xstyle,/ystyle,xthick=3,ythick=3,$
      ytitle='comp4 flux',$
      xtitle='Wavelength ('+cgSymbol("angstrom")+')'
endif

device,/close

if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
  device,file=dir+'decomposed_data/decomp_parameters_comp4.eps',/landscape;xoffset=0,yoffset=0,xsize=11,ysize=8,/inches,/color;,/landscape
  !P.thick=3
  !p.charthick=3
  !p.charsize=1.3
  !p.multi=0 
  resistant_mean,comp4,3,y_mean,y_sigma
  plot,wavelength,comp4,yrange=[0.8*min(comp4),1.2*max(comp4)],$
      xrange=[wavelength[0]-100,wavelength[total_images-1]+100],$
      /xstyle,/ystyle,xthick=3,ythick=3,$
      ytitle='component 4 flux',$
      xtitle='Wavelength ('+cgSymbol("angstrom")+')'
  device,/close
endif

;======================================================================
endif else if galfit_or_galfitm eq 'galfit' then begin



close,70
openw,70,dir+slices_dir+'run_subcomps.sh'

;galfit = file_search(dir+slices_dir+'galfit*.feedme',COUNT=nfiles_gal)  
galfit = file_search(dir+slices_dir+'imgblock_*.fits',COUNT=nfiles_gal)  
  
  
h=headfits(dir+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits')
x_size=sxpar(h,'NAXIS1')
y_size=sxpar(h,'NAXIS2')
step=sxpar(h,'CD3_3')

total_images_used=((nfiles_gal-1)*no_images)+(2*no_images)+1
Re_bulge=fltarr(total_images_used)
Re_disk=fltarr(total_images_used)
n_bulge=fltarr(total_images_used)
wavelength=fltarr(total_images_used)
wavelength[0]=sxpar(h,'WAVELENG')
x=0

for n=0,nfiles_gal-1,1 do begin
;for n=0,1,1 do begin
    ;print,n
  ;identify first and last elements for each imgblock
  if n eq 0 then xx=0 else xx=1
  if n eq nfiles_gal-1 then yy=no_images_final else yy=no_images
  
  
  nband=yy
  res=read_sersic_results(dir+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'.fits', nband, bd=1)
  
  ;print,n,first_image+(n*(no_images-2))+xx,first_image+(n*(no_images-2))+yy
;  print,n,xx,yy
  for m=xx,yy,1 do begin 
;print,m 
    
    image_slice=first_image+(n*(no_images))+m
    print,n,m,res.RE_GALFIT_BAND_B[m]
    ;stop
    wavelength[x]=wavelength[0]+(step*x)
    Re_disk[x]=res.RE_GALFIT_BAND_D[m]
    if n_comp gt 102 then begin
        Re_bulge[x]=res.RE_GALFIT_BAND_B[m]
        n_bulge[x]=res.RE_GALFIT_BAND_B[m]
    endif
    x+=1
    ;print,image_slice
    close,60
    openw,60,dir+slices_dir+'subcomps_'+string(image_slice,format='(I4.4)')+'.feedme'
    
        printf,60,'==============================================================================='
    printf,60,'# IMAGE and GALFIT CONTROL PARAMETERS';+$
    printf, 60, 'A) ../image_'+string(image_slice,format='(I4.4)')+'.fits             # Input data image (FITS file)'
    printf, 60, 'B) imgblock_'+string(image_slice,format='(I4.4)')+'.fits       # Output data image block
    printf, 60, 'C) none                # Sigma image name (made from data if blank or "none") 
    printf, 60, 'D) ../psf.fits           # Input PSF image and (optional) diffusion kernel
    printf, 60, 'E) 1                   # PSF fine sampling factor relative to data 
    printf, 60, 'F) ../badpix.pl                # Bad pixel mask (FITS image or ASCII coord list)
    printf, 60, 'G) ../EXAMPLE.CONSTRAINTS                # File with parameter constraints (ASCII file)' 
    printf, 60, 'H) 1    '+string(x_size,format='(I2.2)')+'   1  '+string(y_size,format='(I2.2)')+'    # Image region to fit (xmin xmax ymin ymax)'
    printf, 60, 'I) 15    15          # Size of the convolution box (x y)'
    printf, 60, 'J) 15              # Magnitude photometric zeropoint '
    printf, 60, 'K) 1  1        # Plate scale (dx dy)    [arcsec per pixel]'
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
    printf, 60, '  1) '+string(res.SKY_GALFIT_BAND[m])+'        0          #  sky background at center of fitting region [ADUs]'
    printf, 60, '  2) 0      0    #  dsky/dx (sky gradient in x)'
    printf, 60, '  3) 0      0    #  dsky/dy (sky gradient in y)'
    printf, 60, '  Z) 0                      #  output option (0 = resid., 1 = Dont subtract) '
     
   
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
     
    printf, 60, '# Object number: 2'    ;disc
    printf, 60, ' 0) sersic                 #  object type'
    printf, 60, ' 1) '+string(res.X_GALFIT_BAND_D[m],format='(F5.2)')+' '+string(res.Y_GALFIT_BAND_D[0],format='(F5.2)')+'  0 0  #  position x, y'
    printf, 60, ' 3) '+string(res.MAG_GALFIT_BAND_D[m],format='(F5.2)')+'        0   #  Integrated magnitude' 
    printf, 60, ' 4) '+string(res.RE_GALFIT_BAND_D[m],format='(F6.2)')+'   0    #  R_e (half-light radius)   [pix]'
    printf, 60, ' 5) '+string(res.N_GALFIT_BAND_D[m],format='(F4.2)')+'             0    #  Sersic index n (de Vaucouleurs n=4) '
    printf, 60, ' 9) '+string(res.Q_GALFIT_BAND_D[m],format='(F4.2)')+'        0    #  axis ratio (b/a)  '
    printf, 60, '10) '+string(res.PA_GALFIT_BAND_D[m],format='(F6.2)')+'   0    #  position angle (PA) [deg: Up=0, Left=90]'
    printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
    ;printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'
   
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
    if n_comp gt 102 then begin
      printf, 60, ' # Object number: 3   '    ;bulge
      printf, 60, '   0) sersic                 #  object type'
      printf, 60, ' 1) '+string(res.X_GALFIT_BAND_B[m],format='(F5.2)')+' '+string(res.Y_GALFIT_BAND_B[0],format='(F5.2)')+'  0 0  #  position x, y'
      printf, 60, ' 3) '+string(res.MAG_GALFIT_BAND_B[m],format='(F5.2)')+'        0   #  Integrated magnitude' 
      printf, 60, ' 4) '+string(res.RE_GALFIT_BAND_B[m],format='(F6.2)')+'   0    #  R_e (half-light radius)   [pix]'
      printf, 60, ' 5) '+string(res.N_GALFIT_BAND_B[m],format='(F4.2)')+'             0    #  Sersic index n (de Vaucouleurs n=4) '
      printf, 60, ' 9) '+string(res.Q_GALFIT_BAND_B[m],format='(F4.2)')+'        0    #  axis ratio (b/a)  '
      printf, 60, '10) '+string(res.PA_GALFIT_BAND_B[m],format='(F6.2)')+'   0    #  position angle (PA) [deg: Up=0, Left=90]'
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
    ;  printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'
   
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
    endif
;    printf, 60, '0) psf                # object type'
;    printf, 60, ' 1) '+string(res.X_GALFIT_BAND_B[m],format='(F5.2)')+' '+string(res.Y_GALFIT_BAND_B[0],format='(F5.2)')+'   0 0  #  position x, y'
;    printf, 60, '3) '+string(res.mag_galfit_band_psf[m],format='(F5.2)')+'      0      # total magnitude   '  
;    printf, 60, 'Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
; 
 
 
 
    printf, 60, '================================================================================'
 
 
 
    
    ;printf,70,'/Users/ejohnsto/Galfit_files/galfit galfit_'+string(image_slice,format='(I4.4)')+'.feedme'
    
    ;printf,70,'/Users/ejohnsto/Galfit_files/galfit -o3 galfit_'+string(image_slice,format='(I4.4)')+'.feedme'
    printf,70,'/Users/ejohnsto/Galfit_files/galfit -o3  subcomps_'+string(image_slice,format='(I4.4)')+'.feedme'
    printf,70,'mv subcomps.fits subcomps_'+string(image_slice,format='(I4.4)')+'.fits'
    ;printf,70,galfitm+'  subcomps_'+string(image_slice,format='(I4.4)')+'.feedme'
    ;printf,70,'mv imgblock_'+string(image_slice,format='(I4.4)')+'.fits subcomps_'+string(image_slice,format='(I4.4)')+'.fits'
    
    close,60
  endfor
  
endfor

close,70

;set_plot,'ps'
;device,file=dir+'decomposed_data/decomp_parameters.eps',xoffset=0,yoffset=0,xsize=11,ysize=8,/inches,/color;,/landscape
;!P.thick=3
;!p.charthick=3
;!p.charsize=1.3
;!p.multi=[0,1,3]
;plot,wavelength,Re_disk,yrange=[0,5],$
;    xrange=[wavelength[0]-100,wavelength[total_images_used-1]+100],$
;    /xstyle,/ystyle,xthick=3,ythick=3,$
;    ytitle='Disk Re',$
;    xtitle='Wavelength ('+cgSymbol("angstrom")+')'
;plot,wavelength,Re_bulge,yrange=[3,8],$
;    xrange=[wavelength[0]-100,wavelength[total_images_used-1]+100],$
;    /xstyle,/ystyle,xthick=3,ythick=3,$
;    ytitle='Bulge Re',$
;    xtitle='Wavelength ('+cgSymbol("angstrom")+')'
;plot,wavelength,n_bulge,yrange=[0,5],$
;    xrange=[wavelength[0]-100,wavelength[total_images_used-1]+100],$
;    /xstyle,/ystyle,xthick=3,ythick=3,$
;    ytitle='Bulge n',$
;    xtitle='Wavelength ('+cgSymbol("angstrom")+')'
;
;device,/close
  
endif else begin
  
  error_message,'subcomps.pro ==> Please choose Galfit or Galfitm'
  
endelse

end