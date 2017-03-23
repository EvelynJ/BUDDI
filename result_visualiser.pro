 
 
pro result_visualiser,setup,info,start_wavelength,end_wavelength,wavelength,MANGA=manga,CALIFA=califa
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
  comp4_x=setup.comp4_x
  comp4_y=setup.comp4_y
  no_slices=setup.no_slices








first_image=info[0]
final_image=info[1] 
no_bins=info[2]
images_per_bin=info[3]
start_wavelength=info[4]
end_wavelength=info[5]
total_images=final_image-first_image+1
no_images_final=total_images mod no_slices


wavelength=10^(wavelength)
;stop
; read in datacubes
fits_read,root+decomp+decomp_dir+'original.fits',original_datacube,h_orig
fits_read,root+decomp+decomp_dir+'bestfit.fits',bestfit_datacube,h_bestfit
fits_read,root+decomp+decomp_dir+'residuals.fits',residual_datacube,h_resid

fits_read,root+decomp+median_dir+'badpix.fits',badpix,h_bp

fits_read,root+decomp+decomp_dir+'disk.fits',disk_datacube,h_disk
fits_read,root+decomp+decomp_dir+'residual_sky.fits',residual_sky_datacube,h_sky
if n_comp ge 1100 then fits_read,root+decomp+decomp_dir+'bulge.fits',bulge_datacube,h_bulge
  if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then fits_read,root+decomp+decomp_dir+'comp3.fits',comp3_datacube,h_comp3
  if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then fits_read,root+decomp+decomp_dir+'comp4.fits',comp4_datacube,h_comp4

npix=sxpar(h_disk,'NAXIS3')
bulge_1D=fltarr(npix)
disk_1D=fltarr(npix)
comp3_1D=fltarr(npix)
comp4_1D=fltarr(npix)
orig_1D=fltarr(npix)
bestfit_1D=fltarr(npix)
resid_1D=fltarr(npix)
sky=fltarr(npix)
resid_sky_1D=fltarr(npix)

bulge_Re=fltarr(npix)
disk_Re=fltarr(npix)
wavelength=fltarr(npix)
for n=0,npix-1,1 do wavelength[n]=10^(sxpar(h_disk,'CRVAL3')+n*sxpar(h_disk,'CD3_3'))

;read in S/N values to identify good pixels for measuring flux
;readcol,root+galaxy_ref+'_S_N_array.txt',format='F,F,F,F',Xpix,Ypix,signal,noise,comment='#'
result=file_test(root+galaxy_ref+'_voronoi_2d_binning_output.txt')

if result eq 1 then begin
  readcol,root+galaxy_ref+'_voronoi_2d_binning_output.txt',format='F,F,F',Xpix,Ypix,BINpix,comment='#'
  
  
  ;calculate integrated bulge and disc flux
  ;first identify pixels affected by foreground/background stars so that their flux isn't included
  
  
  
  if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
    x_temp=-long(x_centre-(comp4_x-1))
    y_temp=-long(y_centre-(comp4_y-1))
    temperoo_x=where(Xpix eq x_temp)
    temperoo_y=where(Ypix eq y_temp)
    match, temperoo_x,temperoo_y,x_element,y_element
    bad_bin=BINpix[temperoo_x[x_element]]
  endif else begin
    bad_bin=-99
    BINpix[*]=-999
  endelse
  
  FOR n=0,n_elements(Xpix)-1,1 do begin
      if BINpix[n] ne bad_bin and badpix[Xpix[n]+x_centre,Ypix[n]+y_centre] eq 0 then begin
        if n_comp ge 1100  then bulge_1D[*]+=bulge_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then comp3_1D[*]+=comp3_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011  then comp4_1D[*]+=comp4_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        disk_1D+=disk_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        orig_1D+=original_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        bestfit_1D+=bestfit_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        resid_1D+=residual_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        resid_sky_1D+=residual_sky_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        
        ;print,original_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        ;print,bulge_datacube[Xpix[n],Ypix[n],*]
      endif
  ENDFOR
endif else begin
  ;***not working yet!!!***
  for aaa=0,npix-1,1 do begin
    disk_1D[aaa]=total(disk_datacube[*,*,aaa]*badpix)
    orig_1D[aaa]=total(original_datacube[*,*,aaa]*badpix)
    bestfit_1D[aaa]=total(bestfit_datacube[*,*,aaa]*badpix)
    resid_1D[aaa]=total(residual_datacube[*,*,aaa]*badpix)
    resid_sky_1D[aaa]=total(residual_sky_datacube[*,*,aaa]*badpix)
    if n_comp ge 1100  then bulge_1D[aaa]=total(bulge_datacube[*,*,aaa]*badpix)
  endfor
endelse
;convert magnitude units into flux units where necessary, Magzp is 15 in the feedme files
;bulge_1D=10^((bulge_1D-15)/(-2.5))
;disk_1D=10^((disk_1D-15)/(-2.5))


;create wavelength array
;wavelength=fltarr(npix)
;wavelength0=sxpar(h_disk,'CRVAL3')
;step=sxpar(h_disk,'CD3_3')
;if keyword_set(manga) then begin
;  for m=0,npix-1,1 do wavelength[m]=10^(wavelength0+(step*m))
;  wavelength_log=alog10(wavelength)
;endif else if keyword_set(califa) then begin
;  for m=0,npix-1,1 do wavelength[m]=exp(wavelength0+(step*m))
;  wavelength_log=alog(wavelength)
;endif

temp = file_search(root+decomp+slices_dir+'galfitm_*.feedme',COUNT=nfiles)

for n=0,nfiles-1,1 do begin
  result=file_test(root+decomp+slices_dir+'subcomps_'+string(n,format='(I4.4)')+'.fits')
  if n ne nfiles-1 then nbands=no_slices else nbands=no_images_final   
  a=n*no_slices
  b=a+nbands-1
  if result eq 0 then begin
;    bulge_1D[a:b]=-99;bulge_1D[0]
;    disk_1D[a:b]=-99;disk_1D[0]
    bulge_Re[a:b]=-99
    disk_Re[a:b]=-99
  endif else begin
;    if n ne nfiles-1 then nbands=no_slices else nbands=no_images_final 
    
    
nband=nbands
  if n_comp eq 1000 then res=read_sersic_results_2comp(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1100 then res=read_sersic_results_2comp(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1001 and comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1001 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  
  else if n_comp eq 1010  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1010  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1110  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1110  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) 


    exptime=sxpar(h_disk,'EXPTIME')
;    print,a,b
;    disk_1D[a:b]=10^((res.mag_galfit_band_d[0:nbands-1]-15)/(-2.5))*exptime
    disk_Re[a:b]=res.Re_galfit_band_d[0:nbands-1]
    sky[a:b]=res.sky_galfit_band[0:nbands-1]*n_elements(Xpix)
    
    if n_comp ge 1100 then begin
;      bulge_1D[a:b]=10^((res.mag_galfit_band_b[0:nbands-1]-15)/(-2.5))*exptime
      bulge_Re[a:b]=res.Re_galfit_band_b[0:nbands-1]
    endif
;    res=read_sersic_results(root+decomp+slices_dir+'models/subcomps_'+string(n,format='(I4.4)')+'.fits', nband, bd=1)
    
;    bulge_1D[n-first_image]=10^((res.mag_galfit_band_b-15)/(-2.5))*exptime
;    disk_1D[n-first_image]=10^((res.mag_galfit_band_d-15)/(-2.5))*exptime
    ;print,bulge_1D[n-944],res.mag_galfit_b
  endelse
endfor    

if n_comp eq 1000 or n_comp eq 1001 then begin
  disk_1D_orig=disk_1D
  disk_1d=disk_1d+sky 
endif else if n_comp eq 1100 or n_comp eq 1101 then begin
  disk_1D_orig=disk_1D
  bulge_1D_orig=bulge_1D
  disk_1D=disk_1d+(0.5*sky)
  bulge_1D=bulge_1D+(0.5*sky)
endif else if n_comp eq 1111  or n_comp eq 1110 then begin
  disk_1D_orig=disk_1D
  bulge_1D_orig=bulge_1D
  comp3_1D_orig=comp3_1D
  disk_1D=disk_1d+(0.333*sky)
  bulge_1D=bulge_1D+(0.333*sky)
  comp3_1D=comp3_1D+(0.333*sky)  
endif else begin
  disk_1D_orig=disk_1D
  comp3_1D_orig=comp3_1D
  disk_1D=disk_1d+(0.5*sky)
  comp3_1D=comp3_1D+(0.5*sky)  
endelse



;obtain new bulge spectrum within 3Re

  if n_comp ge 1100 then res2=read_sersic_results_2comp(root+decomp+binned_dir+'imgblock.fits', no_bins, bd=1) $
    else res2=read_sersic_results_2comp(root+decomp+binned_dir+'imgblock.fits', no_bins, bd=0)
  readcol,root+decomp+slices_dir+'info.txt',format='X,F',info
  PA=res2.PA_GALFIT_BAND_D[0:nbands-1]
  no_slices=info[2]
  h_temp=headfits(root+decomp+binned_dir+'image_'+string(0,format='(I4.4)')+'.fits')
  wave1=sxpar(h_temp,'WAVELENG')
  h_temp=headfits(root+decomp+binned_dir+'image_'+string(no_slices-1,format='(I4.4)')+'.fits')
  wave2=sxpar(h_temp,'WAVELENG')
  wavelength_Re=5500
  measure_circular_radius,indgen(n_elements(Xpix)),Xpix,Ypix,0,0,PA[0], radii
if n_comp ge 1100 then begin
  Re_b=chebeval(wavelength_Re,res2.RE_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
  bulge_1D_small=fltarr(n_elements(bulge_1D))
;stop
  FOR n=0,n_elements(Xpix)-1,1 do begin
    if n_comp ge 1100 and radii[n] le 3*Re_b then bulge_1D_small[*]+=bulge_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
  endfor
endif

;print,bulge_1D
;plot results for comparison
result = FILE_TEST(root+decomp+decomp_dir+'masked_flux/', /DIRECTORY)
if result eq 0 then file_mkdir,root+decomp+decomp_dir+'masked_flux/'

set_plot,'ps'
device,file=root+decomp+decomp_dir+'masked_flux/Spectra_integrated.eps',/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
!P.thick=2
!p.charthick=3
!p.charsize=1.3
!p.multi=0;[0,1,4]
;start_wavelength=4600
;end_wavelength=10000;7000

plot,wavelength,disk_1D,/NODATA,yrange=[-0.1,2.5],$
    xrange=[start_wavelength-100,end_wavelength+100],$
    /xstyle,/ystyle,xthick=2,ythick=2,$;ytickinterval=30,$
   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux',title=galaxy_ref



if n_comp ge 1100 then oplot,wavelength,(bulge_1D/median(orig_1D)),color=cgcolor('blue');/10000;+90
oplot,wavelength,(disk_1D/median(orig_1D)),color=cgcolor('red');/10000;+60
oplot,wavelength,(orig_1D/median(orig_1D));/10000;+30
if n_comp ge 1100 then oplot,wavelength,((bulge_1D+disk_1D)/median(orig_1D)),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=cgcolor('red')
oplot,wavelength,(resid_1D/median(orig_1D)),color=cgcolor('olive');/10000,color=cgcolor('green')

if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength,(comp3_1D/median(orig_1D)),color=cgcolor('skyblue')

if n_comp eq 1000 or n_comp eq 1001 then al_legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Residuals'],linestyle=[0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1100 or n_comp eq 1101 then al_legend,['Integrated spectrum from datacube','Bulge + Disc','Bulge','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1010 or n_comp eq 1011 then al_legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('skyblue'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1110 or n_comp eq 1111 then al_legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('skyblue'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top

 




plot,wavelength,disk_1D,/NODATA,yrange=[-0.1,2.5],$
    xrange=[start_wavelength-100,end_wavelength+100],$
    /xstyle,/ystyle,xthick=2,ythick=2,$;ytickinterval=30,$
   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux',title=galaxy_ref



if n_comp ge 1100 then oplot,wavelength,(bulge_1D_orig/median(orig_1D)),color=cgcolor('blue');/10000;+90
oplot,wavelength,(disk_1D_orig/median(orig_1D)),color=cgcolor('red');/10000;+60
oplot,wavelength,(orig_1D/median(orig_1D));/10000;+30
if n_comp ge 1100 then oplot,wavelength,((bulge_1D_orig+disk_1D_orig)/median(orig_1D)),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=cgcolor('red')
oplot,wavelength,(resid_1D/median(orig_1D)),color=cgcolor('olive');/10000,color=cgcolor('green')

if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength,(comp3_1D/median(orig_1D)),color=cgcolor('skyblue')

if n_comp eq 1000 or n_comp eq 1001 then al_legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Residuals'],linestyle=[0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1100 or n_comp eq 1101 then al_legend,['Integrated spectrum from datacube','Bulge + Disc','Bulge','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1010 or n_comp eq 1011 then al_legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('skyblue'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1110 or n_comp eq 1111 then al_legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('skyblue'),cgcolor('red'),cgcolor('olive')],charsize=1.2,box=0,/left,/top



if n_comp ge 1100 then print,'bulge',bulge_1D[100]
print,'disk',disk_1D[100]
print,'original',orig_1D[100]
print,'bestfit',bestfit_1D[100]
print,'residuals',resid_1D[100]

;!p.multi=0
device,/close

;update header with correct parameters for 1D spectra.
h0=headfits(root+decomp+decomp_dir+'original.fits')
tempf=mrdfits(root+decomp+decomp_dir+'original.fits',1,h_flux)
delvarx,tempf
wavelength0=sxpar(h0,'WAVELENG') 
step=sxpar(h0,'CD3_3')                 
if n_comp ne 1111 and n_comp ne 1011 and n_comp ne 1110 and n_comp ne 1010 then xx=2 else xx=3
for bd=1,xx,1 do begin
    if n_comp lt 1100 and bd eq 1 then bd=2
    if bd eq 1 then h=h_bulge else if bd eq 2 then h=h_disk else if bd eq 3 then h=h_comp3
    sxdelpar,h0,'CRVAL1'
    sxdelpar,h0,'CRVAL2'
    sxdelpar,h0,'CDELT1'
    sxdelpar,h0,'CDELT2'
    sxdelpar,h0,'CD1_1'
    sxdelpar,h0,'CD2_2'
    
    sxaddpar,h0,'CUNIT1','Angstrom'
    temp=sxpar(h0,'CRVAL3')
    sxaddpar,h0,'CRVAL1',temp
    temp=sxpar(h0,'CD3_3')
    sxaddpar,h0,'CD1_1',temp
    sxaddpar,h0,'CTYPE1','WAVE-LOG'
    sxaddpar,h0,'CRPIX1',1
    sxaddpar,h0,'CRPIX2',1
    
    sxdelpar,h0,'CRVAL3'
    sxdelpar,h0,'CD3_3'
    sxdelpar,h0,'CTYPE3'
    sxdelpar,h0,'CRPIX3'
    sxdelpar,h0,'CUNIT3'
    sxdelpar,h0,'CUNIT2'

    if bd eq 3 then begin
      fits_write,root+decomp+decomp_dir+'masked_flux/comp3_1D.fits',comp3_1D,extname='FLUX'
      modfits,root+decomp+decomp_dir+'masked_flux/comp3_1D.fits',0,h0
      modfits,root+decomp+decomp_dir+'masked_flux/comp3_1D.fits',1,h_flux,extname='FLUX'
    endif
    if bd eq 2 then begin
      fits_write,root+decomp+decomp_dir+'masked_flux/disk_1D.fits',disk_1D,extname='FLUX'
      modfits,root+decomp+decomp_dir+'masked_flux/disk_1D.fits',0,h0
      modfits,root+decomp+decomp_dir+'masked_flux/disk_1D.fits',1,h_flux,extname='FLUX'
      fits_write,root+decomp+decomp_dir+'masked_flux/disk_Re.fits',disk_Re,extname='FLUX'
      modfits,root+decomp+decomp_dir+'masked_flux/disk_Re.fits',0,h0
      modfits,root+decomp+decomp_dir+'masked_flux/disk_Re.fits',1,h_flux,extname='FLUX'
      fits_write,root+decomp+decomp_dir+'masked_flux/residual_sky_1D.fits',resid_sky_1D,extname='FLUX'
      modfits,root+decomp+decomp_dir+'masked_flux/residual_sky_1D.fits',0,h0
      modfits,root+decomp+decomp_dir+'masked_flux/residual_sky_1D.fits',1,h_flux,extname='FLUX'
    endif
    if bd eq 1 then begin
      fits_write,root+decomp+decomp_dir+'masked_flux/bulge_1D.fits',bulge_1D,extname='FLUX'
      modfits,root+decomp+decomp_dir+'masked_flux/bulge_1D.fits',0,h0
      modfits,root+decomp+decomp_dir+'masked_flux/bulge_1D.fits',1,h_flux,extname='FLUX'
      fits_write,root+decomp+decomp_dir+'masked_flux/bulge_1D_small.fits',bulge_1D_small,extname='FLUX'
      modfits,root+decomp+decomp_dir+'masked_flux/bulge_1D_small.fits',0,h0
      modfits,root+decomp+decomp_dir+'masked_flux/bulge_1D_small.fits',1,h_flux,extname='FLUX'
      fits_write,root+decomp+decomp_dir+'masked_flux/bulge_Re.fits',bulge_Re,extname='FLUX'
      modfits,root+decomp+decomp_dir+'masked_flux/bulge_Re.fits',0,h0
      modfits,root+decomp+decomp_dir+'masked_flux/bulge_Re.fits',1,h_flux,extname='FLUX'
    endif
endfor



set_plot,'ps'
;device,file='/Users/ejohnsto/Dropbox/papers/Paper4/decomposed_spectra_1D.eps',xsize=19.5,ysize=10,/portrait;,/landscape
device,file=root+decomp+decomp_dir+'masked_flux/Spectra_integrated_2.eps',/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
!P.thick=3
!p.charthick=3
!p.charsize=1.0
!p.multi=0;[0,1,4]
;start_wavelength=4600
;end_wavelength=10300;6900



plot,wavelength,disk_1D,/NODATA,yrange=[-0.1,1.8],$
    xrange=[start_wavelength-100,end_wavelength+100],$
    /xstyle,/ystyle,xthick=3,ythick=3,$;ytickinterval=30,$
   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux';,title=galaxy_ref



if n_comp ge 1100 then oplot,wavelength,(bulge_1D/median(orig_1D)),color=cgcolor('red');/10000;+90
oplot,wavelength,(disk_1D/median(orig_1D)),color=cgcolor('blue');/10000;+60
oplot,wavelength,(orig_1D/median(orig_1D));/10000;+30
oplot,wavelength,(bestfit_1D/median(orig_1D)),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
;oplot,wavelength,((bulge_1D+disk_1D)/median(orig_1D)),color=cgcolor('orange')
;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=cgcolor('red')
oplot,wavelength,(resid_1D/median(orig_1D)),color=cgcolor('olive');/10000,color=cgcolor('green')
oplot,wavelength,(resid_sky_1D/median(orig_1D)),color=cgcolor('dark grey');/10000,color=cgcolor('green')

if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength,(comp3_1D/median(orig_1D)),color=cgcolor('skyblue')

if n_comp eq 1000 or n_comp eq 1001 then al_legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Sky','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1100 or n_comp eq 1101 then al_legend,['Integrated spectrum from datacube','Bulge + Disc','Bulge','Disc','Sky','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1010 or n_comp eq 1011 then al_legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Sky','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1110 or n_comp eq 1111 then al_legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Sky','Residuals'],linestyle=[0,0,0,0,0,0,0],$
  colors=[cgcolor('black'),cgcolor('purple'),cgcolor('red'),cgcolor('skyblue'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top

;a=(Redshift+1)*(Redshift+1)
;Velocity=3e5*((a-1)/(a+1))
Ha_new=6563+(6563*Redshift)
Hb_new=4861+(4861*Redshift)
Hg_new=4341+(4341*Redshift)
Hd_new=4102+(4102*Redshift)
Mgb_new=5177+(5177*Redshift)
Fe5335_new=5335+(5335*Redshift)
Fe5270_new=5270+(5270*Redshift)
Na_new= 5895.92+(5895.92*Redshift)
CaII_new=3934+(3934*Redshift)
CaII2_new= 8542.09+(8542.09*Redshift)
Pae_new= 9546+(9546*Redshift)
Paz_new= 9229+(9229*Redshift)



a=where(wavelength ge Ha_new)
b=where(wavelength ge Hb_new)
c=where(wavelength ge Hg_new)
d=where(wavelength ge Hd_new)
e=where(wavelength ge Mgb_new)
f=where(wavelength ge Fe5335_new)
g=where(wavelength ge Fe5270_new)
h=where(wavelength ge Na_new)
i=where(wavelength ge CaII_new)
j=where(wavelength ge CaII2_new)
k=where(wavelength ge Pae_new)
l=where(wavelength ge Paz_new)

if a[0]-30 gt 0 and a[0]+30 lt n_elements(orig_1D)-1 then $
  aa=mean((orig_1D[a[0]-30:a[0]+30]/median(orig_1D)))+0.2 $
  else aa=-1000
if b[0]-30 gt 0 and b[0]+30 lt n_elements(orig_1D)-1 then $
  bb=mean((orig_1D[b[0]-30:b[0]+30]/median(orig_1D)))+0.2 $
  else bb=-1000
if c[0]-30 gt 0 and c[0]+30 lt n_elements(orig_1D)-1 then $
  cc=mean((orig_1D[c[0]-30:c[0]+30]/median(orig_1D)))+0.2 $
  else cc=-1000
if d[0]-30 gt 0 and d[0]+30 lt n_elements(orig_1D)-1 then $
  dd=mean((orig_1D[d[0]-30:d[0]+30]/median(orig_1D)))+0.2 $
  else dd=-1000
if e[0]-30 gt 0 and e[0]+30 lt n_elements(orig_1D)-1 then $
  ee=mean((orig_1D[e[0]-30:e[0]+30]/median(orig_1D)))+0.2 $
  else ee=-1000
if f[0]-30 gt 0 and f[0]+30 lt n_elements(orig_1D)-1 then $
  ff=mean((orig_1D[f[0]-30:f[0]+30]/median(orig_1D)))+0.2 $
  else ff=-1000
if g[0]-30 gt 0 and g[0]+30 lt n_elements(orig_1D)-1 then $
  gg=mean((orig_1D[g[0]-30:g[0]+30]/median(orig_1D)))+0.2 $
  else gg=-1000
if h[0]-30 gt 0 and h[0]+30 lt n_elements(orig_1D)-1 then $
  hh=mean((orig_1D[h[0]-30:h[0]+30]/median(orig_1D)))+0.2 $
  else hh=-1000
if i[0]-30 gt 0 and i[0]+30 lt n_elements(orig_1D)-1 then $
  ii=mean((orig_1D[i[0]-30:i[0]+30]/median(orig_1D)))+0.2 $
  else ii=-1000
if j[0]-30 gt 0 and j[0]+30 lt n_elements(orig_1D)-1 then $
  jj=mean((orig_1D[j[0]-30:j[0]+30]/median(orig_1D)))+0.2 $
  else jj=-1000
if k[0]-30 gt 0 and k[0]+30 lt n_elements(orig_1D)-1 then $
  kk=mean((orig_1D[k[0]-30:k[0]+30]/median(orig_1D)))+0.2 $
  else kk=-1000
if l[0]-30 gt 0 and l[0]+30 lt n_elements(orig_1D)-1 then $
  ll=mean((orig_1D[l[0]-30:l[0]+30]/median(orig_1D)))+0.2 $
  else ll=-1000


xyouts,Ha_new-30,aa,'H'+greek('alpha'),charsize=0.9
xyouts,Hb_new-30,bb,'H'+greek('beta'),charsize=0.9
xyouts,Hg_new-30,cc,'H'+greek('gamma'),charsize=0.9
xyouts,Hd_new-30,dd,'H'+greek('delta'),charsize=0.9
xyouts,Mgb_new-30,ee,'Mg',charsize=0.9
xyouts,Fe5335_new-30,ff,'Fe',charsize=0.9
xyouts,Fe5270_new-30,gg,'Fe',charsize=0.9
xyouts,Na_new-30,hh,'Na D',charsize=0.9
;xyouts,CaII_new-30,0.95,'Ca II',charsize=0.9
xyouts,CaII_new-30,ii,'Ca II',charsize=0.9
xyouts,CaII2_new-30,jj,'Ca II',charsize=0.9

;!p.multi=0
device,/close




end

