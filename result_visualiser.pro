


pro result_visualiser,root,directory,decomp,galaxy_ref,slices_dir,info,x_centre,y_centre,start_wavelength,end_wavelength,wavelength,Redshift,n_comp,comp3_type,comp4_type,comp4_x,comp4_y,no_slices,MANGA=manga,CALIFA=califa
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
fits_read,root+directory+decomp+'decomposed_data/original.fits',original_datacube,h_orig
fits_read,root+directory+decomp+'decomposed_data/bestfit.fits',bestfit_datacube,h_bestfit
fits_read,root+directory+decomp+'decomposed_data/residuals.fits',residual_datacube,h_resid

fits_read,root+directory+decomp+'decomposed_data/disk.fits',disk_datacube,h_disk
if n_comp ge 1100 then fits_read,root+directory+decomp+'decomposed_data/bulge.fits',bulge_datacube,h_bulge
  if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then fits_read,root+directory+decomp+'decomposed_data/comp3.fits',comp3_datacube,h_comp3
  if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then fits_read,root+directory+decomp+'decomposed_data/comp4.fits',comp4_datacube,h_comp4

npix=sxpar(h_disk,'NAXIS3')
bulge_1D=fltarr(npix)
disk_1D=fltarr(npix)
comp3_1D=fltarr(npix)
comp4_1D=fltarr(npix)
orig_1D=fltarr(npix)
bestfit_1D=fltarr(npix)
resid_1D=fltarr(npix)
sky=fltarr(npix)

bulge_Re=fltarr(npix)
disk_Re=fltarr(npix)
wavelength=fltarr(npix)
for n=0,npix-1,1 do wavelength[n]=10^(sxpar(h_disk,'CRVAL3')+n*sxpar(h_disk,'CD3_3'))

;read in S/N values to identify good pixels for measuring flux
;readcol,root+directory+galaxy_ref+'_S_N_array.txt',format='F,F,F,F',Xpix,Ypix,signal,noise,comment='#'
result=file_test(root+directory+galaxy_ref+'_voronoi_2d_binning_output.txt')

if result eq 1 then begin
  readcol,root+directory+galaxy_ref+'_voronoi_2d_binning_output.txt',format='F,F,F',Xpix,Ypix,BINpix,comment='#'
  
  
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
      if BINpix[n] ne bad_bin then begin
        if n_comp ge 1100  then bulge_1D[*]+=bulge_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then comp3_1D[*]+=comp3_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011  then comp4_1D[*]+=comp4_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        disk_1D+=disk_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        orig_1D+=original_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        bestfit_1D+=bestfit_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        resid_1D+=residual_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        ;print,original_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
        ;print,bulge_datacube[Xpix[n],Ypix[n],*]
      endif
  ENDFOR
endif else begin
  ;***not working yet!!!***
  for aaa=0,npix-1,1 do begin
    disk_1D[aaa]=total(disk_datacube[*,*,aaa])
    orig_1D[aaa]=total(original_datacube[*,*,aaa])
    bestfit_1D[aaa]=total(bestfit_datacube[*,*,aaa])
    resid_1D[aaa]=total(residual_datacube[*,*,aaa])
    if n_comp ge 1100  then bulge_1D[aaa]=total(bulge_datacube[*,*,aaa])
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

temp = file_search(root+directory+decomp+slices_dir+'galfitm_*.feedme',COUNT=nfiles)

for n=0,nfiles-1,1 do begin
  result=file_test(root+directory+decomp+slices_dir+'subcomps_'+string(n,format='(I4.4)')+'.fits')
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
    
    

  if n_comp eq 1000 then res=read_sersic_results_2comp(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1100 then res=read_sersic_results_2comp(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'psf' then res=read_sersic_results_3psf(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1001 and comp4_type eq 'psf' then res=read_sersic_results_3psf(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1001 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  
  else if n_comp eq 1010  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1010  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1110  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1110  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=1) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+directory+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits', nband, bd=0) 


    exptime=sxpar(h_disk,'EXPTIME')
;    print,a,b
;    disk_1D[a:b]=10^((res.mag_galfit_band_d[0:nbands-1]-15)/(-2.5))*exptime
    disk_Re[a:b]=res.Re_galfit_band_d[0:nbands-1]
    sky[a:b]=res.sky_galfit_band[0:nbands-1]*n_elements(Xpix)
    
    if n_comp ge 1100 then begin
;      bulge_1D[a:b]=10^((res.mag_galfit_band_b[0:nbands-1]-15)/(-2.5))*exptime
      bulge_Re[a:b]=res.Re_galfit_band_b[0:nbands-1]
    endif
;    res=read_sersic_results(root+directory+decomp+slices_dir+'models/subcomps_'+string(n,format='(I4.4)')+'.fits', nband, bd=1)
    
;    bulge_1D[n-first_image]=10^((res.mag_galfit_band_b-15)/(-2.5))*exptime
;    disk_1D[n-first_image]=10^((res.mag_galfit_band_d-15)/(-2.5))*exptime
    ;print,bulge_1D[n-944],res.mag_galfit_b
  endelse
endfor    

if n_comp eq 1000 or n_comp eq 1001 then disk_1d=disk_1d+sky $
else if n_comp eq 1100 or n_comp eq 1101 then begin
  disk_1D=disk_1d+(0.5*sky)
  bulge_1D=bulge_1D+(0.5*sky)
endif else if n_comp eq 1111  or n_comp eq 1110 then begin
  disk_1D=disk_1d+(0.333*sky)
  bulge_1D=bulge_1D+(0.333*sky)
  comp3_1D=comp3_1D+(0.333*sky)  
endif else begin
  disk_1D=disk_1d+(0.5*sky)
  comp3_1D=comp3_1D+(0.5*sky)  
endelse



;obtain new bulge spectrum within 3Re
res2=read_sersic_results_2comp(root+directory+decomp+'binned_images/imgblock.fits', nband, bd=1)
readcol,root+directory+decomp+slices_dir+'info.txt',format='X,F',info
PA=res2.PA_GALFIT_BAND_D[0:nbands-1]
no_slices=info[2]
h_temp=headfits(root+directory+decomp+'binned_images/image_'+string(0,format='(I4.4)')+'.fits')
wave1=sxpar(h_temp,'WAVELENG')
h_temp=headfits(root+directory+decomp+'binned_images/image_'+string(no_slices-1,format='(I4.4)')+'.fits')
wave2=sxpar(h_temp,'WAVELENG')
wavelength_Re=5500
Re_b=chebeval(wavelength_Re,res2.RE_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
bulge_1D_small=fltarr(n_elements(bulge_1D))
;stop
measure_circular_radius,indgen(n_elements(Xpix)),Xpix,Ypix,0,0,PA[0], radii
FOR n=0,n_elements(Xpix)-1,1 do begin
    if n_comp ge 1100 and radii[n] le 3*Re_b then bulge_1D_small[*]+=bulge_datacube[Xpix[n]+x_centre,Ypix[n]+y_centre,*]
endfor


;print,bulge_1D
;plot results for comparison
set_plot,'ps'
device,file=root+directory+decomp+'decomposed_data/Spectra_integrated.eps',/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
!P.thick=3
!p.charthick=3
!p.charsize=1.3
!p.multi=0;[0,1,4]
;start_wavelength=4600
;end_wavelength=7000

plot,wavelength,disk_1D,/NODATA,yrange=[-0.1,2.5],$
    xrange=[start_wavelength-100,end_wavelength+100],$
    /xstyle,/ystyle,xthick=3,ythick=3,$;ytickinterval=30,$
   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux',title=galaxy_ref



if n_comp ge 1100 then oplot,wavelength,(bulge_1D/median(orig_1D)),color=fsc_color('blue');/10000;+90
oplot,wavelength,(disk_1D/median(orig_1D)),color=fsc_color('red');/10000;+60
oplot,wavelength,(orig_1D/median(orig_1D));/10000;+30
oplot,wavelength,(bestfit_1D/median(orig_1D)),color=fsc_color('purple');/10000;+30,color=fsc_color('red')
;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=fsc_color('red')
oplot,wavelength,(resid_1D/median(orig_1D)),color=fsc_color('olive');/10000,color=fsc_color('green')

if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength,(comp3_1D/median(orig_1D)),color=fsc_color('skyblue')

if n_comp eq 1000 or n_comp eq 1001 then legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Residuals'],linestyle=[0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('red'),fsc_color('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1100 or n_comp eq 1101 then legend,['Integrated spectrum from datacube','Bulge + Disc','Bulge','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('blue'),fsc_color('red'),fsc_color('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1010 or n_comp eq 1011 then legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('skyblue'),fsc_color('red'),fsc_color('olive')],charsize=1.2,box=0,/left,/top

if n_comp eq 1110 or n_comp eq 1111 then legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('blue'),fsc_color('skyblue'),fsc_color('red'),fsc_color('olive')],charsize=1.2,box=0,/left,/top



if n_comp ge 1100 then print,'bulge',bulge_1D[100]
print,'disk',disk_1D[100]
print,'original',orig_1D[100]
print,'bestfit',bestfit_1D[100]
print,'residuals',resid_1D[100]

;!p.multi=0
device,/close

;update header with correct parameters for 1D spectra.
if n_comp ne 1111 and n_comp ne 1011 and n_comp ne 1110 and n_comp ne 1010 then xx=2 else xx=3
for bd=1,xx,1 do begin
    if n_comp lt 1100 and bd eq 1 then bd=2
    if bd eq 1 then h=h_bulge else if bd eq 2 then h=h_disk else if bd eq 3 then h=h_comp3
    sxdelpar,h,'CRVAL1'
    sxdelpar,h,'CRVAL2'
    sxdelpar,h,'CDELT1'
    sxdelpar,h,'CDELT2'
    sxdelpar,h,'CD1_1'
    sxdelpar,h,'CD2_2'
    
    sxaddpar,h,'CUNIT1','Angstrom'
    temp=sxpar(h,'CRVAL3')
    sxaddpar,h,'CRVAL1',temp
    temp=sxpar(h,'CD3_3')
    sxaddpar,h,'CD1_1',temp
    sxaddpar,h,'CTYPE1','WAVE-LOG'
    sxaddpar,h,'CRPIX1',1
    sxaddpar,h,'CRPIX2',1
    
    sxdelpar,h,'CRVAL3'
    sxdelpar,h,'CD3_3'
    sxdelpar,h,'CTYPE3'
    sxdelpar,h,'CRPIX3'
    sxdelpar,h,'CUNIT3'
    sxdelpar,h,'CUNIT2'
print,bd
    if bd eq 3 then fits_write,root+directory+decomp+'decomposed_data/comp3_1D.fits',comp3_1D,h
    if bd eq 2 then fits_write,root+directory+decomp+'decomposed_data/disk_1D.fits',disk_1D,h
    if bd eq 1 then fits_write,root+directory+decomp+'decomposed_data/bulge_1D.fits',bulge_1D,h
    if bd eq 1 then fits_write,root+directory+decomp+'decomposed_data/bulge_1D_small.fits',bulge_1D_small,h
    if bd eq 2 then fits_write,root+directory+decomp+'decomposed_data/disk_Re.fits',disk_Re,h
    if bd eq 1 then fits_write,root+directory+decomp+'decomposed_data/bulge_Re.fits',bulge_Re,h

endfor



set_plot,'ps'
;device,file='/Users/ejohnsto/Dropbox/papers/Paper4/decomposed_spectra_1D.eps',xsize=19.5,ysize=10,/portrait;,/landscape
device,file=root+directory+decomp+'decomposed_data/Spectra_integrated_2.eps',/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
!P.thick=3
!p.charthick=3
!p.charsize=1.0
!p.multi=0;[0,1,4]
;start_wavelength=4600
end_wavelength=6900



plot,wavelength,disk_1D,/NODATA,yrange=[-0.1,1.8],$
    xrange=[start_wavelength-100,end_wavelength+100],$
    /xstyle,/ystyle,xthick=3,ythick=3,$;ytickinterval=30,$
   ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux';,title=galaxy_ref



if n_comp ge 1100 then oplot,wavelength,(bulge_1D/median(orig_1D)),color=fsc_color('red');/10000;+90
oplot,wavelength,(disk_1D/median(orig_1D)),color=fsc_color('blue');/10000;+60
oplot,wavelength,(orig_1D/median(orig_1D));/10000;+30
oplot,wavelength,(bestfit_1D/median(orig_1D)),color=fsc_color('purple');/10000;+30,color=fsc_color('red')
;oplot,wavelength,((bulge_1D+disk_1D)-median(bulge_1D+disk_1D))/10+10,color=fsc_color('red')
oplot,wavelength,(resid_1D/median(orig_1D)),color=fsc_color('olive');/10000,color=fsc_color('green')

if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then oplot,wavelength,(comp3_1D/median(orig_1D)),color=fsc_color('skyblue')

if n_comp eq 1000 or n_comp eq 1001 then legend,['Integrated spectrum from datacube','Bulge + Disc','Disc','Residuals'],linestyle=[0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('blue'),fsc_color('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1100 or n_comp eq 1101 then legend,['Integrated spectrum from datacube','Bulge + Disc','Bulge','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('red'),fsc_color('blue'),fsc_color('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1010 or n_comp eq 1011 then legend,['Integrated spectrum from datacube','Centre + Disc','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('skyblue'),fsc_color('blue'),fsc_color('olive')],charsize=0.9,box=0,/left,/top

if n_comp eq 1110 or n_comp eq 1111 then legend,['Integrated spectrum from datacube','Centre + Bulge + Disc','Bulge','Centre','Disc','Residuals'],linestyle=[0,0,0,0,0,0],$
  colors=[fsc_color('black'),fsc_color('purple'),fsc_color('red'),fsc_color('skyblue'),fsc_color('blue'),fsc_color('olive')],charsize=0.9,box=0,/left,/top

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


;oplot,[Ha_new,Ha_new],[-5,5],linestyle=1
;oplot,[Hb_new,Hb_new],[-5,5],linestyle=1
;oplot,[Hg_new,Hg_new],[-5,5],linestyle=1
;oplot,[Hd_new,Hd_new],[-5,5],linestyle=1
;oplot,[Mgb_new,Mgb_new],[-5,5],linestyle=1
;oplot,[Fe5335_new,Fe5335_new],[-5,5],linestyle=1
;oplot,[Fe5270_new,Fe5270_new],[-5,5],linestyle=1

a=where(wavelength ge Ha_new)
b=where(wavelength ge Hb_new)
c=where(wavelength ge Hg_new)
d=where(wavelength ge Hd_new)
e=where(wavelength ge Mgb_new)
f=where(wavelength ge Fe5335_new)
g=where(wavelength ge Fe5270_new)
h=where(wavelength ge Na_new)

aa=mean((orig_1D[a[0]-30:a[0]+30]/median(orig_1D)))+0.2
bb=mean((orig_1D[b[0]-30:b[0]+30]/median(orig_1D)))+0.2
cc=mean((orig_1D[c[0]-30:c[0]+30]/median(orig_1D)))+0.2
dd=mean((orig_1D[d[0]-30:d[0]+30]/median(orig_1D)))+0.2
ee=mean((orig_1D[e[0]-30:e[0]+30]/median(orig_1D)))+0.2
ff=mean((orig_1D[f[0]-30:f[0]+30]/median(orig_1D)))+0.2
gg=mean((orig_1D[g[0]-30:g[0]+30]/median(orig_1D)))+0.2
hh=mean((orig_1D[h[0]-30:h[0]+30]/median(orig_1D)))+0.2


xyouts,Ha_new-30,aa,'H'+greek('alpha'),charsize=0.9
xyouts,Hb_new-30,bb,'H'+greek('beta'),charsize=0.9
xyouts,Hg_new-30,cc,'H'+greek('gamma'),charsize=0.9
xyouts,Hd_new-30,dd,'H'+greek('delta'),charsize=0.9
xyouts,Mgb_new-30,ee,'Mg',charsize=0.9
xyouts,Fe5335_new-30,ff,'Fe',charsize=0.9
xyouts,Fe5270_new-30,gg,'Fe',charsize=0.9
xyouts,Na_new-30,hh,'Na D',charsize=0.9

;!p.multi=0
device,/close




end

