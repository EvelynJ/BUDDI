pro bin2d_display_pixels_pos, x, y, counts, pixelSize,pos
  COMPILE_OPT IDL2, HIDDEN

  ; Plots colored pixels with the same pixels size

  PLOT, [MIN(x)-pixelSize,MAX(x)+pixelSize], [MIN(y)-pixelSize,MAX(y)+pixelSize], $
    /NODATA, /XSTYLE, /YSTYLE, XTITLE='pixels', YTITLE='pixels', /ISO, color=cgcolor('black'),$
    xrange=[MIN(x)-pixelSize,MAX(x)+pixelSize],yrange=[MIN(y)-pixelSize,MAX(y)+pixelSize],position=pos
  x1 = [-0.5, -0.5, +0.5, +0.5, -0.5] * pixelSize
  y1 = [+0.5, -0.5, -0.5, +0.5, +0.5] * pixelSize
  color = bytscl(counts)
  FOR j=0, N_ELEMENTS(counts)-1 DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]

END


pro buddi_plot_results, dir,plate_ifu,ONECOMP=onecomp,TWOCOMP=twocomp

if keyword_set(ONECOMP) then decomp_dir='IFU_decomp_1comp/' else decomp_dir='IFU_decomp_2comp/'
fits_read,dir+decomp_dir+'decomposed_data/decomposed_data_ergs/component1_cube.fits',c1_in,h_c1
fits_read,dir+decomp_dir+'decomposed_data/decomposed_data_ergs/original_cube.fits',orig_in
fits_read,dir+decomp_dir+'decomposed_data/decomposed_data_ergs/bestfit_cube.fits',bestfit_in
fits_read,dir+decomp_dir+'decomposed_data/decomposed_data_ergs/residuals_cube.fits',resid_in
fits_read,dir+decomp_dir+'decomposed_data/decomposed_data_ergs/residual_sky_cube.fits',sky_in

fits_read,dir+decomp_dir+'median_image/badpix.fits',bp
bp=abs(bp-1)

s=size(c1_in)

fits_read,dir+decomp_dir+'decomposed_data/decomposed_data_ergs/component1_flux.fits',temp,h_temp
wave0=sxpar(h_temp,'CRVAL1')
step=sxpar(h_temp,'CD1_1')
wave=fltarr(s[3])
for n=0,s[3]-1,1 do wave[n]=wave0+n*step

;create the white light images
c1_im=median(c1_in[*,*,50:-50],dimension=3)
orig_im=median(orig_in[*,*,50:-50],dimension=3)
bestfit_im=median(bestfit_in[*,*,50:-50],dimension=3)
resid_im=median(resid_in[*,*,50:-50],dimension=3)
sky_im=median(sky_in[*,*,50:-50],dimension=3)

if keyword_set(TWOCOMP) then begin
  fits_read,dir+decomp_dir+'decomposed_data/decomposed_data_ergs/component2_cube.fits',c2_in
  c2_im=median(c2_in[*,*,50:-50],dimension=3)
endif




;print the results to a pdf file
set_plot,'ps'
;loadct,0
device,file=dir+decomp_dir+'decomposed_data/decomposed_data_ergs/'+plate_ifu+'_overview.eps',ysize=11,xsize=8,/inches,/color
!P.thick=2
!p.charthick=2
!p.charsize=1
!p.multi=[0,3,6]
lower_limit=0.01
upper_limit=0.95

;range_in_val=max(orig_im)-min(orig_im)
;img_range = BYTSCL(orig_im, MIN=min(orig_im)+(lower_limit*range_in_val), MAX=min(orig_im)+(upper_limit*range_in_val))
;plotimage,orig_im,range=img_range,PIXEL_ASPECT_RATIO=[1,1];,position=[0.8,0.95,0.05,0.2]

;x=findgen(s[1])-(s[1]/2)
;y=findgen(s[2])-(s[2]/2)
;image(orig_im,x,y,axis_style=2)


ncol=s[1]
nrow=s[2]
y=fltarr(ncol)
x=indgen(ncol)
orig_1d=orig_im[*,0]
model_1d=bestfit_im[*,0]
resid_1d=resid_im[*,0]
comp1_1d=c1_im[*,0]
if keyword_set(TWOCOMP) then comp2_1d=c2_im[*,0]

for n=1,nrow-1,1 do begin
  x=[x,indgen(ncol)]
  y=[y,fltarr(ncol)+n]
  orig_1d=[orig_1d,orig_im[*,n]]
  model_1d=[model_1d,bestfit_im[*,n]]
  resid_1d=[resid_1d,resid_im[*,n]]
  comp1_1d=[comp1_1d,c1_im[*,n]]
  if keyword_set(TWOCOMP) then comp2_1d=[comp2_1d,c2_im[*,n]]
endfor
pixelSize=1
bin2d_display_pixels_pos, x,y,orig_1d, pixelSize,[0.05,0.85,0.3,0.99]
bin2d_display_pixels_pos, x,y,model_1d, pixelSize,[0.35,0.85,0.6,0.99]
bin2d_display_pixels_pos, x,y,resid_1d, pixelSize,[0.65,0.85,0.9,0.99]

;plot images of each component
bin2d_display_pixels_pos, x,y,comp1_1d, pixelSize,[0.05,0.65,0.3,0.79]
if keyword_set(TWOCOMP) then bin2d_display_pixels_pos, x,y,comp2_1d, pixelSize,[0.05,0.45,0.3,0.59]



;plot fit parameters
if keyword_set(ONECOMP) then bd=0 else bd=1
nband=12
res_free=read_sersic_results_2comp(dir+decomp_dir+'binned_images/imgblock_free.fits', nband, bd=bd) 
res=read_sersic_results_2comp(dir+decomp_dir+'binned_images/imgblock.fits', nband, bd=bd)
if keyword_set(TWOCOMP) then begin
  decomp_dir2='IFU_decomp_2comp/'
  bd=1
  res2_free=read_sersic_results_2comp(dir+decomp_dir2+'binned_images/imgblock_free.fits', nband, bd=bd)
  res2=read_sersic_results_2comp(dir+decomp_dir2+'binned_images/imgblock.fits', nband, bd=bd)
endif

images = file_search(dir+decomp_dir+'binned_images/image_*.fits',COUNT=nfiles)
wavelength=fltarr(nfiles)
for m=0,nfiles-1,1 do begin
  h=headfits(images[m])
  wavelength[m]=sxpar(h,'WAVELENG')
endfor

symbolsize=0.5
!p.charsize=1.2
x=wavelength

y=res_free.MAG_GALFIT_BAND_D
y2=res.MAG_GALFIT_BAND_D
ymin=max([y,y2])+0.5
ymax=min([y,y2])-0.5
if keyword_set(TWOCOMP) then begin
  y_b=res2_free.MAG_GALFIT_BAND_B
  y2_b=res2.MAG_GALFIT_BAND_B
  ymin=max([y,y_b,y2,y2_b])+0.5
  ymax=min([y,y_b,y2,y2_b])-0.5
endif
plot,x,y,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
  ytitle='m!ITot!N',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,title='Component 1',position=[0.35,0.72,0.625,0.8]
oplot,x,y2,linestyle=0

if keyword_set(TWOCOMP) then begin
  plot,x,y_b,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",YTICKFORMAT="(A1)",$
    /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,title='Component 2',position=[0.625,0.72,0.9,0.8]
  oplot,x,y2_b,linestyle=0
endif




y=res_free.Re_GALFIT_BAND_D
y2=res.Re_GALFIT_BAND_D
ymin=min([y,y2])-2
ymax=max([y,y2])+2
if keyword_set(TWOCOMP) then begin
  y_b=res2_free.Re_GALFIT_BAND_B
  y2_b=res2.Re_GALFIT_BAND_B
  ymin=min([y,y_b,y2,y2_b])-2
  ymax=max([y,y_b,y2,y2_b])+2
endif
plot,x,y,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
  ytitle='Re',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.35,0.64,0.625,0.72]
oplot,x,y2,linestyle=0
if keyword_set(TWOCOMP) then begin
  plot,x,y_b,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
    /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.625,0.64,0.9,0.72]
  oplot,x,y2_b,linestyle=0
endif


y=res_free.N_GALFIT_BAND_D
y2=res.N_GALFIT_BAND_D
ymin=min([y,y2])-0.5
ymax=max([y,y2])+0.5
if keyword_set(TWOCOMP) then begin
  y_b=res2_free.N_GALFIT_BAND_B
  y2_b=res2.N_GALFIT_BAND_B
  ymin=min([y,y_b,y2,y2_b])-0.5
  ymax=max([y,y_b,y2,y2_b])+0.5
endif
plot,x,y,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
  ytitle='n',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.35,0.56,0.625,0.64]
oplot,x,y2,linestyle=0
if keyword_set(TWOCOMP) then begin
  plot,x,y_b,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
    /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.625,0.56,0.9,0.64]
  oplot,x,y2_b,linestyle=0
endif




y=res_free.PA_GALFIT_BAND_D
y2=res.PA_GALFIT_BAND_D
ymin=min([y,y2])-5
ymax=max([y,y2])+5
if keyword_set(TWOCOMP) then begin
  y_b=res2_free.PA_GALFIT_BAND_B
  y2_b=res2.PA_GALFIT_BAND_B
  ymin=min([y,y_b,y2,y2_b])-5
  ymax=max([y,y_b,y2,y2_b])+5
endif
plot,x,y,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
  ytitle='PA',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.35,0.48,0.625,0.56]
oplot,x,y2,linestyle=0
if keyword_set(TWOCOMP) then begin
  plot,x,y_b,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
    /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.625,0.48,0.9,0.56]
  oplot,x,y2_b,linestyle=0
endif



y=res_free.Q_GALFIT_BAND_D
y2=res.Q_GALFIT_BAND_D
ymin=0;min([y,y2])+0.5
ymax=1;max([y,y2])-0.5
if keyword_set(TWOCOMP) then begin
  y_b=res2_free.Q_GALFIT_BAND_B
  y2_b=res2.Q_GALFIT_BAND_B
  ymin=0;min([y,y_b,y2,y2_b])+0.5
  ymax=1;max([y,y_b,y2,y2_b])-0.5
endif
plot,x,y,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50],$
  ytitle='q',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.35,0.40,0.625,0.48],$
  xtitle='Wavelength (AA)'
oplot,x,y2,linestyle=0
if keyword_set(TWOCOMP) then begin
  plot,x,y_b,yrange=[ymin,ymax],xrange=[x[0]-50,x[-1]+50], XTICKFORMAT="(A1)",$
    /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,position=[0.625,0.40,0.9,0.48]
  oplot,x,y2_b,linestyle=0
endif




s=size(orig_in)
bp_cube=fltarr(s[1],s[2],s[3])
for j=0,s[3]-1,1 do bp_cube[*,*,j]=bp

;plot the spectrum for each component
orig_1d_temp=total(orig_in*bp_cube,1) 
orig_1d=total(orig_1d_temp,1)
model_1d_temp=total(bestfit_in*bp_cube,1)
model_1d=total(model_1d_temp,1)
resid_1d_temp=total(resid_in*bp_cube,1)
resid_1d=total(resid_1d_temp,1)
sky_1d_temp=total(sky_in*bp_cube,1)
sky_1d=total(sky_1d_temp,1)  
c1_1d_temp=total(c1_in*bp_cube,1)
c1_1d=total(c1_1d_temp,1)
if keyword_set(TWOCOMP) then begin
  c2_1d_temp=total(c2_in*bp_cube,1)
  c2_1d=total(c2_1d_temp,1)
endif


plot,10^wave,orig_1d/median(orig_1d),yrange=[-0.2,4],$
  xrange=[10^wave[0]-50,10^wave[-1]+50],$
  /xstyle,/ystyle,xthick=3,ythick=3,/NODATA,$;
  xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux',position=[0.05,0.05,0.9,0.35]
;oplot,10^wave,model_1d/median(orig_1d),color=cgcolor('purple')
oplot,10^wave,resid_1d/median(c1_1d),color=cgcolor('olive')
;oplot,10^wave,sky_1d/median(c1_1d),color=cgcolor('grey')
oplot,10^wave,orig_1d/median(orig_1d),color=cgcolor('black')
if keyword_set(ONECOMP) then begin
  oplot,10^wave,resid_1d/median(c1_1d),color=cgcolor('olive')
  oplot,10^wave,sky_1d/median(c1_1d),color=cgcolor('grey')
  oplot,10^wave,c1_1d/median(c1_1d),color=cgcolor('red')
  al_legend,['Integrated spectrum from datacube','Component 1','Sky','Residuals'],linestyle=[0,0,0,0],$
    colors=[cgcolor('black'),cgcolor('red'),cgcolor('grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
endif
if keyword_set(TWOCOMP) then begin
  oplot,10^wave,resid_1d/median(c1_1d+c2_1d),color=cgcolor('olive')
  oplot,10^wave,sky_1d/median(c1_1d+c2_1d),color=cgcolor('grey')
  oplot,10^wave,c1_1d/median(c1_1d+c2_1d),color=cgcolor('red')
  oplot,10^wave,c2_1d/median(c1_1d+c2_1d),color=cgcolor('blue')
  oplot,10^wave,(c1_1d+c2_1d)/median(c1_1d+c2_1d),color=cgcolor('purple') 
;  oplot,10^wave,(c1_1d+c2_1d+sky_1d)/median(c1_1d+c2_1d+sky_1d),color=cgcolor('green')
  al_legend,['Integrated spectrum from datacube','Component 1','Component 2','Component 1 + Component 2','Sky','Residuals'],linestyle=[0,0,0,0,0,0],$
    colors=[cgcolor('black'),cgcolor('red'),cgcolor('blue'),cgcolor('purple'),cgcolor('grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
endif

delvarx,wave,resid_1d,sky_1d,c1_1d,c2_1d,c1_in,c2_in
delvarx,sky_in,resid_in,bestfit_in, orig_in, temp
device,/close

end
