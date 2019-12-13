;code to print out images of the decomposed components from BUDDI
;
;===============================================================
pro bin2d_display_pixels, x, y, counts, pixelSize
  COMPILE_OPT IDL2, HIDDEN

  ; Plots colored pixels with the same pixels size

  PLOT, [MIN(x)-pixelSize,MAX(x)+pixelSize], [MIN(y)-pixelSize,MAX(y)+pixelSize], $
    /NODATA, /XSTYLE, /YSTYLE, XTITLE='pixels', YTITLE='pixels', /ISO, color=cgcolor('black'),$
    xrange=[MIN(x)-pixelSize,MAX(x)+pixelSize],yrange=[MIN(y)-pixelSize,MAX(y)+pixelSize]
  x1 = [-0.5, -0.5, +0.5, +0.5, -0.5] * pixelSize
  y1 = [+0.5, -0.5, -0.5, +0.5, +0.5] * pixelSize
  color = bytscl(counts)
  FOR j=0, N_ELEMENTS(counts)-1 DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]

END


;##############################################################

pro model_image
;dir1='/raid/ejohnston/double_nucleated/0338/BUDDI/IFU_decomp/'
;dir1='/raid/ejohnston/FCC306/BUDDI/IFU_decomp/'
;dir1='/data2/ejohnston/Fornax_new_sample/FCC182/BUDDI/IFU_decomp/'
dir1='/raid/ejohnston/MUSE_S0s/2MIG1814/IFU_decomp/'
dir=dir1+'decomposed_data5/'
binned_dir='binned_images5/'

ncomp=1111
bulge_type='sersic'
comp3_type='sersic'
comp4_type='psf'

;read in the decomposed datacubes for each compomnent
fits_read,dir+'original_cube.fits',original,h_orig
fits_read,dir+'bestfit_cube.fits',model,h_mod
fits_read,dir+'residuals_cube.fits',residual,h_res
fits_read,dir+'component1_cube.fits',comp1,h_c1

if ncomp ge 1100 then fits_read,dir+'component2_cube.fits',comp2,h_c2
if ncomp ge 1110 then fits_read,dir+'component3_cube.fits',comp3,h_c3
if ncomp ge 1111 then fits_read,dir+'component4_cube.fits',comp4,h_c4
;fits_read,dir+'component5_cube.fits',comp5,h_c5


;create median images for each component 
im_orig=median(original[*,*,10:-10],dimension=3)
im_model=median(model[*,*,10:-10],dimension=3)
im_resid=median(residual[*,*,10:-10],dimension=3)
im_comp1=median(comp1[*,*,10:-10],dimension=3)
if ncomp ge 1100 then im_comp2=median(comp2[*,*,10:-10],dimension=3)
if ncomp ge 1110 then im_comp3=median(comp3[*,*,10:-10],dimension=3)
if ncomp ge 1111 then im_comp4=median(comp4[*,*,10:-10],dimension=3)
;im_comp5=median(comp5[*,*,10:-10],dimension=3)

spawn,'mkdir '+dir+'images/'
fits_write,dir+'images/orig.fits',im_orig,h_orig
fits_write,dir+'images/model.fits',im_model,h_orig
fits_write,dir+'images/resid.fits',im_resid,h_orig
fits_write,dir+'images/comp1.fits',im_comp1,h_orig
if ncomp ge 1100 then fits_write,dir+'images/comp2.fits',im_comp2,h_orig
if ncomp ge 1110 then fits_write,dir+'images/comp3.fits',im_comp3,h_orig
if ncomp ge 1111 then fits_write,dir+'images/comp4.fits',im_comp4,h_orig
;fits_write,dir+'images/comp5.fits',im_comp5,h_orig


;print out the binned fit parameters to a file for each component
if ncomp eq 1000 then res=read_sersic_results_2comp(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
else if ncomp eq 1100  and  bulge_type eq 'sersic' then res=read_sersic_results_2comp(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
else if ncomp eq 1100  and  bulge_type eq 'psf' then res=read_sersic_results_2comp_p(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
else if ncomp eq 1101 and  comp4_type eq 'psf' then res=read_sersic_results_3psf(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
else if ncomp eq 1101 and  comp4_type eq 'sersic' then res=read_sersic_results_3sersic(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
else if ncomp eq 1001 and  comp4_type eq 'psf' then res=read_sersic_results_3psf(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
else if ncomp eq 1001 and  comp4_type eq 'sersic' then res=read_sersic_results_3sersic(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
;else if ncomp eq 1110  and  comp3_type eq 'psf' then res=read_sersic_results_2comp_p(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
;else if ncomp eq 1110  and  comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'psf' and  comp3_type eq 'psf' then res=read_sersic_results_3psf_p(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'psf' and  comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'sersic' and  comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'sersic' and  comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
else if ncomp eq 1110 and bulge_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
else if ncomp eq 1110 and bulge_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_p(dir1+binned_dir+'imgblock.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1110 and bulge_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(dir1+binned_dir+'imgblock.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1110 and bulge_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(dir1+binned_dir+'imgblock.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1111 and bulge_type eq 'psf' and comp3_type eq 'psf' and comp4_type eq 'psf' then res=read_sersic_results_4psf_p_p(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'psf' and comp4_type eq 'psf' then res=read_sersic_results_4psf_p(dir1+binned_dir+'imgblock.fits', nband, bd=1) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'psf' and comp4_type eq 'sersic' then res=read_sersic_results_3psf_s(dir1+binned_dir+'imgblock.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'sersic' and comp4_type eq 'psf' then res=read_sersic_results_4sersic_p(dir1+binned_dir+'imgblock.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'sersic' and comp4_type eq 'sersic' then res=read_sersic_results_3sersic_s(dir1+binned_dir+'imgblock.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1011 and  comp4_type eq 'psf' and  comp3_type eq 'psf' then res=read_sersic_results_3psf_p(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
else if ncomp eq 1011 and  comp4_type eq 'psf' and  comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
else if ncomp eq 1011 and  comp4_type eq 'sersic' and  comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
else if ncomp eq 1011 and  comp4_type eq 'sersic' and  comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(dir1+binned_dir+'imgblock.fits', nband, bd=0)

if ncomp eq 1000 then res=read_sersic_results_2comp(dir1+binned_dir+'imgblock.fits', nband, bd=0) $
else if ncomp eq 1100  and  bulge_type eq 'sersic' then res_free=read_sersic_results_2comp(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
else if ncomp eq 1100  and  bulge_type eq 'psf' then res_free=read_sersic_results_2comp_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=0) $
else if ncomp eq 1101 and  comp4_type eq 'psf' then res_free=read_sersic_results_3psf(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
else if ncomp eq 1101 and  comp4_type eq 'sersic' then res_free=read_sersic_results_3sersic(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
else if ncomp eq 1001 and  comp4_type eq 'psf' then res_free=read_sersic_results_3psf(dir1+binned_dir+'imgblock_free.fits', nband, bd=0) $
else if ncomp eq 1001 and  comp4_type eq 'sersic' then res_free=read_sersic_results_3sersic(dir1+binned_dir+'imgblock_free.fits', nband, bd=0) $
;else if ncomp eq 1110  and  comp3_type eq 'psf' then res_free=read_sersic_results_2comp_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
;else if ncomp eq 1110  and  comp3_type eq 'sersic' then res_free=read_sersic_results_2comp_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'psf' and  comp3_type eq 'psf' then res_free=read_sersic_results_3psf_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'psf' and  comp3_type eq 'sersic' then res_free=read_sersic_results_3psf_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'sersic' and  comp3_type eq 'psf' then res_free=read_sersic_results_3sersic_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
;else if ncomp eq 1111 and  comp4_type eq 'sersic' and  comp3_type eq 'sersic' then res_free=read_sersic_results_3sersic_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
else if ncomp eq 1110 and bulge_type eq 'psf' and comp3_type eq 'psf' then res_free=read_sersic_results_3psf_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
else if ncomp eq 1110 and bulge_type eq 'psf' and comp3_type eq 'sersic' then res_free=read_sersic_results_3sersic_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1110 and bulge_type eq 'sersic' and comp3_type eq 'psf' then res_free=read_sersic_results_2comp_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1110 and bulge_type eq 'sersic' and comp3_type eq 'sersic' then res_free=read_sersic_results_2comp_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1111 and bulge_type eq 'psf' and comp3_type eq 'psf' and comp4_type eq 'psf' then res_free=read_sersic_results_4psf_p_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'psf' and comp4_type eq 'psf' then res_free=read_sersic_results_4psf_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'psf' and comp4_type eq 'sersic' then res_free=read_sersic_results_3psf_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'sersic' and comp4_type eq 'psf' then res_free=read_sersic_results_4sersic_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1111 and bulge_type eq 'sersic' and comp3_type eq 'sersic' and comp4_type eq 'sersic' then res_free=read_sersic_results_3sersic_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=1,boxy=boxy_yn) $
else if ncomp eq 1011 and  comp4_type eq 'psf' and  comp3_type eq 'psf' then res_free=read_sersic_results_3psf_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=0) $
else if ncomp eq 1011 and  comp4_type eq 'psf' and  comp3_type eq 'sersic' then res_free=read_sersic_results_3psf_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=0) $
else if ncomp eq 1011 and  comp4_type eq 'sersic' and  comp3_type eq 'psf' then res_free=read_sersic_results_3sersic_p(dir1+binned_dir+'imgblock_free.fits', nband, bd=0) $
else if ncomp eq 1011 and  comp4_type eq 'sersic' and  comp3_type eq 'sersic' then res_free=read_sersic_results_3sersic_s(dir1+binned_dir+'imgblock_free.fits', nband, bd=0)

wave=fltarr(10)
for j=1,10,1 do begin
  h=headfits(dir1+binned_dir+'image_'+string(j,format='(i4.4)')+'.fits')
  wave[j-1]=sxpar(h,'WAVELENG')
endfor

openw,01,dir+'images/binned_free_fits.txt'
openw,02,dir+'images/binned_fits.txt'
printf,01,'#component  wavelength  mag  Re  n  q  PA'
printf,02,'#component  wavelength  mag  Re  n  q  PA'
for j=1,10,1 do begin
  printf,01,'1',wave[j-1],res_free.MAG_GALFIT_BAND_D[j],res_free.RE_GALFIT_BAND_D[j],res_free.N_GALFIT_BAND_D[j],res_free.Q_GALFIT_BAND_D[j],res_free.PA_GALFIT_BAND_D[j],format='(a3,6f10.2)'
  printf,02,'1',wave[j-1],res.MAG_GALFIT_BAND_D[j],res.RE_GALFIT_BAND_D[j],res.N_GALFIT_BAND_D[j],res.Q_GALFIT_BAND_D[j],res.PA_GALFIT_BAND_D[j],format='(a3,6f10.2)'
endfor
for j=1,10,1 do begin
  if ncomp ge 1100 and bulge_type eq 'sersic' then printf,01,'2',wave[j-1],res_free.MAG_GALFIT_BAND_B[j],res_free.RE_GALFIT_BAND_B[j],res_free.N_GALFIT_BAND_B[j],res_free.Q_GALFIT_BAND_B[j],res_free.PA_GALFIT_BAND_B[j],format='(a3,6f10.2)'
  if ncomp ge 1100 and bulge_type eq 'sersic' then printf,02,'2',wave[j-1],res.MAG_GALFIT_BAND_B[j],res.RE_GALFIT_BAND_B[j],res.N_GALFIT_BAND_B[j],res.Q_GALFIT_BAND_B[j],res.PA_GALFIT_BAND_B[j],format='(a3,6f10.2)'
  if ncomp ge 1100 and bulge_type eq 'psf' then printf,01,'2',wave[j-1],res_free.MAG_GALFIT_BAND_B[j],0,0,0,0,format='(a3,6f10.2)'
  if ncomp ge 1100 and bulge_type eq 'psf' then printf,02,'2',wave[j-1],res.MAG_GALFIT_BAND_B[j],0,0,0,0,format='(a3,6f10.2)'
endfor
for j=1,10,1 do begin
  if ncomp ge 1110 and comp3_type eq 'sersic' then printf,01,'3',wave[j-1],res_free.MAG_GALFIT_BAND_COMP3[j],res_free.RE_GALFIT_BAND_COMP3[j],res_free.N_GALFIT_BAND_COMP3[j],res_free.Q_GALFIT_BAND_COMP3[j],res_free.PA_GALFIT_BAND_COMP3[j],format='(a3,6f10.2)'
  if ncomp ge 1110 and comp3_type eq 'sersic' then printf,02,'3',wave[j-1],res.MAG_GALFIT_BAND_COMP3[j],res.RE_GALFIT_BAND_COMP3[j],res.N_GALFIT_BAND_COMP3[j],res.Q_GALFIT_BAND_COMP3[j],res.PA_GALFIT_BAND_COMP3[j],format='(a3,6f10.2)'
  if ncomp ge 1110 and comp3_type eq 'psf' then printf,01,'3',wave[j-1],res_free.MAG_GALFIT_BAND_COMP3[j],0,0,0,0,format='(a3,6f10.2)'
  if ncomp ge 1110 and comp3_type eq 'psf' then printf,02,'3',wave[j-1],res.MAG_GALFIT_BAND_COMP3[j],0,0,0,0,format='(a3,6f10.2)'
endfor
for j=1,10,1 do begin
  if ncomp eq 1111 and comp4_type eq 'sersic' then printf,01,'4',wave[j-1],res_free.MAG_GALFIT_BAND_COMP4[j],res_free.RE_GALFIT_BAND_COMP4[j],res_free.N_GALFIT_BAND_COMP4[j],res_free.Q_GALFIT_BAND_COMP4[j],res_free.PA_GALFIT_BAND_COMP4[j],format='(a3,6f10.2)'
  if ncomp eq 1111 and comp4_type eq 'sersic' then printf,02,'4',wave[j-1],res.MAG_GALFIT_BAND_COMP4[j],res.RE_GALFIT_BAND_COMP4[j],res.N_GALFIT_BAND_COMP4[j],res.Q_GALFIT_BAND_COMP4[j],res.PA_GALFIT_BAND_COMP4[j],format='(a3,6f10.2)'
  if ncomp eq 1111 and comp4_type eq 'psf' then printf,01,'4',wave[j-1],res_free.MAG_GALFIT_BAND_COMP4[j],0,0,0,0,format='(a3,6f10.2)'
  if ncomp eq 1111 and comp4_type eq 'psf' then printf,02,'4',wave[j-1],res.MAG_GALFIT_BAND_COMP4[j],0,0,0,0,format='(a3,6f10.2)'
endfor

close,01
close,02


;obtain spectra for each component within 1Re



;print the results to a pdf file
set_plot,'ps'
device,file=dir+'model_images.eps',ysize=11,xsize=8,/inches,/color;,/landscape
!P.thick=2
!p.charthick=3
!p.charsize=1.3
!p.multi=[0,3,4]

;;im1=image(im_orig)
;RESULT = BYTSCL(im_orig, MIN=0, MAX=100)
;plotimage,RESULT

s=size(im_orig)
ncol=s[1]
nrow=s[2]
y=fltarr(ncol)
x=indgen(ncol)
orig_1d=im_orig[*,0]
model_1d=im_model[*,0]
resid_1d=im_resid[*,0]
comp1_1d=im_comp1[*,0]
if ncomp ge 1100 then comp2_1d=im_comp2[*,0]
if ncomp ge 1110 then comp3_1d=im_comp3[*,0]
if ncomp ge 1111 then comp4_1d=im_comp4[*,0]
;comp5_1d=im_comp5[*,0]

for n=1,nrow-1,1 do begin
  x=[x,indgen(ncol)]
  y=[y,fltarr(ncol)+n]
  orig_1d=[orig_1d,im_orig[*,n]]
  model_1d=[model_1d,im_model[*,n]]
  resid_1d=[resid_1d,im_resid[*,n]]
  comp1_1d=[comp1_1d,im_comp1[*,n]]
  if ncomp ge 1100 then comp2_1d=[comp2_1d,im_comp2[*,n]]
  if ncomp ge 1110 then comp3_1d=[comp3_1d,im_comp3[*,n]]
  if ncomp ge 1111 then comp4_1d=[comp4_1d,im_comp4[*,n]]
  ;comp5_1d=[comp5_1d,im_comp5[*,n]]
endfor
pixelSize=1
;TV, im_orig

lower_limit=0.01
upper_limit=0.95

s=size(orig_1d)
arr=(orig_1d)     ;convert to a single column array

no_NAN=where(finite(arr))  ;identify NAN values
arr_no_NAN=arr[no_NAN]     ;remove NAN values
new_order=sort(arr_no_NAN)   ;resort the real values into ascending order
arr_sorted=arr_no_NAN[new_order]
num=n_elements(arr_sorted)
lim1=arr_sorted[lower_limit*num]
lim2=arr_sorted[upper_limit*num]

orig_1d=bytscl(orig_1d,min=lim1,max=lim2)
model_1d=bytscl(model_1d,min=lim1,max=lim2)
comp1_1d=bytscl(comp1_1d,min=lim1,max=lim2)
if ncomp ge 1100 then comp2_1d=bytscl(comp2_1d,min=lim1,max=lim2)
if ncomp ge 1110 then comp3_1d=bytscl(comp3_1d,min=lim1,max=lim2)
if ncomp ge 1111 then comp4_1d=bytscl(comp4_1d,min=lim1,max=lim2)
;comp5_1d=bytscl(comp5_1d,min=lim1,max=lim2)



lower_limit=0.01
upper_limit=0.99
s=size(resid_1d)
arr=(resid_1d)     ;convert to a single column array

no_NAN=where(finite(arr))  ;identify NAN values
arr_no_NAN=arr[no_NAN]     ;remove NAN values
new_order=sort(arr_no_NAN)   ;resort the real values into ascending order
arr_sorted=arr_no_NAN[new_order]
num=n_elements(arr_sorted)
lim1=arr_sorted[lower_limit*num]
lim2=arr_sorted[upper_limit*num]


resid_1d=bytscl(resid_1d,min=lim1,max=lim2)
bin2d_display_pixels, x,y,orig_1d, pixelSize
bin2d_display_pixels, x,y,model_1d, pixelSize
bin2d_display_pixels, x,y,resid_1d, pixelSize
bin2d_display_pixels, x,y,comp1_1d, pixelSize
if ncomp ge 1100 then bin2d_display_pixels, x,y,comp2_1d, pixelSize
if ncomp ge 1110 then bin2d_display_pixels, x,y,comp3_1d, pixelSize
if ncomp ge 1111 then bin2d_display_pixels, x,y,comp4_1d, pixelSize
;bin2d_display_pixels, x,y,comp5_1d, pixelSize

device,/close
stop
end