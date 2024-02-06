; quick code to read in decomposed_spectra from BUDDI output
; and convert back into units of erg/s/cm2/AA. Code then
; repeats the plots and spectra created in result_visualiser 
; to the directory decomposed_data_ergs
; 

pro buddi_manga_final_flux,input_file


read_input, input_file, setup



;*** Set up directory structure for all output files
root=setup.root
decomp=setup.decomp
decomp_dir=setup.decomp_dir
psf_file=setup.psf_file
stellib_dir=setup.stellib_dir
n_comp=setup.n_comp
file=setup.file

cubes=file_search(root+decomp+decomp_dir+'*_cube.fits',count=nfiles_cube)
spec=file_search(root+decomp+decomp_dir+'*_flux.fits',count=nfiles_1D)

if file_test(root+decomp+decomp_dir+'decomposed_data_ergs/',/DIRECTORY) eq 0 then $
  file_mkdir,root+decomp+decomp_dir+'decomposed_data_ergs/'

cut=10   ;remove the last n elements in the spectra, usually due to poor fits

for j=0,nfiles_1D-1,1 do begin
  ;read in datacube
  fits_read,spec[j],input,h
  s=size(input)
  wave=fltarr(s[1])
  wave_lin=fltarr(s[1])
  wave0=sxpar(h,'CRVAL1')
  step=sxpar(h,'CD1_1')
  if step eq 0 then step=sxpar(h,'CDELT1')
  
  ;convert flux units
  for zz=0,s[1]-1,1 do begin
    ;convert from Jy to 10^-17 erg/s/cm2/AA
    wave[zz]=wave0+zz*step
    wave_lin[zz]=10^wave[zz]
    input[zz]=input[zz]/(3.34e4*wave_lin[zz]*wave_lin[zz]*1e-17) 
  endfor

  ;print out datacube, after modifying header
  sxaddpar,h,'BUNIT','1E-17 erg/s/cm^2/Ang/spaxel', 'Specific intensity (per spaxel)' 
 

  split = STRSPLIT(spec[j], '/', /EXTRACT)

  fits_write,root+decomp+decomp_dir+'decomposed_data_ergs/'+split[-1],input[0:-cut],h
  
;  ;linearly rebin the data
;  lin10_rebin, [wave[0],wave[-1]], input, input_lin, lin_wave
;  split2=STRSPLIT(split[6], '.', /EXTRACT)
;  h_lin=h
;  
;  sxaddpar,h_lin,'CRVAL1',lin_wave[0]
;  sxaddpar,h_lin,'CDELT1',lin_wave[1]-lin_wave[0]
;  sxaddpar,h_lin,'CD1_1',lin_wave[1]-lin_wave[0]
;
;  fits_write,root+decomp+decomp_dir+'decomposed_data_ergs/'+split2[0]+'_linear.fits',input_lin,h_lin
  
  if split[-1] eq 'component1_flux.fits' then comp1=input;_lin
  if split[-1] eq 'component2_flux.fits' then comp2=input;_lin
  if split[-1] eq 'component3_flux.fits' then comp3=input;_lin
  if split[-1] eq 'component4_flux.fits' then comp4=input;_lin
endfor


s_cube=size(input);_lin)


for j=0,nfiles_cube-1,1 do begin
  ;read in datacube
  fits_read,cubes[j],input,h
  s=size(input)
  wave=fltarr(s[3])
  wave_lin=fltarr(s[3])

  
  ;convert flux units
  for zz=0,s[3]-1,1 do begin
    ;convert from Jy to 10^-17 erg/s/cm2/AA
    wave[zz]=sxpar(h,'CRVAL3')+zz*sxpar(h,'CD3_3')
    wave_lin[zz]=10^wave[zz]
    input[*,*,zz]=input[*,*,zz]/(3.34e4*wave_lin[zz]*wave_lin[zz]*1e-17)
  endfor
  
  ;print out datacube, after modifying header
;  temp=MRDFITS(cubes[j], 1, h_flux)
;  temp=MRDFITS(cubes[j], 0, h)
    h=headfits(cubes[j])
    if keyword_set(MANGA) then begin
      sxaddpar,h,'BUNIT','1E-17 erg/s/cm^2/Ang/spaxel', 'Specific intensity (per spaxel)'
      sxaddpar,h_flux,'BUNIT','1E-17 erg/s/cm^2/Ang/spaxel', 'Specific intensity (per spaxel)'
    endif else begin
      sxaddpar,h,'BUNIT','1E-20 erg/s/cm^2/Ang/spaxel', 'Specific intensity (per spaxel)'
      sxaddpar,h_flux,'BUNIT','1E-20 erg/s/cm^2/Ang/spaxel', 'Specific intensity (per spaxel)'      
    endelse
  
  split = STRSPLIT(cubes[j], '/', /EXTRACT)
  ;print,split[-1]

  fits_write,root+decomp+decomp_dir+'decomposed_data_ergs/'+split[-1],input[*,*,0:-cut],h_flux,extname='FLUX'
  modfits,root+decomp+decomp_dir+'decomposed_data_ergs/'+split[-1],0,h
;  modfits,root+decomp+decomp_dir+'decomposed_data_ergs/'+split[-1],1,h_flux,extname='FLUX'



  
;  linear_cube=fltarr(s[1],s[2],s_cube[1])
;  split2=STRSPLIT(split[6], '.', /EXTRACT)
;  for col=0,s[1]-1,1 do begin
;    for row=0,s[2]-1,1 do begin
;      ;linearly rebin the data
;      temp=fltarr(s[3])
;      temp[*]=input[col,row,*]
;      lin10_rebin, [wave[0],wave[-1]], temp, spec_lin, lin_wave
;      linear_cube[col,row,*]=spec_lin
;;      if col eq 150 and row eq 150 then stop
;    endfor
;  endfor
;;  h_lin=h_flux
;  h_lin=headfits(root+decomp+decomp_dir+'bestfit_cube.fits',EXTEN=1)
;  sxaddpar,h_lin,'CRVAL3',lin_wave[0]
;  sxaddpar,h_lin,'CDELT3',lin_wave[1]-lin_wave[0]
;  sxaddpar,h_lin,'CD3_3',lin_wave[1]-lin_wave[0]
;  sxaddpar,h,'CRVAL3',lin_wave[0]
;  sxaddpar,h,'CDELT3',lin_wave[1]-lin_wave[0]
;  sxaddpar,h,'CD3_3',lin_wave[1]-lin_wave[0]
;  
;  fits_write,root+decomp+decomp_dir+'decomposed_data_ergs/'+split2[0]+'_linear.fits',linear_cube,h_lin;,extname='FLUX'
;;  modfits,root+decomp+decomp_dir+'decomposed_data_ergs/'+split2[0]+'_linear.fits',0,h
;;  modfits,root+decomp+decomp_dir+'decomposed_data_ergs/'+split2[0]+'_linear.fits',1,h_lin,extname='FLUX'
;stop
  fits_read,root+decomp+setup.median_dir+'badpix.fits',badpix,h_bp
  badpix=abs(badpix-1)
  
  if split[-1] eq 'component1_cube.fits' then begin
    comp1_cube=fltarr(zz)
    for n=0,zz-1,1 do comp1_cube[n]=total(input[*,*,n]*badpix[*,*])
  endif
  if split[-1] eq 'component2_cube.fits' then begin
    comp2_cube=fltarr(zz)
    for n=0,zz-1,1 do comp2_cube[n]=total(input[*,*,n]*badpix[*,*])
  endif
  if split[-1] eq 'component3_cube.fits' then begin
    comp3_cube=fltarr(zz)
    for n=0,zz-1,1 do comp3_cube[n]=total(input[*,*,n]*badpix[*,*])
  endif
  if split[-1] eq 'original_cube.fits' then begin
    orig=fltarr(zz)
    for n=0,zz-1,1 do orig[n]=total(input[*,*,n]*badpix[*,*])
    
  endif
  if split[-1] eq 'bestfit_cube.fits' then begin
    bestfit=fltarr(zz)
    for n=0,zz-1,1 do bestfit[n]=total(input[*,*,n]*badpix[*,*])
  endif
  if split[-1] eq 'residuals_cube.fits' then begin
    resid=fltarr(zz)
    for n=0,zz-1,1 do resid[n]=total(input[*,*,n]*badpix[*,*])
  endif
  if split[-1] eq 'residual_sky_cube.fits' then begin
    sky=fltarr(zz)
    for n=0,zz-1,1 do sky[n]=total(input[*,*,n]*badpix[*,*])
  endif

endfor





set_plot,'ps'
device,file=root+decomp+decomp_dir+'decomposed_data_ergs/Spectra_integrated.eps',xsize=10,ysize=7,/inches,/color;,/landscape
!P.thick=2
!p.charthick=3
!p.charsize=1.3
!p.multi=0;[0,1,4]
;start_wavelength=4600
;end_wavelength=10000;7000

wavelength=10^wave
plot,wavelength[0:-cut],comp1[0:-cut],/NODATA,yrange=[-0.1,2.1],$
  xrange=[wavelength[0]-100,wavelength[-1]+100],$
  /xstyle,/ystyle,xthick=2,ythick=2,$;ytickinterval=30,$
  ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
  xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux',title=galaxy_ref
  
  
  oplot,wavelength,(orig[0:-cut]/median(orig[0:-cut]));/10000;+30
  oplot,wavelength,(resid[0:-cut]/median(orig[0:-cut])),color=cgcolor('olive');/10000,color=cgcolor('green')
  oplot,wavelength,(sky[0:-cut]/median(orig[0:-cut])),color=cgcolor('dark grey');/10000,color=cgcolor('green')
  ;oplot,wavelength,(bestfit[0:-cut]/median(orig[0:-cut])),color=cgcolor('purple');/10000,color=cgcolor('green')

  if setup.n_comp eq 1000 then begin
    oplot,wavelength,(comp1[0:-cut]/median(comp1[0:-cut])),color=cgcolor('blue');/10000;+90
    al_legend,['Integrated spectrum from datacube','Comp1','Sky','Residuals'],linestyle=[0,0,0,0],$
      colors=[cgcolor('black'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
  endif
  if setup.n_comp eq 1100 then begin
    oplot,wavelength,(comp1[0:-cut]/median(comp1+comp2[0:-cut])),color=cgcolor('blue');/10000;+90
    oplot,wavelength,(comp2[0:-cut]/median(comp1+comp2[0:-cut])),color=cgcolor('red');/10000;+90
    oplot,wavelength,((comp1[0:-cut]+comp2[0:-cut])/median(comp1+comp2[0:-cut])),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
    al_legend,['Integrated spectrum from datacube','Comp1+Comp2','Comp1','Comp2','Sky','Residuals'],linestyle=[0,0,0,0,0,0],$
      colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('red'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
  endif

 
  
  
  
    

  Ha_new=6563+(6563*setup.Redshift)
  Hb_new=4861+(4861*setup.Redshift)
  Hg_new=4341+(4341*setup.Redshift)
  Hd_new=4102+(4102*setup.Redshift)
  Mgb_new=5177+(5177*setup.Redshift)
  Fe5335_new=5335+(5335*setup.Redshift)
  Fe5270_new=5270+(5270*setup.Redshift)
  Na_new= 5895.92+(5895.92*setup.Redshift)
  CaII_new=3934+(3934*setup.Redshift)
  CaII2_new= 8542.09+(8542.09*setup.Redshift)
  Pae_new= 9546+(9546*setup.Redshift)
  Paz_new= 9229+(9229*setup.Redshift)
  
  
  
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
  
  if a[0]-30 gt 0 and a[0]+30 lt n_elements(orig)-1 then $
    aa=mean((orig[a[0]-30:a[0]+30]/median(orig)))+0.2 $
  else aa=-1000
  if b[0]-30 gt 0 and b[0]+30 lt n_elements(orig)-1 then $
    bb=mean((orig[b[0]-30:b[0]+30]/median(orig)))+0.2 $
  else bb=-1000
  if c[0]-30 gt 0 and c[0]+30 lt n_elements(orig)-1 then $
    cc=mean((orig[c[0]-30:c[0]+30]/median(orig)))+0.2 $
  else cc=-1000
  if d[0]-30 gt 0 and d[0]+30 lt n_elements(orig)-1 then $
    dd=mean((orig[d[0]-30:d[0]+30]/median(orig)))+0.2 $
  else dd=-1000
  if e[0]-30 gt 0 and e[0]+30 lt n_elements(orig)-1 then $
    ee=mean((orig[e[0]-30:e[0]+30]/median(orig)))+0.2 $
  else ee=-1000
  if f[0]-30 gt 0 and f[0]+30 lt n_elements(orig)-1 then $
    ff=mean((orig[f[0]-30:f[0]+30]/median(orig)))+0.2 $
  else ff=-1000
  if g[0]-30 gt 0 and g[0]+30 lt n_elements(orig)-1 then $
    gg=mean((orig[g[0]-30:g[0]+30]/median(orig)))+0.2 $
  else gg=-1000
  if h[0]-30 gt 0 and h[0]+30 lt n_elements(orig)-1 then $
    hh=mean((orig[h[0]-30:h[0]+30]/median(orig)))+0.2 $
  else hh=-1000
  if i[0]-30 gt 0 and i[0]+30 lt n_elements(orig)-1 then $
    ii=mean((orig[i[0]-30:i[0]+30]/median(orig)))+0.2 $
  else ii=-1000
  if j[0]-30 gt 0 and j[0]+30 lt n_elements(orig)-1 then $
    jj=mean((orig[j[0]-30:j[0]+30]/median(orig)))+0.2 $
  else jj=-1000
  if k[0]-30 gt 0 and k[0]+30 lt n_elements(orig)-1 then $
    kk=mean((orig[k[0]-30:k[0]+30]/median(orig)))+0.2 $
  else kk=-1000
  if l[0]-30 gt 0 and l[0]+30 lt n_elements(orig)-1 then $
    ll=mean((orig[l[0]-30:l[0]+30]/median(orig)))+0.2 $
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
;
;  device,file=root+decomp+decomp_dir+'decomposed_data_ergs/Spectra_integrated_masked.eps',/landscape;,xsize=11,ysize=8,/inches,/color;,/landscape
;  !P.thick=2
;  !p.charthick=3
;  !p.charsize=1.3
;  !p.multi=0;[0,1,4]
;  ;start_wavelength=4600
;  ;end_wavelength=10000;7000
;  
;  
;  wavelength=10^wave
;  plot,wavelength,comp1_cube,/NODATA,yrange=[-0.1,2.1],$
;    xrange=[wavelength[0]-100,wavelength[-1]+100],$
;    /xstyle,/ystyle,xthick=2,ythick=2,$;ytickinterval=30,$
;    ; ytickname=['Residuals','Galaxy + !CBest Fit','Disc','Bulge'],$
;    xtitle='Wavelength ('+cgSymbol("angstrom")+')',ytitle='Relative Flux',title=galaxy_ref
;    
;    
;  oplot,wavelength,(orig/median(orig));/10000;+30
;  oplot,wavelength,(resid/median(orig)),color=cgcolor('olive');/10000,color=cgcolor('green')
;  oplot,wavelength,(sky/median(orig)),color=cgcolor('dark grey');/10000,color=cgcolor('green')
;  
;  if setup.n_comp eq 1000 then begin
;    oplot,wavelength,(comp1_cube/median(comp1_cube)),color=cgcolor('blue');/10000;+90
;    al_legend,['Integrated spectrum from datacube','Comp1','Comp1','Sky','Residuals'],linestyle=[0,0,0,0,0],$
;      colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;  endif
;  if setup.n_comp eq 1100 then begin
;    oplot,wavelength,(comp1_cube/median(comp1_cube+comp2_cube)),color=cgcolor('blue');/10000;+90
;    oplot,wavelength,(comp2_cube/median(comp1_cube+comp2_cube)),color=cgcolor('red');/10000;+90
;    oplot,wavelength,((comp1_cube+comp2_cube)/median(comp1_cube+comp2_cube)),color=cgcolor('purple');/10000;+30,color=cgcolor('red')
;    al_legend,['Integrated spectrum from datacube','Comp1 + Comp2','Comp1','Comp2','Sky','Residuals'],linestyle=[0,0,0,0,0,0],$
;      colors=[cgcolor('black'),cgcolor('purple'),cgcolor('blue'),cgcolor('red'),cgcolor('dark grey'),cgcolor('olive')],charsize=0.9,box=0,/left,/top
;  endif
;
;  
;  
;  
;  
;  
;  
;
;;  xyouts,Ha_new-30,aa,'H'+greek('alpha'),charsize=0.9
;;  xyouts,Hb_new-30,bb,'H'+greek('beta'),charsize=0.9
;;  xyouts,Hg_new-30,cc,'H'+greek('gamma'),charsize=0.9
;;  xyouts,Hd_new-30,dd,'H'+greek('delta'),charsize=0.9
;;  xyouts,Mgb_new-30,ee,'Mg',charsize=0.9
;;  xyouts,Fe5335_new-30,ff,'Fe',charsize=0.9
;;  xyouts,Fe5270_new-30,gg,'Fe',charsize=0.9
;;  xyouts,Na_new-30,hh,'Na D',charsize=0.9
;;  ;xyouts,CaII_new-30,0.95,'Ca II',charsize=0.9
;;  xyouts,CaII_new-30,ii,'Ca II',charsize=0.9
;;  xyouts,CaII2_new-30,jj,'Ca II',charsize=0.9
;   
;
;  device,/close
;  


delvarx,wavelength,comp1_cube,comp2_cube,input,orig
delvarx,sky,resid,bestfit





end
