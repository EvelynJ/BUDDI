; 
; This code will print out the GalfitM feedme files for multi- 
; band fits for IFU data. Note that GALFITM is being used fo all 
; fits to maintain consistency in the results and in the codes.
; 
;  BINNED keyword- to be used with the binned images
;  SLICES keyword- to be used with individual image slices
; 
; scale
;
pro galfitm_multiband,setup,info,x,y,scale,$
  estimates_bulge,estimates_disk,estimates_comp3,estimates_comp4,$
  rep,BINNED=binned,SLICES=slices,FILE=file,HEADER=header
  
  
  root=setup.root
  decomp=setup.decomp
  binned_dir=setup.binned_dir
  slices_dir=setup.slices_dir
  median_dir=setup.median_dir
  galaxy_ref=setup.galaxy_ref
  psf_cube=setup.psf_file
  x_centre=fix(setup.x_centre-1)             ;x position of centre of galaxy, -1 to convert to position in array
  y_centre=fix(setup.y_centre-1)             ;y position of centre of galaxy, -1 to convert to position in array
  magzpt_in=setup.magzpt
  n_comp=setup.n_comp
  no_slices=setup.no_slices
  stars_file=setup.stars_file
  galfitm=setup.galfitm
  disk_re_polynomial_in=setup.disk_re_polynomial
  disk_mag_polynomial_in=setup.disk_mag_polynomial
  disk_n_polynomial_in=setup.disk_n_polynomial
  bulge_re_polynomial_in=setup.bulge_re_polynomial
  bulge_mag_polynomial_in=setup.bulge_mag_polynomial
  bulge_n_polynomial_in=setup.bulge_n_polynomial
  comp3_re_polynomial_in=setup.comp3_re_polynomial
  comp3_mag_polynomial_in=setup.comp3_mag_polynomial
  comp3_n_polynomial_in=setup.comp3_n_polynomial



  
output=root+decomp
first_image=info[0]
final_image=info[1]
no_bins=info[2]
images_per_bin=info[3]
start_wavelength=info[4]
end_wavelength=info[5]

ugriz=[4770,6231,7625,9134]
PSF_files=['Gpsf.fits','Rpsf.fits','Ipsf.fits','Zpsf.fits']

if estimates_disk[0] eq 0 then disk_type='sersic' else disk_type='psf'
if estimates_bulge[0] eq 0 then bulge_type='sersic' else bulge_type='psf'
if estimates_comp3[0] eq 0 then comp3_type='sersic' else comp3_type='psf'
if estimates_comp4[0] eq 0 then comp4_type='sersic' else comp4_type='psf'


if keyword_set(binned) then begin
  if keyword_set(file) then input='file'    ;'header' or input 'file'
  if keyword_set(header) then input='header'    ;'header' or input 'file'
  ;sometimes GalfitM cannot manage a good fit to the median image, 
  ;generally resulting in a sersic index of 20. This value is too 
  ;far wrong for Galfitm to constrain it when fitting the binned 
  ;images.
  
  x1=no_bins
  nband=1;info[2]
  
  if median_dir eq binned_dir then imgblock='imgblock_free' $
    else if median_dir ne binned_dir and n_comp eq 1000 then imgblock='imgblock_single' $
    else imgblock='imgblock_double'
  
  if n_comp eq 1000 then res=read_sersic_results_2comp(output+median_dir+imgblock+'.fits', nband, bd=0) $
  else if n_comp eq 1100 then res=read_sersic_results_2comp(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'psf' then res=read_sersic_results_3psf(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1001 and comp4_type eq 'psf' then res=read_sersic_results_3psf(output+median_dir+imgblock+'.fits', nband, bd=0) $
  else if n_comp eq 1001 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(output+median_dir+imgblock+'.fits', nband, bd=0) $
  
  else if n_comp eq 1010  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(output+median_dir+imgblock+'.fits', nband, bd=0) $
  else if n_comp eq 1010  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(output+median_dir+imgblock+'.fits', nband, bd=0) $
  else if n_comp eq 1110  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1110  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(output+median_dir+imgblock+'.fits', nband, bd=1) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(output+median_dir+imgblock+'.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(output+median_dir+imgblock+'.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(output+median_dir+imgblock+'.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(output+median_dir+imgblock+'.fits', nband, bd=0) 

  if disk_mag_polynomial_in lt 0 then disk_mag_polynomial=x1 else disk_mag_polynomial=disk_mag_polynomial_in
  if bulge_mag_polynomial_in lt 0 then bulge_mag_polynomial=x1 else bulge_mag_polynomial=bulge_mag_polynomial_in
  if comp3_mag_polynomial_in lt 0 and comp3_type eq 'sersic' then comp3_mag_polynomial=x1+1 else comp3_mag_polynomial=comp3_mag_polynomial_in

  if disk_re_polynomial_in lt 0 then disk_re_polynomial=x1 else disk_re_polynomial=disk_re_polynomial_in
  if disk_n_polynomial_in lt 0 then disk_n_polynomial=x1 else disk_n_polynomial=disk_n_polynomial_in
  if bulge_re_polynomial_in lt 0 then bulge_re_polynomial=x1 else bulge_re_polynomial=bulge_re_polynomial_in
  if bulge_n_polynomial_in lt 0 then bulge_n_polynomial=x1 else bulge_n_polynomial=bulge_n_polynomial_in
  
  if comp3_re_polynomial_in lt 0 and comp3_type eq 'sersic' then comp3_re_polynomial=x1 else comp3_re_polynomial=bulge_n_polynomial_in
  if comp3_n_polynomial_in lt 0 and comp3_type eq 'sersic' then comp3_n_polynomial=x1 else comp3_n_polynomial=comp3_n_polynomial_in


  
  
;  test= file_search(output+binned_dir+'imgblock*',COUNT=nfiles1)
;  if nfiles1 gt 0 then spawn,'rm '+output+binned_dir+'imgblock*'
  
  ;make string arrays for each galfitm input parameter
  n=0
  fits_read,output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits',tempy,h
;  h=headfits(output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits')
  x_size=sxpar(h,'NAXIS1')
  y_size=sxpar(h,'NAXIS2')


  if disk_n_polynomial_in eq 0 then res.N_GALFIT_BAND_D[*]=estimates_disk[3]
  if disk_Re_polynomial_in eq 0 then res.RE_GALFIT_BAND_D[*]=estimates_disk[2]
  if bulge_n_polynomial_in eq 0 then res.N_GALFIT_BAND_B[*]=estimates_bulge[3]
  if bulge_Re_polynomial_in eq 0 then res.RE_GALFIT_BAND_B[*]=estimates_bulge[2]
  
  if input eq 'header' and rep ne 1 then begin
      ;determine the psf file to be used first.
      temp=abs(ugriz-sxpar(h,'WAVELENG'))
      m=where(temp eq min(temp))
      ;psf_temp=PSF_files[m]

      file='image_'+string(n,format='(I4.4)')+'.fits'
      band=string(n,format='(I3.3)')
      wavelength=string(sxpar(h,'WAVELENG'),format='(F09.3)')
      psf='PSF/'+string(n,format='(I4.4)')+'.fits';psf_temp
      badpix='badpix_end.fits'
      temp2=file_search(output+binned_dir+'sigma*.fits',COUNT=nfiles_sig)
      if nfiles_sig gt 0 then sigma='sigma_0000.fits' else sigma=0
      
      magzpt=string(magzpt_in,format='(F04.1)')
      sky=string(mean(res.SKY_GALFIT_BAND),format='(F010.0)')
      sky_grad='0.0'
      x_D=string(mean(res.X_GALFIT_BAND_D),format='(F07.2)')
      y_D=string(mean(res.Y_GALFIT_BAND_D),format='(F07.2)')
      mag_D=string(mean(res.MAG_GALFIT_BAND_D),format='(F05.2)')
      Re_D=string(mean(res.RE_GALFIT_BAND_D),format='(F07.2)')
      n_D=string(mean(res.N_GALFIT_BAND_D),format='(F06.2)')
      q_D=string(mean(res.Q_GALFIT_BAND_D),format='(F04.2)')
      pa_D=string(mean(res.PA_GALFIT_BAND_D),format='(F06.2)')
      
      if n_comp ge 1100 then begin
        x_B=string(mean(res.X_GALFIT_BAND_B),format='(F07.2)')
        y_B=string(mean(res.Y_GALFIT_BAND_B),format='(F07.2)')
        mag_B=string(mean(res.MAG_GALFIT_BAND_B),format='(F05.2)')
        Re_B=string(mean(res.RE_GALFIT_BAND_B),format='(F07.2)')
        n_B=string(mean(res.N_GALFIT_BAND_B),format='(F06.2)')
        q_B=string(mean(res.Q_GALFIT_BAND_B),format='(F04.2)')
        pa_B=string(mean(res.PA_GALFIT_BAND_B),format='(F06.2)')
      endif
      
      
      if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin        
        x_comp3=string(mean(res.x_galfit_band_comp3),format='(F07.2)')
        y_comp3=string(mean(res.y_galfit_band_comp3),format='(F07.2)')
        mag_comp3=string(mean(res.mag_galfit_band_comp3),format='(F05.2)')
        if comp3_type eq 'sersic' then begin
            Re_comp3=string(mean(res.RE_galfit_band_comp3),format='(F06.2)')
            n_comp3=string(mean(res.N_galfit_band_comp3),format='(F06.2)')
            q_comp3=string(mean(res.Q_galfit_band_comp3),format='(F04.2)')
            pa_comp3=string(mean(res.PA_galfit_band_comp3),format='(F06.2)')
        endif  
      endif
      if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
        x_comp4=string(mean(res.x_galfit_band_comp4),format='(F07.2)')
        y_comp4=string(mean(res.y_galfit_band_comp4),format='(F07.2)')
        mag_comp4=string(mean(res.mag_galfit_band_comp4),format='(F05.2)')
        if comp4_type eq 'sersic' then begin
            Re_comp4=string(mean(res.RE_galfit_band_comp4),format='(F06.2)')
            n_comp4=string(mean(res.N_galfit_band_comp4),format='(F06.2)')
            q_comp4=string(mean(res.Q_galfit_band_comp4),format='(F04.2)')
            pa_comp4=string(mean(res.PA_galfit_band_comp4),format='(F06.2)')
        endif  
      endif


      for n=1,no_bins-1,1 do begin
        temp=abs(ugriz-sxpar(h,'WAVELENG'))
        m=where(temp eq min(temp))
        psf_temp=PSF_files[m]
        
        
        fits_read,output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits',crap,h
;        h=headfits(output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits')
        file+=',image_'+string(n,format='(I4.4)')+'.fits'
        band+=','+string(n,format='(I3.3)')
        wavelength+=','+string(sxpar(h,'WAVELENG'),format='(F09.3)')
        psf+=',PSF/'+string(n,format='(I4.4)')+'.fits'
        ;if n ne no_bins-1 then badpix=badpix+',badpix.fits' else badpix=badpix+',badpix_end.fits'
        
        result=file_search(root+decomp+binned_dir+'badpix*.fits',COUNT=nfiles_bp)
        if n ne no_bins-1 and nfiles_bp gt 2 then badpix=badpix+',badpix_'+string(n,format='(I4.4)')+'.fits' $
        else if n ne no_bins-1 and nfiles_bp le 2 then badpix=badpix+',badpix.fits' $
        else badpix=badpix+',badpix_end.fits'

        if nfiles_sig gt 0 then sigma=sigma+',sigma_'+string(n,format='(I4.4)')+'.fits'         
        magzpt+=','+string(magzpt_in,format='(F04.1)')
        sky+=','+string(mean(res.SKY_GALFIT_BAND),format='(F010.0)')
        sky_grad+=',0.0'
        x_D+=','+string(mean(res.X_GALFIT_BAND_D),format='(F07.2)')
        y_D+=','+string(mean(res.Y_GALFIT_BAND_D),format='(F07.2)')
        mag_D+=','+string(mean(res.MAG_GALFIT_BAND_D),format='(F05.2)')
        Re_D+=','+string(mean(res.RE_GALFIT_BAND_D),format='(F07.2)')
        n_D+=','+string(mean(res.N_GALFIT_BAND_D),format='(F05.2)')
        q_D+=','+string(mean(res.Q_GALFIT_BAND_D),format='(F04.2)')
        pa_D+=','+string(mean(res.PA_GALFIT_BAND_D),format='(F06.2)')
        
        ;insert parameters for bulge
        if n_comp ge 1100 then begin
          x_B+=','+string(mean(res.X_GALFIT_BAND_B),format='(F07.2)')
          y_B+=','+string(mean(res.Y_GALFIT_BAND_B),format='(F07.2)')
          mag_B+=','+string(mean(res.MAG_GALFIT_BAND_B),format='(F05.2)')
          Re_B+=','+string(mean(res.RE_GALFIT_BAND_B),format='(F07.2)')
          n_B+=','+string(mean(res.N_GALFIT_BAND_B),format='(F05.2)')
          q_B+=','+string(mean(res.Q_GALFIT_BAND_B),format='(F04.2)')
          pa_B+=','+string(mean(res.PA_GALFIT_BAND_B),format='(F06.2)')
        endif
        
        ;insert parameters for 3rd galaxy component
        if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
          mag_comp3+=','+string(mean(res.mag_galfit_band_comp3),format='(F05.2)')
          x_comp3+=','+string(mean(res.x_galfit_band_comp3),format='(F07.2)')
          y_comp3+=','+string(mean(res.y_galfit_band_comp3),format='(F07.2)')
          if comp3_type eq 'sersic' then begin
            Re_comp3+=','+string(mean(res.RE_galfit_band_comp3),format='(F06.2)')
            n_comp3+=','+string(mean(res.N_galfit_band_comp3),format='(F06.2)')
            q_comp3+=','+string(mean(res.Q_galfit_band_comp3),format='(F04.2)')
            pa_comp3+=','+string(mean(res.PA_galfit_band_comp3),format='(F06.2)') 
          endif
        endif
        ;insert parameters for 4th component
        if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
          mag_comp4+=','+string(mean(res.mag_galfit_band_comp4),format='(F05.2)')
          x_comp4+=','+string(mean(res.x_galfit_band_comp4),format='(F07.2)')
          y_comp4+=','+string(mean(res.y_galfit_band_comp4),format='(F07.2)')
          if comp4_type eq 'sersic' then begin
            Re_comp4+=','+string(mean(res.RE_galfit_band_comp4),format='(F06.2)')
            n_comp4+=','+string(mean(res.N_galfit_band_comp4),format='(F06.2)')
            q_comp4+=','+string(mean(res.Q_galfit_band_comp4),format='(F04.2)')
            pa_comp4+=','+string(mean(res.PA_galfit_band_comp4),format='(F06.2)') 
          endif
        endif
          
      endfor  
      
    endif else if input eq 'header' and rep eq 1 then begin
      ;determine the psf file to be used first.
      temp=abs(ugriz-sxpar(h,'WAVELENG'))
      m=where(temp eq min(temp))
      psf_temp=PSF_files[m]

      file='image_'+string(n,format='(I4.4)')+'.fits'
      band=string(n,format='(I3.3)')
      wavelength=string(sxpar(h,'WAVELENG'),format='(F09.3)')
      psf='PSF/'+string(n,format='(I4.4)')+'.fits'
      ;psf='psf_'+string(n,format='(I4.4)')+'.fits'
      badpix='badpix_end.fits'
      temp2=file_search(output+binned_dir+'sigma*.fits',COUNT=nfiles_sig)
      if nfiles_sig gt 0 then sigma='sigma_0000.fits' else sigma=0

      
      magzpt=string(magzpt_in,format='(F04.1)')
      sky=string((res.SKY_GALFIT_BAND),format='(F010.0)')
      sky_grad='0.0'
      x_D=string((res.X_GALFIT_BAND_D),format='(F07.2)')
      y_D=string((res.Y_GALFIT_BAND_D),format='(F07.2)')
      mag_D=string((res.MAG_GALFIT_BAND_D),format='(F05.2)')
      Re_D=string((res.RE_GALFIT_BAND_D),format='(F07.2)')
      n_D=string((res.N_GALFIT_BAND_D),format='(F06.2)')
      q_D=string((res.Q_GALFIT_BAND_D),format='(F04.2)')
      pa_D=string((res.PA_GALFIT_BAND_D),format='(F06.2)')
      
      if n_comp ge 1100 then begin
        x_B=string((res.X_GALFIT_BAND_B),format='(F07.2)')
        y_B=string((res.Y_GALFIT_BAND_B),format='(F07.2)')
        mag_B=string((res.MAG_GALFIT_BAND_B),format='(F05.2)')
        Re_B=string((res.RE_GALFIT_BAND_B),format='(F07.2)')
        n_B=string((res.N_GALFIT_BAND_B),format='(F06.2)')
        q_B=string((res.Q_GALFIT_BAND_B),format='(F04.2)')
        pa_B=string((res.PA_GALFIT_BAND_B),format='(F06.2)')
      endif
      
      
      if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin        
        x_comp3=string((res.x_galfit_band_comp3),format='(F07.2)')
        y_comp3=string((res.y_galfit_band_comp3),format='(F07.2)')
        mag_comp3=string((res.mag_galfit_band_comp3),format='(F05.2)')
        if comp3_type eq 'sersic' then begin
            Re_comp3=string((res.RE_galfit_band_comp3),format='(F06.2)')
            n_comp3=string((res.N_galfit_band_comp3),format='(F06.2)')
            q_comp3=string((res.Q_galfit_band_comp3),format='(F04.2)')
            pa_comp3=string((res.PA_galfit_band_comp3),format='(F06.2)')
        endif  
      endif
      if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
        x_comp4=string((res.x_galfit_band_comp4),format='(F07.2)')
        y_comp4=string((res.y_galfit_band_comp4),format='(F07.2)')
        mag_comp4=string((res.mag_galfit_band_comp4),format='(F05.2)')
        if comp4_type eq 'sersic' then begin
            Re_comp4=string((res.RE_galfit_band_comp4),format='(F06.2)')
            n_comp4=string((res.N_galfit_band_comp4),format='(F06.2)')
            q_comp4=string((res.Q_galfit_band_comp4),format='(F04.2)')
            pa_comp4=string((res.PA_galfit_band_comp4),format='(F06.2)')
        endif  
      endif


      for n=1,no_bins-1,1 do begin
        temp=abs(ugriz-sxpar(h,'WAVELENG'))
        m=where(temp eq min(temp))
        psf_temp=PSF_files[m]
        
        
        fits_read,output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits',crap,h
;        h=headfits(output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits')
        file+=',image_'+string(n,format='(I4.4)')+'.fits'
        band+=','+string(n,format='(I3.3)')
        wavelength+=','+string(sxpar(h,'WAVELENG'),format='(F09.3)')
        ;print,n,sxpar(h,'WAVELENG')
        psf+=',PSF/'+string(n,format='(I4.4)')+'.fits'
;        psf+=',psf_'+string(n,format='(I4.4)')+'.fits'
        
        result=file_search(root+decomp+binned_dir+'badpix*.fits',COUNT=nfiles_bp)
        if n ne no_bins-1 and nfiles_bp gt 2 then badpix=badpix+',badpix_'+string(n,format='(I4.4)')+'.fits' $
        else if n ne no_bins-1 and nfiles_bp le 2 then badpix=badpix+',badpix.fits' $
        else badpix=badpix+',badpix_end.fits'
        if nfiles_sig gt 0 then sigma=sigma+',sigma_'+string(n,format='(I4.4)')+'.fits' else sigma=0

        
        magzpt+=','+string(magzpt_in,format='(F04.1)')
        sky+=','+string((res.SKY_GALFIT_BAND),format='(F010.0)')
        sky_grad+=',0.0'
        x_D+=','+string((res.X_GALFIT_BAND_D),format='(F07.2)')
        y_D+=','+string((res.Y_GALFIT_BAND_D),format='(F07.2)')
        mag_D+=','+string((res.MAG_GALFIT_BAND_D),format='(F05.2)')
        Re_D+=','+string((res.RE_GALFIT_BAND_D),format='(F07.2)')
        n_D+=','+string((res.N_GALFIT_BAND_D),format='(F05.2)')
        q_D+=','+string((res.Q_GALFIT_BAND_D),format='(F04.2)')
        pa_D+=','+string((res.PA_GALFIT_BAND_D),format='(F06.2)')
        
        ;insert parameters for bulge
        if n_comp ge 1100 then begin
          x_B+=','+string((res.X_GALFIT_BAND_B),format='(F07.2)')
          y_B+=','+string((res.Y_GALFIT_BAND_B),format='(F07.2)')
          mag_B+=','+string((res.MAG_GALFIT_BAND_B),format='(F05.2)')
          Re_B+=','+string((res.RE_GALFIT_BAND_B),format='(F07.2)')
          n_B+=','+string((res.N_GALFIT_BAND_B),format='(F05.2)')
          q_B+=','+string((res.Q_GALFIT_BAND_B),format='(F04.2)')
          pa_B+=','+string((res.PA_GALFIT_BAND_B),format='(F06.2)')
        endif
        
        ;insert parameters for 3rd galaxy component
        if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
          mag_comp3+=','+string((res.mag_galfit_band_comp3),format='(F05.2)')
          x_comp3+=','+string((res.x_galfit_band_comp3),format='(F07.2)')
          y_comp3+=','+string((res.y_galfit_band_comp3),format='(F07.2)')
          if comp3_type eq 'sersic' then begin
            Re_comp3+=','+string((res.RE_galfit_band_comp3),format='(F06.2)')
            n_comp3+=','+string((res.N_galfit_band_comp3),format='(F06.2)')
            q_comp3+=','+string((res.Q_galfit_band_comp3),format='(F04.2)')
            pa_comp3+=','+string((res.PA_galfit_band_comp3),format='(F06.2)') 
          endif
        endif
        ;insert parameters for 4th component
        if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
          mag_comp4+=','+string((res.mag_galfit_band_comp4),format='(F05.2)')
          x_comp4+=','+string((res.x_galfit_band_comp4),format='(F07.2)')
          y_comp4+=','+string((res.y_galfit_band_comp4),format='(F07.2)')
          if comp4_type eq 'sersic' then begin
            Re_comp4+=','+string((res.RE_galfit_band_comp4),format='(F06.2)')
            n_comp4+=','+string((res.N_galfit_band_comp4),format='(F06.2)')
            q_comp4+=','+string((res.Q_galfit_band_comp4),format='(F04.2)')
            pa_comp4+=','+string((res.PA_galfit_band_comp4),format='(F06.2)') 
          endif
        endif
          
      endfor  
      
    endif else if input eq 'file' then begin
      
      file='image_'+string(n,format='(I4.4)')+'.fits'
      band=string(n,format='(I3.3)')
      wavelength=string(sxpar(h,'WAVELENG'),format='(F09.3)')
      psf='PSF/'+string(n,format='(I4.4)')+'.fits'
;      psf='psf_'+string(n,format='(I4.4)')+'.fits'
      badpix='badpix_end.fits'
      temp2=file_search(output+binned_dir+'sigma*.fits',COUNT=nfiles_sig)
      if nfiles_sig gt 0 then sigma='sigma_0000.fits' else sigma=0

      
      magzpt=string(magzpt_in,format='(F04.1)')
      sky=string(0,format='(F010.0)')
      sky_grad='0.0'
      x_D=string(x,format='(F07.2)')
      y_D=string(y,format='(F07.2)')
      mag_D=string(estimates_disk[1],format='(F05.2)')
      Re_D=string(estimates_disk[2],format='(F07.2)')
      n_D=string(estimates_disk[3],format='(F06.2)')
      q_D=string(estimates_disk[4],format='(F04.2)')
      pa_D=string(estimates_disk[5],format='(F06.2)')

      if n_comp ge 1100 then begin
          x_B=string(x,format='(F07.2)')
          y_B=string(y,format='(F07.2)')
          mag_B=string(estimates_bulge[1],format='(F05.2)')
          Re_B=string(estimates_bulge[2],format='(F07.2)')
          n_B=string(estimates_bulge[3],format='(F06.2)')
          q_B=string(estimates_bulge[4],format='(F04.2)')
          pa_B=string(estimates_bulge[5],format='(F06.2)')
      endif
      
       if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
          x_comp3=string(x,format='(F07.2)')
          y_comp3=string(y,format='(F07.2)')
          mag_comp3=string(estimates_comp3[1],format='(F05.2)')
          if comp3_type eq 'sersic' then begin
            Re_comp3=string(estimates_comp3[2],format='(F06.2)')
            n_comp3=string(estimates_comp3[3],format='(F06.2)')
            q_comp3=string(estimates_comp3[4],format='(F04.2)')
            pa_comp3=string(estimates_comp3[5],format='(F06.2)') 
          endif             
       endif
       if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
          x_comp4=string(estimates_comp4[1],format='(F07.2)')
          y_comp4=string(estimates_comp4[2],format='(F07.2)')
          mag_comp4=string(estimates_comp4[3],format='(F05.2)')
          if comp4_type eq 'sersic' then begin
            Re_comp4=string(estimates_comp4[4],format='(F06.2)')
            n_comp4=string(estimates_comp4[5],format='(F06.2)')
            q_comp4=string(estimates_comp4[6],format='(F04.2)')
            pa_comp4=string(estimates_comp4[7],format='(F06.2)') 
          endif             
       endif    
      
      

    
      for n=1,no_bins-1,1 do begin
        fits_read,output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits',crap,h
;        h=headfits(output+binned_dir+'image_'+string(n,format='(I4.4)')+'.fits')
        file+=',image_'+string(n,format='(I4.4)')+'.fits'
        band+=','+string(n,format='(I3.3)')
        wavelength+=','+string(sxpar(h,'WAVELENG'),format='(F09.3)')
        ;print,n,sxpar(h,'WAVELENG')
        psf+=',PSF/'+string(n,format='(I4.4)')+'.fits'
;        psf+=',psf_'+string(n,format='(I4.4)')+'.fits'
        result=file_search(root+decomp+binned_dir+'badpix*.fits',COUNT=nfiles_bp)
        if n ne no_bins-1 and nfiles_bp gt 2 then badpix=badpix+',badpix_'+string(n,format='(I4.4)')+'.fits' $
        else if n ne no_bins-1 and nfiles_bp le 2 then badpix=badpix+',badpix.fits' $
        else badpix=badpix+',badpix_end.fits'
        if nfiles_sig gt 0 then sigma=sigma+',sigma_'+string(n,format='(I4.4)')+'.fits'

        ;badpix+=',badpix.fits'
        magzpt+=','+string(magzpt_in,format='(F04.1)')
        sky+=','+string(0,format='(F010.0)')
        sky_grad+=',0.0'
        x_D+=','+string(x,format='(F07.2)')
        y_D+=','+string(y,format='(F07.2)')
        mag_D+=','+string(estimates_disk[1],format='(F05.2)')
        Re_D+=','+string(estimates_disk[2],format='(F07.2)')
        n_D+=','+string(estimates_disk[3],format='(F05.2)')
        q_D+=','+string(estimates_disk[4],format='(F04.2)')
        pa_D+=','+string(estimates_disk[5],format='(F06.2)')

        if n_comp ge 1100 then begin
            x_B+=','+string(x,format='(F07.2)')
            y_B+=','+string(y,format='(F07.2)')
            mag_B+=','+string(estimates_bulge[1],format='(F05.2)')
            Re_B+=','+string(estimates_bulge[2],format='(F07.2)')
            n_B+=','+string(estimates_bulge[3],format='(F05.2)')
            q_B+=','+string(estimates_bulge[4],format='(F04.2)')
            pa_B+=','+string(estimates_bulge[5],format='(F06.2)')
        endif 
        if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
          x_comp3+=','+string(x,format='(F07.2)')
          y_comp3+=','+string(y,format='(F07.2)')
          mag_comp3+=','+string(estimates_comp3[1],format='(F05.2)')
          if comp3_type eq 'sersic' then begin
            Re_comp3+=','+string(estimates_comp3[2],format='(F06.2)')
            n_comp3+=','+string(estimates_comp3[3],format='(F06.2)')
            q_comp3+=','+string(estimates_comp3[4],format='(F04.2)')
            pa_comp3+=','+string(estimates_comp3[5],format='(F06.2)')  
          endif            
        endif
        if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
          x_comp4+=','+string(estimates_comp4[1],format='(F07.2)')
          y_comp4+=','+string(estimates_comp4[2],format='(F07.2)')
          mag_comp4+=','+string(estimates_comp4[3],format='(F05.2)')
          if comp4_type eq 'sersic' then begin
            Re_comp4+=','+string(estimates_comp4[4],format='(F06.2)')
            n_comp4+=','+string(estimates_comp4[5],format='(F06.2)')
            q_comp4+=','+string(estimates_comp4[6],format='(F04.2)')
            pa_comp4+=','+string(estimates_comp4[7],format='(F06.2)')  
          endif            
        endif
        

;        endif  
      endfor  

      
      

    endif else begin
        Message,'Please select an input method for the Galfitm constraints'
    endelse
     

    if rep ne 1 then begin
      close,60
      openw,60,output+binned_dir+'galfitm.feedme'
      
      printf,60,'==============================================================================='
      printf,60,'# IMAGE and GALFIT CONTROL PARAMETERS';+$
      printf, 60, 'A) '+file+'             # Input data image (FITS file)'
      printf, 60, 'A1) '+band+'             # Band labels (can be omitted if fitting a single band)'
      printf, 60, 'A2) '+wavelength+'             # Band wavelengths'
      printf, 60, 'B) imgblock.fits       # Output data image block
      printf, 60, 'C) '+sigma+'                # Sigma image name (made from data if blank or "none") 
      printf, 60, 'D) '+psf+'           # Input PSF image and (optional) diffusion kernel
      printf, 60, 'E) 1                   # PSF fine sampling factor relative to data 
      printf, 60, 'F) '+badpix+'                # Bad pixel mask (FITS image or ASCII coord list)
      
      if n_comp ne 1000 and n_comp ne 1001 then printf, 60, 'G) galfitm.constraints                # File with parameter constraints (ASCII file)'  $
        else printf, 60, 'G) none                # File with parameter constraints (ASCII file)'  
      printf, 60, 'H) 1    '+string(x_size,format='(I4.4)')+'   1  '+string(y_size,format='(I4.4)')+'    # Image region to fit (xmin xmax ymin ymax)'
      printf, 60, 'I) '+string(x_size,format='(I4.4)')+'    '+string(x_size,format='(I4.4)')+'      # Size of the convolution box (x y)'
      printf, 60, 'J) '+magzpt+'              # Magnitude photometric zeropoint '
      printf, 60, 'K) '+string(scale[0],format='(F4.2)')+'    '+string(scale[1],format='(F4.2)')+'        # Plate scale (dx dy)    [arcsec per pixel]'
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
      printf, 60, '  1) '+sky+'   '+string(no_bins)+' band #  sky background at center of fitting region [ADUs]'
      printf, 60, '  2) '+sky_grad+'      0 band  #  dsky/dx (sky gradient in x)'
      printf, 60, '  3) '+sky_grad+'      0 band  #  dsky/dy (sky gradient in y)'
      printf, 60, '  Z) 0                      #  output option (0 = resid., 1 = Dont subtract) '
       
     
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      
       
      printf, 60, '# Object number: 2'    ;disc
      printf, 60, ' 0) '+disk_type+'                 #  object type'
      printf, 60, ' 1) '+x_D+'   1 band  #  position x, y'
      printf, 60, ' 2) '+y_D+'   1 band  #  position x, y'
      printf, 60, ' 3) '+mag_D+'        '+string(no_bins)+' band  #  Integrated magnitude' 
      printf, 60, ' 4) '+Re_D+'   '+string(disk_re_polynomial)+' band  #  R_e (half-light radius)   [pix]'
      printf, 60, ' 5) '+n_D+'             '+string(disk_n_polynomial)+' band  #  Sersic index n (de Vaucouleurs n=4) '
      printf, 60, ' 9) '+Q_D+'        1 band  #  axis ratio (b/a)  '
      printf, 60, '10) '+PA_D+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
      ;printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'
     print,disk_n_polynomial
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      
      if n_comp ge 1100 then begin
        printf, 60, ' # Object number: 3   '    ;bulge
        printf, 60, ' 0) '+bulge_type+'                 #  object type'
        printf, 60, ' 1) '+x_B+'   1 band  #  position x, y'
        printf, 60, ' 2) '+y_B+'   1 band  #  position x, y'
        printf, 60, ' 3) '+mag_B+'        '+string(no_bins)+' band  #  Integrated magnitude' 
        printf, 60, ' 4) '+Re_B+'   '+string(bulge_re_polynomial)+' band  #  R_e (half-light radius)   [pix]'
        printf, 60, ' 5) '+n_B+'             '+string(bulge_n_polynomial)+' band  #  Sersic index n (de Vaucouleurs n=4) '
        printf, 60, ' 9) '+q_B+'        1 band  #  axis ratio (b/a)  '
        printf, 60, '10) '+pa_B+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
        printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
        ;printf, 60, 'C0) 1.2         1      # traditional diskyness(-)/boxyness(+)'
    
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
      endif 
      
      if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
        
        printf, 60, ' 0) '+comp3_type+'                # object type'
        printf, 60, ' 1) '+x_comp3+'   1   #  position x, y'
        printf, 60, ' 2) '+y_comp3+'   1   #  position x, y'
        printf, 60, ' 3) '+mag_comp3+'       '+string(no_bins)+'       # total magnitude   '  
        if comp3_type eq 'sersic' then begin
          printf, 60, ' 4) '+Re_comp3+'   '+string(comp3_re_polynomial)+' band  #  R_e (half-light radius)   [pix]'
          printf, 60, ' 5) '+n_comp3+'    '+string(comp3_n_polynomial)+' band  #  Sersic index n (de Vaucouleurs n=4) '
          printf, 60, ' 9) '+q_comp3+'        1 band  #  axis ratio (b/a)  '
          printf, 60, '10) '+pa_comp3+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
        endif
        printf, 60, ' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
      endif
  
      if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
        printf, 60, ' 0) '+comp4_type+'                # object type'
        printf, 60, ' 1) '+x_comp4+'   1   #  position x, y'
        printf, 60, ' 2) '+y_comp4+'   1   #  position x, y'
        printf, 60, ' 3) '+mag_comp4+'       '+string(no_bins)+'       # total magnitude   '  
        if comp4_type eq 'sersic' then begin
          printf, 60, ' 4) '+Re_comp4+'   1 band  #  R_e (half-light radius)   [pix]'
          printf, 60, ' 5) '+n_comp4+'             1 band  #  Sersic index n (de Vaucouleurs n=4) '
          printf, 60, ' 9) '+q_comp4+'        1 band  #  axis ratio (b/a)  '
          printf, 60, '10) '+pa_comp4+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
        endif
        printf, 60, ' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
      endif
  
      if file_test(root+stars_file) eq 1 then begin
        readcol,root+stars_file,format='f,f,f',x_star,y_star,mag_star,comment='#',/SILENT
        for j=0,n_elements(x_star)-1,1 do begin
          x_pos=string(x_star[j],format='(F07.2)')
          y_pos=string(y_star[j],format='(F07.2)')
          mag_pos=string(mag_star[j],format='(F05.2)')
          
          for n=1,no_bins-1,1 do begin
            x_pos+=','+string(x_star[j],format='(F07.2)')
            y_pos+=','+string(y_star[j],format='(F07.2)')
            mag_pos+=','+string(mag_star[j],format='(F05.2)')
          endfor
          
          printf, 60, ' # Object number:  '+string(j)
          printf, 60, ' 0) psf                 #  object type'
          printf, 60, ' 1) '+x_pos+'   '+string(no_bins)+' band  #  position x'
          printf, 60, ' 2) '+y_pos+'   '+string(no_bins)+' band  #  position y'
          printf, 60, ' 3) '+mag_pos+'   '+string(no_bins)+' band  #  Integrated magnitude'
          printf, 60, ' '
          printf, 60, ' '
          printf, 60, ' '
        endfor
      endif
  
  
      
      printf, 60, '================================================================================'
      
      
      
      close,60
  
  
     
;================================================================================'
    
    endif else begin
      close,60
      openw,60,output+binned_dir+'galfitm.feedme'
      
      printf,60,'==============================================================================='
      printf,60,'# IMAGE and GALFIT CONTROL PARAMETERS';+$
      printf, 60, 'A) '+file+'             # Input data image (FITS file)'
      printf, 60, 'A1) '+band+'             # Band labels (can be omitted if fitting a single band)'
      printf, 60, 'A2) '+wavelength+'             # Band wavelengths'
      printf, 60, 'B) imgblock.fits       # Output data image block
      printf, 60, 'C) '+sigma+'                # Sigma image name (made from data if blank or "none") 
      printf, 60, 'D) '+psf+'           # Input PSF image and (optional) diffusion kernel
      printf, 60, 'E) 1                   # PSF fine sampling factor relative to data 
      printf, 60, 'F) '+badpix+'                # Bad pixel mask (FITS image or ASCII coord list)
      if n_comp ne 1000 and n_comp ne 1001 then printf, 60, 'G) galfitm.constraints                # File with parameter constraints (ASCII file)'  $
        else printf, 60, 'G) none                # File with parameter constraints (ASCII file)'  
      printf, 60, 'H) 1    '+string(x_size,format='(I4.4)')+'   1  '+string(y_size,format='(I4.4)')+'    # Image region to fit (xmin xmax ymin ymax)'
      printf, 60, 'I) '+string(x_size,format='(I4.4)')+'    '+string(x_size,format='(I4.4)')+'          # Size of the convolution box (x y)'
      printf, 60, 'J) '+magzpt+'              # Magnitude photometric zeropoint '
      printf, 60, 'K) '+string(scale[0],format='(F4.2)')+'    '+string(scale[1],format='(F4.2)')+'        # Plate scale (dx dy)    [arcsec per pixel]'
      printf, 60, 'O) regular             # Display type (regular, curses, both)'
      printf, 60, 'P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps'
     
      printf, 60, ' '
      printf, 60, ' '
      ; 
  ;    x1=60
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
      printf, 60, '  1) '+sky+'   '+string(no_bins)+' band #  sky background at center of fitting region [ADUs]'
      printf, 60, '  2) '+sky_grad+'      0 band  #  dsky/dx (sky gradient in x)'
      printf, 60, '  3) '+sky_grad+'      0 band  #  dsky/dy (sky gradient in y)'
      printf, 60, '  Z) 0                      #  output option (0 = resid., 1 = Dont subtract) '
       
     
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
       
      printf, 60, '# Object number: 2'    ;disc
      printf, 60, ' 0) '+disk_type+'                 #  object type'
      printf, 60, ' 1) '+x_D+'   1 band  #  position x, y'
      printf, 60, ' 2) '+y_D+'   1 band  #  position x, y'
      printf, 60, ' 3) '+mag_D+'        '+string(no_bins,format='(I3.3)')+' band  #  Integrated magnitude' 
  ;   k_n_polynomial=0
      printf, 60, ' 4) '+Re_D+'   '+string(no_bins,format='(I3.3)')+' band  #  R_e (half-light radius)   [pix]'
      printf, 60, ' 5) '+n_D+'             '+string(no_bins,format='(I3.3)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
      printf, 60, ' 9) '+Q_D+'        '+string(no_bins,format='(I3.3)')+' band  #  axis ratio (b/a)  '
      printf, 60, '10) '+PA_D+'   '+string(no_bins,format='(I3.3)')+' band  #  position angle (PA) [deg: Up=0, Left=90]'
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
      ;printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'
     
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
     
      if n_comp ge 1100 then begin
        printf, 60, ' # Object number: 3   '    ;bulge
        printf, 60, ' 0) '+bulge_type+'                 #  object type'
        printf, 60, ' 1) '+x_B+'   1 band  #  position x, y'
        printf, 60, ' 2) '+y_B+'   1 band  #  position x, y'
        printf, 60, ' 3) '+mag_B+'        '+string(no_bins,format='(I3.3)')+' band  #  Integrated magnitude' 
        printf, 60, ' 4) '+Re_B+'   '+string(no_bins,format='(I3.3)')+' band  #  R_e (half-light radius)   [pix]'
        printf, 60, ' 5) '+n_B+'             '+string(no_bins,format='(I3.3)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
        printf, 60, ' 9) '+q_B+'        '+string(no_bins,format='(I3.3)')+' band  #  axis ratio (b/a)  '
        printf, 60, '10) '+pa_B+'   '+string(no_bins,format='(I3.3)')+' band  #  position angle (PA) [deg: Up=0, Left=90]'
        printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
        ;printf, 60, 'C0) 1.2         '+string(no_bins,format='(I3.3)')+'      # traditional diskyness(-)/boxyness(+)'
      ;  printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
      endif  
        
      if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
        printf, 60, ' 0) '+comp3_type+'                # object type'
        printf, 60, ' 1) '+x_comp3+'   1   #  position x, y'
        printf, 60, ' 2) '+y_comp3+'   1   #  position x, y'
        printf, 60, ' 3) '+mag_comp3+'       '+string(no_bins,format='(I3.3)')+'       # total magnitude   '  
        if comp3_type eq 'sersic' then begin
          printf, 60, ' 4) '+Re_comp3+'   '+string(no_bins,format='(I3.3)')+' band  #  R_e (half-light radius)   [pix]'
          printf, 60, ' 5) '+n_comp3+'             '+string(no_bins,format='(I3.3)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
          printf, 60, ' 9) '+q_comp3+'        '+string(no_bins,format='(I3.3)')+' band  #  axis ratio (b/a)  '
          printf, 60, '10) '+pa_comp3+'   '+string(no_bins,format='(I3.3)')+' band  #  position angle (PA) [deg: Up=0, Left=90]'
        endif
        printf, 60, ' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
      endif    
      if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
        printf, 60, ' 0) '+comp4_type+'                # object type'
        printf, 60, ' 1) '+x_comp4+'   1   #  position x, y'
        printf, 60, ' 2) '+y_comp4+'   1   #  position x, y'
        printf, 60, ' 3) '+mag_comp4+'       '+string(no_bins,format='(I3.3)')+'       # total magnitude   '  
        if comp4_type eq 'sersic' then begin
          printf, 60, ' 4) '+Re_comp4+'   1 band  #  R_e (half-light radius)   [pix]'
          printf, 60, ' 5) '+n_comp4+'             1 band  #  Sersic index n (de Vaucouleurs n=4) '
          printf, 60, ' 9) '+q_comp4+'        1 band  #  axis ratio (b/a)  '
          printf, 60, '10) '+pa_comp4+'   1 band  #  position angle (PA) [deg: Up=0, Left=90]'
        endif
        printf, 60, ' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
      endif
      
      if file_test(root+stars_file) eq 1 then begin
        readcol,root+stars_file,format='f,f,f',x_star,y_star,mag_star,comment='#',/SILENT
        for j=0,n_elements(x_star)-1,1 do begin
          x_pos=string(x_star[j],format='(F07.2)')
          y_pos=string(y_star[j],format='(F07.2)')
          mag_pos=string(mag_star[j],format='(F05.2)')
          
          for n=1,no_bins-1,1 do begin
            x_pos+=','+string(x_star[j],format='(F07.2)')
            y_pos+=','+string(y_star[j],format='(F07.2)')
            mag_pos+=','+string(mag_star[j],format='(F05.2)')
          endfor
          
          printf, 60, ' # Object number:  '+string(j)
          printf, 60, ' 0) psf                 #  object type'
          printf, 60, ' 1) '+x_pos+'   '+string(no_bins)+' band  #  position x'
          printf, 60, ' 2) '+y_pos+'   '+string(no_bins)+' band  #  position y'
          printf, 60, ' 3) '+mag_pos+'   '+string(no_bins)+' band  #  Integrated magnitude'
          printf, 60, ' '
          printf, 60, ' '
          printf, 60, ' '
        endfor
      endif

      printf, 60, '================================================================================'
      
    
      close,60
      
    endelse

endif

;=======================================================================================
;=======================================================================================
;=======================================================================================

if keyword_set(slices) then begin
  nband=no_bins
;  res=read_sersic_results(output+binned_dir+'imgblock.fits', nband, bd=1)
  if n_comp eq 1000 then res=read_sersic_results_2comp(output+binned_dir+'imgblock.fits', nband, bd=0) $
  else if n_comp eq 1100 then res=read_sersic_results_2comp(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'psf' then res=read_sersic_results_3psf(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1101 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1001 and comp4_type eq 'psf' then res=read_sersic_results_3psf(output+binned_dir+'imgblock.fits', nband, bd=0) $
  else if n_comp eq 1001 and comp4_type eq 'sersic' then res=read_sersic_results_3sersic(output+binned_dir+'imgblock.fits', nband, bd=0) $
  
  else if n_comp eq 1010  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(output+binned_dir+'imgblock.fits', nband, bd=0) $
  else if n_comp eq 1010  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(output+binned_dir+'imgblock.fits', nband, bd=0) $
  else if n_comp eq 1110  and comp3_type eq 'psf' then res=read_sersic_results_2comp_p(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1110  and comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1111 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(output+binned_dir+'imgblock.fits', nband, bd=1) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'psf' then res=read_sersic_results_3psf_p(output+binned_dir+'imgblock.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'psf' and comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(output+binned_dir+'imgblock.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(output+binned_dir+'imgblock.fits', nband, bd=0) $
  else if n_comp eq 1011 and comp4_type eq 'sersic' and comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(output+binned_dir+'imgblock.fits', nband, bd=0) 
  
;  test= file_search(output+slices_dir+'imgblock*',COUNT=nfiles1)
;  if nfiles1 gt 0 then spawn, 'rm '+output+slices_dir+'imgblock_*'
  
  
  close,70
  openw,70,output+slices_dir+'run_galfitm.sh'
  
  n=0
  fits_read,output+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits',crap,h_first
;  h_first=headfits(output+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits')
  x_size=sxpar(h_first,'NAXIS1')
  y_size=sxpar(h_first,'NAXIS2')
  
  no_images=no_slices                     ;no of image slices in each set
  total_images=final_image-first_image+1  ;total no of image slices
  no_loops=long(total_images/no_images)   ;calculate number of sets of images (i.e. no of feedme files to create)


  x1=no_images-1
  n_poly=0


   
  ;use cheb polynomials to calculate parameter values at each wavelength
;  wavelength_binned = fltarr(no_bins)
;  for bin eq 0,no_bins-1,1 do begin
;      h_temp=headfits(output+binned_dir+'image_'+string(bin,format='(I4.4)')+'.fits')
;      wavelength_binned[bin]=sxpar(h_temp,'WAVELENG')
;  endfor
  fits_read,output+binned_dir+'image_'+string(0,format='(I4.4)')+'.fits',crap,h_temp
;  h_temp=headfits(output+binned_dir+'image_'+string(0,format='(I4.4)')+'.fits')
  wave1=sxpar(h_temp,'WAVELENG')
  fits_read,output+binned_dir+'image_'+string(no_bins-1,format='(I4.4)')+'.fits',crap,h_temp  
;  h_temp=headfits(output+binned_dir+'image_'+string(no_bins-1,format='(I4.4)')+'.fits')
  wave2=sxpar(h_temp,'WAVELENG')
  wavelength_slices = fltarr(total_images)
  for bin=first_image,final_image,1 do begin
      fits_read,output+slices_dir+'image_'+string(bin,format='(I4.4)')+'.fits',crap,h_temp
;      h_temp=headfits(output+slices_dir+'image_'+string(bin,format='(I4.4)')+'.fits')
      wavelength_slices[bin-first_image]=sxpar(h_temp,'WAVELENG')
  endfor
  
  ;re_BAND_everywhere = cheb_interpolate(re_galfit_cheb_values, wavelengths_from_binned_fit(including start and finish), ALL_wavelengths, ordmax)
  ;ordmax is the order of polynomial that you used (Galfit number, e.g. 2 in your case, I think)
  
  
  wavelength_binned=fltarr(no_bins)
  for m=0,no_bins-1,1 do begin
    fits_read,output+binned_dir+'image_'+string(m,format='(I4.4)')+'.fits',crap,h
;    h=headfits(output+binned_dir+'image_'+string(m,format='(I4.4)')+'.fits')
    wavelength_binned[m]=sxpar(h,'WAVELENG')
  endfor
  
   ;y = chebeval(x, p, interval=[0d,10d])
  x_disk_all=chebeval(wavelength_slices,res.X_GALFIT_CHEB_D,INTERVAL=[wave1,wave2])
  y_disk_all=chebeval(wavelength_slices,res.Y_GALFIT_CHEB_D,INTERVAL=[wave1,wave2])
  
  mag_d_temp=res.MAG_GALFIT_BAND_D
  mag_d_temp[0]=mag_d_temp[1]
  mag_d_temp[n_elements(mag_d_temp)-1]=mag_d_temp[n_elements(mag_d_temp)-2]
  
  if disk_mag_polynomial_in le no_images/2 and disk_mag_polynomial_in ge 0 then mag_disk_all=chebeval(wavelength_slices,res.MAG_GALFIT_CHEB_D,INTERVAL=[wave1,wave2]) $
    else mag_disk_all=linear_interpolate(wavelength_slices,wavelength_binned,mag_d_temp)
  if disk_re_polynomial_in le no_images/2 and disk_re_polynomial_in ge 0 then Re_disk_all=chebeval(wavelength_slices,res.RE_GALFIT_CHEB_D,INTERVAL=[wave1,wave2]) $
    else Re_disk_all=linear_interpolate(wavelength_slices,wavelength_binned,res.RE_GALFIT_BAND_D)
  ;n_disk_all=cheb_interpolate(res.X_GALFIT_CHEB_D,wavelength_binned,wavelength_slices,0)
  q_disk_all=chebeval(wavelength_slices,res.Q_GALFIT_CHEB_D,INTERVAL=[wave1,wave2])
  pa_disk_all=chebeval(wavelength_slices,res.PA_GALFIT_CHEB_D,INTERVAL=[wave1,wave2])
  ;note= cheb_interpolate doesn't work for 0th order!
  n_disk_all=fltarr(total_images)
  n_disk_all[*]=res.N_GALFIT_BAND_D[0]

  if n_comp ge 1100 then begin
    x_bulge_all=chebeval(wavelength_slices,res.X_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
    y_bulge_all=chebeval(wavelength_slices,res.Y_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
;    mag_bulge_all=chebeval(wavelength_slices,res.MAG_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
;    Re_bulge_all=chebeval(wavelength_slices,res.RE_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])

    mag_b_temp=res.MAG_GALFIT_BAND_B
    mag_b_temp[0]=mag_b_temp[1]
    mag_b_temp[n_elements(mag_b_temp)-1]=mag_b_temp[n_elements(mag_b_temp)-2]

    if bulge_mag_polynomial_in le no_images/2 and bulge_mag_polynomial_in ge 0 then mag_bulge_all=chebeval(wavelength_slices,res.MAG_GALFIT_CHEB_B,INTERVAL=[wave1,wave2]) $
      else mag_bulge_all=linear_interpolate(wavelength_slices,wavelength_binned,mag_b_temp)
    if bulge_re_polynomial_in le no_images/2 and bulge_re_polynomial_in ge 0 then Re_bulge_all=chebeval(wavelength_slices,res.RE_GALFIT_CHEB_B,INTERVAL=[wave1,wave2]) $
      else Re_bulge_all=linear_interpolate(wavelength_slices,wavelength_binned,res.RE_GALFIT_BAND_B)
    n_bulge_all=chebeval(wavelength_slices,res.N_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
    q_bulge_all=chebeval(wavelength_slices,res.Q_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
    pa_bulge_all=chebeval(wavelength_slices,res.PA_GALFIT_CHEB_B,INTERVAL=[wave1,wave2])
    ;note= cheb_interpolate doesn't work for 0th order!
  ;  n_bulge_all=fltarr(total_images)
  ;  n_bulge_all[*]=res.X_GALFIT_BAND_D[0]
  endif
  if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
    x_comp3_all=chebeval(wavelength_slices,res.x_galfit_cheb_comp3,INTERVAL=[wave1,wave2])
    y_comp3_all=chebeval(wavelength_slices,res.y_galfit_cheb_comp3,INTERVAL=[wave1,wave2])

    mag_comp3_temp=res.MAG_GALFIT_BAND_COMP3
    mag_comp3_temp[0]=mag_comp3_temp[1]
    mag_comp3_temp[n_elements(mag_comp3_temp)-1]=mag_comp3_temp[n_elements(mag_comp3_temp)-2]

    mag_comp3_all=linear_interpolate(wavelength_slices,wavelength_binned,mag_comp3_temp)
    if comp3_type eq 'sersic' then begin
      Re_comp3_all=chebeval(wavelength_slices,res.re_galfit_cheb_comp3,INTERVAL=[wave1,wave2])
      n_comp3_all=chebeval(wavelength_slices,res.n_galfit_cheb_comp3,INTERVAL=[wave1,wave2])
      q_comp3_all=chebeval(wavelength_slices,res.q_galfit_cheb_comp3,INTERVAL=[wave1,wave2])
      pa_comp3_all=chebeval(wavelength_slices,res.pa_galfit_cheb_comp3,INTERVAL=[wave1,wave2])   
    endif   
  endif
  
  if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
    x_comp4_all=chebeval(wavelength_slices,res.x_galfit_cheb_comp4,INTERVAL=[wave1,wave2])
    y_comp4_all=chebeval(wavelength_slices,res.y_galfit_cheb_comp4,INTERVAL=[wave1,wave2])

    mag_comp4_temp=res.MAG_GALFIT_BAND_COMP4
    mag_comp4_temp[0]=mag_comp4_temp[1]
    mag_comp4_temp[n_elements(mag_comp4_temp)-1]=mag_comp4_temp[n_elements(mag_comp4_temp)-2]

    mag_comp4_all=linear_interpolate(wavelength_slices,wavelength_binned,mag_comp4_temp)
    if comp4_type eq 'sersic' then begin
      Re_comp4_all=chebeval(wavelength_slices,res.re_galfit_cheb_comp4,INTERVAL=[wave1,wave2])
      n_comp4_all=chebeval(wavelength_slices,res.n_galfit_cheb_comp4,INTERVAL=[wave1,wave2])
      q_comp4_all=chebeval(wavelength_slices,res.q_galfit_cheb_comp4,INTERVAL=[wave1,wave2])
      pa_comp4_all=chebeval(wavelength_slices,res.pa_galfit_cheb_comp4,INTERVAL=[wave1,wave2])  
    endif    
  endif
  
  ;NOTE- for this loop, running a total of no_loops+1. This is because 
  ;no_loops is an integer, giving the total number of loops with the full
  ;50 image slices. The final loop will catch all the remaining image slices
  for loop=0,no_loops,1 do begin
    ;insert details for first image slice used to obtain polynomials
    temp=abs(ugriz-(wavelength_slices[loop*no_images]))
    m=where(temp eq min(temp))
    psf_temp=PSF_files[m]
    
    n=0
    x0=first_image+(loop*no_images)
    fits_read,output+slices_dir+'image_'+string(x0,format='(I4.4)')+'.fits',crap,h
;    h=headfits(output+slices_dir+'image_'+string(x0,format='(I4.4)')+'.fits')
    file='image_'+string(x0,format='(I4.4)')+'.fits'
    band=string(n,format='(I3.3)')
    ;wavelength=string(sxpar(h,'WAVELENG'),format='(F08.3)')
;    psf='psf/'+string(x0,format='(I4.4)')+'.fits'
    psf='PSF/'+string(x0,format='(I4.4)')+'.fits'
    temp1=file_search(output+slices_dir+'badpix/badpix*.fits',COUNT=nfiles_bp)
    if nfiles_bp gt 2 then badpix='badpix/badpix_'+string(x0,format='(I4.4)')+'.fits' else badpix='badpix.fits'
    
    temp2=file_search(output+slices_dir+'sigma/sigma*.fits',COUNT=nfiles_sig)
    if nfiles_sig gt 0 then sigma='sigma/sigma_'+string(x0,format='(I4.4)')+'.fits' else sigma='none'
    
    magzpt=string(magzpt_in,format='(F04.1)')
    ;sky=string(res.SKY_GALFIT_BAND[loop],format='(F08.4)')
    resistant_mean, res.SKY_GALFIT_BAND,3,mean_sky
    sky=string(mean_sky,format='(F010.0)')
    sky_grad='0.0'
;    mag_B=string(res.MAG_GALFIT_BAND_B[loop],format='(F05.2)')
;    mag_D=string(res.MAG_GALFIT_BAND_D[loop],format='(F05.2)')
    ;mag_psf='25.0'    ;string(res.psf_galfit_band,format='(F05.2)')
    
      
    wavelength=string(wavelength_slices[loop*no_images],format='(F09.3)')
    x_D=string(x_disk_all[x0-first_image],format='(F07.2)')
    y_D=string(y_disk_all[x0-first_image],format='(F07.2)')
    mag_D=string(mag_disk_all[x0-first_image],format='(F05.2)')
    Re_D=string(Re_disk_all[x0-first_image],format='(F07.2)')
    n_D=string(n_disk_all[x0-first_image],format='(F06.3)')
    q_D=string(q_disk_all[x0-first_image],format='(F04.2)')
    pa_D=string(pa_disk_all[x0-first_image],format='(F07.2)')

    if n_comp ge 1100 then begin
      x_B=string(x_bulge_all[x0-first_image],format='(F07.2)')
      y_B=string(y_bulge_all[x0-first_image],format='(F07.2)')
      mag_B=string(mag_bulge_all[x0-first_image],format='(F05.2)')
      Re_B=string(Re_bulge_all[x0-first_image],format='(F07.2)')
      n_B=string(n_bulge_all[x0-first_image],format='(F06.3)')
      q_B=string(q_bulge_all[x0-first_image],format='(F04.2)')
      pa_B=string(pa_bulge_all[x0-first_image],format='(F07.2)')
    endif  
    
    if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
      x_comp3=string(x_comp3_all[x0-first_image],format='(F07.2)')
      y_comp3=string(y_comp3_all[x0-first_image],format='(F07.2)')
      mag_comp3=string(median(mag_comp3_all),format='(F05.2)');string(mag_comp3_all[x0-first_image],format='(F05.2)')
      if comp3_type eq 'sersic' then begin
        Re_comp3=string(re_comp3_all[x0-first_image],format='(F07.3)')
        n_comp3=string(n_comp3_all[x0-first_image],format='(F06.3)')
        q_comp3=string(q_comp3_all[x0-first_image],format='(F04.2)')
        pa_comp3=string(pa_comp3_all[x0-first_image],format='(F07.2)')
      endif
    endif
    if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
      x_comp4=string(x_comp4_all[x0-first_image],format='(F07.2)')
      y_comp4=string(y_comp4_all[x0-first_image],format='(F07.2)')
      mag_comp4=string(median(mag_comp4_all),format='(F05.2)');string(mag_comp4_all[x0-first_image],format='(F05.2)')
      if comp4_type eq 'sersic' then begin
        Re_comp4=string(re_comp4_all[x0-first_image],format='(F07.3)')
        n_comp4=string(n_comp4_all[x0-first_image],format='(F06.3)')
        q_comp4=string(q_comp4_all[x0-first_image],format='(F04.2)')
        pa_comp4=string(pa_comp4_all[x0-first_image],format='(F07.2)')
      endif
    endif
    ;use band values from binned images to calculate the values  
    ;for individual image slices by linear interpolation
    ;insert details for batches of 50 images within the range
    if loop ne no_loops then x1=no_images-1
    if loop eq no_loops then x1=(total_images mod no_images)-1    
    ;mod command => find remainder- this will be the number of images for the final galfitm file
    
    if disk_mag_polynomial_in lt 0 then disk_mag_polynomial=x1+1 else disk_mag_polynomial=disk_mag_polynomial_in
    if bulge_mag_polynomial_in lt 0 then bulge_mag_polynomial=x1+1 else bulge_mag_polynomial=bulge_mag_polynomial_in
    if comp3_mag_polynomial_in lt 0 and comp3_type eq 'sersic' then comp3_mag_polynomial=x1+1 else comp3_mag_polynomial=comp3_mag_polynomial_in

    if disk_re_polynomial_in lt 0 then disk_re_polynomial=x1+1 else disk_re_polynomial=0;disk_re_polynomial_in
    if disk_n_polynomial_in lt 0 then disk_n_polynomial=x1+1 else disk_n_polynomial=0;disk_n_polynomial_in
    if bulge_re_polynomial_in lt 0 then bulge_re_polynomial=x1+1 else bulge_re_polynomial=0;bulge_re_polynomial_in
    if bulge_n_polynomial_in lt 0 then bulge_n_polynomial=x1+1  else bulge_n_polynomial=0;bulge_n_polynomial_in
;    if disk_re_polynomial_in lt 0 and mean(Re_bulge_all) ge 3 then disk_re_polynomial=x1+1 else disk_re_polynomial=0
;    if disk_n_polynomial_in lt 0 and mean(n_bulge_all) ge 3 then disk_n_polynomial=x1+1 else disk_n_polynomial=0
;    if bulge_re_polynomial_in lt 0 and mean(Re_disk_all) ge 3 then bulge_re_polynomial=x1+1 else bulge_re_polynomial=0
;    if bulge_n_polynomial_in lt 0 and mean(n_disk_all) ge 3 then bulge_n_polynomial=x1+1 else bulge_n_polynomial=0
    
    if comp3_re_polynomial_in lt 0 and comp3_type eq 'sersic' then comp3_re_polynomial=x1+1 else comp3_re_polynomial=0
    if comp3_n_polynomial_in lt 0 and comp3_type eq 'sersic' then comp3_n_polynomial=x1+1 else comp3_n_polynomial=0
  ;  if comp4_re_polynomial lt 0 and comp4_type eq 'sersic' then comp4_re_polynomial=x1+1
  ;  if comp4_mag_polynomial lt 0 and comp4_type eq 'sersic' then comp4_mag_polynomial=x1+1
  ;  if comp4_n_polynomial lt 0 and comp4_type eq 'sersic' then comp4_n_polynomial=x1+1

    for n=1,x1,1 do begin
      x2=x0+n
      temp=abs(ugriz-(wavelength_slices[x0-first_image+n]))
      m=where(temp eq min(temp))
      psf_temp=PSF_files[m]
      
      fits_read,output+slices_dir+'image_'+string(x2,format='(I4.4)')+'.fits',crap,h
;      h=headfits(output+slices_dir+'image_'+string(x2,format='(I4.4)')+'.fits')
      file+=',image_'+string(x2,format='(I4.4)')+'.fits'
      band+=','+string(n,format='(I3.3)')
      ;wavelength=wavelength+','+string(sxpar(h,'WAVELENG'),format='(F08.3)')
      psf+=',PSF/'+string(x2,format='(I4.4)')+'.fits'
;      psf+=',psf/'+string(x2,format='(I4.4)')+'.fits'
      
      if n ne no_bins-1 and nfiles_bp gt 1 then badpix=badpix+',badpix/badpix_'+string(x2,format='(I4.4)')+'.fits' $
      else if n ne no_bins-1 and nfiles_bp le 1 then badpix=badpix+',badpix.fits' $
      else badpix=badpix+',badpix/badpix_end.fits'
      
      if nfiles_sig gt 0 then sigma=sigma+',sigma/sigma_'+string(x2,format='(I4.4)')+'.fits'
      
      ;badpix+=',badpix.fits'
      magzpt+=','+string(magzpt_in,format='(F04.1)')
;      sky=sky+','+string(res.SKY_GALFIT_BAND[loop],format='(F08.4)')
      sky+=','+string(mean_sky,format='(F010.0)')
      sky_grad+=',0.0'
;      mag_B+=','+string(res.MAG_GALFIT_BAND_B[loop],format='(F05.2)')
;      mag_D+=','+string(res.MAG_GALFIT_BAND_D[loop],format='(F05.2)')
      ;mag_psf=mag_psf+',25.0';','+string(res.psf_galfit_band,format='(F05.2)')
      
      wavelength+=','+string(wavelength_slices[x0-first_image+n],format='(F09.3)')
      x_D+=','+string(x_disk_all[x0-first_image+n],format='(F07.2)')
      y_D+=','+string(y_disk_all[x0-first_image+n],format='(F07.2)')
      mag_D+=','+string(mag_disk_all[x0-first_image+n],format='(F05.2)')
      Re_D+=','+string(Re_disk_all[x0-first_image+n],format='(F07.2)')
      n_D+=','+string(n_disk_all[x0-first_image+n],format='(F06.3)')
      q_D+=','+string(q_disk_all[x0-first_image+n],format='(F04.2)')
      pa_D+=','+string(pa_disk_all[x0-first_image+n],format='(F07.2)')
;print,loop,n,mag_disk_all[x0-first_image+n],Re_disk_all[x0-first_image+n]
      if n_comp ge 1100 then begin
        x_B+=','+string(x_bulge_all[x0-first_image+n],format='(F07.2)')
        y_B+=','+string(y_bulge_all[x0-first_image+n],format='(F07.2)')
        mag_B+=','+string(mag_bulge_all[x0-first_image+n],format='(F05.2)')
        Re_B+=','+string(Re_bulge_all[x0-first_image+n],format='(F07.2)')
        n_B+=','+string(n_bulge_all[x0-first_image+n],format='(F06.3)')
        q_B+=','+string(q_bulge_all[x0-first_image+n],format='(F04.2)')
        pa_B+=','+string(pa_bulge_all[x0-first_image+n],format='(F07.2)')
      endif
      if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
        x_comp3+=','+string(x_comp3_all[x0-first_image+n],format='(F07.2)')
        y_comp3+=','+string(y_comp3_all[x0-first_image+n],format='(F07.2)')
        mag_comp3+=','+string(median(mag_comp3_all),format='(F05.2)');string(mag_comp3_all[x0-first_image+n],format='(F05.2)')
        if comp3_type eq 'sersic' then begin
          Re_comp3+=','+string(re_comp3_all[x0-first_image+n],format='(F07.3)')
          n_comp3+=','+string(n_comp3_all[x0-first_image+n],format='(F06.3)')
          q_comp3+=','+string(q_comp3_all[x0-first_image+n],format='(F04.2)')
          pa_comp3+=','+string(pa_comp3_all[x0-first_image+n],format='(F07.2)')
        endif
      endif
      
      if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
        x_comp4+=','+string(x_comp4_all[x0-first_image+n],format='(F07.2)')
        y_comp4+=','+string(y_comp4_all[x0-first_image+n],format='(F07.2)')
        mag_comp4+=','+string(median(mag_comp4_all),format='(F05.2)');string(mag_comp4_all[x0-first_image+n],format='(F05.2)')
        if comp4_type eq 'sersic' then begin
          Re_comp4+=','+string(re_comp4_all[x0-first_image+n],format='(F07.3)')
          n_comp4+=','+string(n_comp4_all[x0-first_image+n],format='(F06.3)')
          q_comp4+=','+string(q_comp4_all[x0-first_image+n],format='(F04.2)')
          pa_comp4+=','+string(pa_comp4_all[x0-first_image+n],format='(F07.2)')
        endif
      endif

    endfor   
    
 
    close,60
    openw,60,output+slices_dir+'galfitm_'+string(loop,format='(I4.4)')+'.feedme'
    
    printf,60,'==============================================================================='
    printf,60,'# IMAGE and GALFIT CONTROL PARAMETERS';+$
    printf, 60, 'A) '+file+'             # Input data image (FITS file)'
    printf, 60, 'A1) '+band+'             # Band labels (can be omitted if fitting a single band)'
    printf, 60, 'A2) '+wavelength+'             # Band wavelengths'
    printf, 60, 'B) imgblock_'+string(loop,format='(I4.4)')+'.fits       # Output data image block
    printf, 60, 'C) '+sigma+'                # Sigma image name (made from data if blank or "none") 
    printf, 60, 'D) '+psf+'           # Input PSF image and (optional) diffusion kernel
    printf, 60, 'E) 1                   # PSF fine sampling factor relative to data 
    printf, 60, 'F) '+badpix+'                # Bad pixel mask (FITS image or ASCII coord list)
    printf, 60, 'G) galfitm.constraints                # File with parameter constraints (ASCII file)' 
    printf, 60, 'H) 1    '+string(x_size,format='(I4.4)')+'   1  '+string(y_size,format='(I4.4)')+'    # Image region to fit (xmin xmax ymin ymax)'
    printf, 60, 'I) '+string(x_size,format='(I4.4)')+'    '+string(x_size,format='(I4.4)')+'          # Size of the convolution box (x y)'
    printf, 60, 'J) '+magzpt+'              # Magnitude photometric zeropoint '
    printf, 60, 'K) 1  1        # Plate scale (dx dy)    [arcsec per pixel]'
    printf, 60, 'O) regular             # Display type (regular, curses, both)'
    printf, 60, 'P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps'
   
    printf, 60, ' '
    printf, 60, ' '
    ; 
;    x1=60
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
    printf, 60, '  1) '+sky+'   '+string(x1+1)+' band #  sky background at center of fitting region [ADUs]'
    printf, 60, '  2) '+sky_grad+'      0 band  #  dsky/dx (sky gradient in x)'
    printf, 60, '  3) '+sky_grad+'      0 band  #  dsky/dy (sky gradient in y)'
    printf, 60, '  Z) 0                      #  output option (0 = resid., 1 = Dont subtract) '
     
   
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
     
    printf, 60, '# Object number: 2'    ;disc
    printf, 60, ' 0) '+disk_type+'                 #  object type'
    printf, 60, ' 1) '+x_D+'   0 band  #  position x, y'
    printf, 60, ' 2) '+y_D+'   0 band  #  position x, y'
    printf, 60, ' 3) '+mag_D+'        '+string(disk_mag_polynomial,format='(I3.3)')+' band  #  Integrated magnitude' 
;    if disk_re_polynomial eq no_images then disk_re_polynomial=x1 else disk_re_polynomial=0
;    if disk_n_polynomial eq no_images then disk_n_polynomial=x1 else disk_n_polynomial=0
disk_n_polynomial=0
    printf, 60, ' 4) '+Re_D+'   '+string(disk_re_polynomial,format='(I3.3)')+' band  #  R_e (half-light radius)   [pix]'
    printf, 60, ' 5) '+n_D+'             '+string(disk_n_polynomial,format='(I3.3)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
    printf, 60, ' 9) '+Q_D+'        0 band  #  axis ratio (b/a)  '
    printf, 60, '10) '+PA_D+'   0 band  #  position angle (PA) [deg: Up=0, Left=90]'
    printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
    ;printf, 60, '#C0) 0.1         1      # traditional diskyness(-)/boxyness(+)'
   
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
    printf, 60, ' '
   
    if n_comp ge 1100 then begin
      printf, 60, ' # Object number: 3   '    ;bulge
      printf, 60, ' 0) '+bulge_type+'                 #  object type'
      printf, 60, ' 1) '+x_B+'   0 band  #  position x, y'
      printf, 60, ' 2) '+y_B+'   0 band  #  position x, y'
      printf, 60, ' 3) '+mag_B+'        '+string(bulge_mag_polynomial,format='(I3.3)')+' band  #  Integrated magnitude' 
;      if bulge_re_polynomial eq no_images then bulge_re_polynomial=x1 else bulge_re_polynomial=0
;      if bulge_n_polynomial eq no_images then bulge_n_polynomial=x1 else bulge_n_polynomial=0
      printf, 60, ' 4) '+Re_B+'   '+string(bulge_re_polynomial,format='(I3.3)')+' band  #  R_e (half-light radius)   [pix]'
      printf, 60, ' 5) '+n_B+'             '+string(bulge_n_polynomial,format='(I3.3)')+' band  #  Sersic index n (de Vaucouleurs n=4) '
      printf, 60, ' 9) '+q_B+'        0 band  #  axis ratio (b/a)  '
      printf, 60, '10) '+pa_B+'   0 band  #  position angle (PA) [deg: Up=0, Left=90]'
      printf, 60, ' Z) 0                      #  output option (0 = resid., 1 = Dont subtract)' 
      ;printf, 60, 'C0) 1.17         0      # traditional diskyness(-)/boxyness(+)'
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
    endif  
      
    if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
      printf, 60, ' 0) '+comp3_type+'                # object type'
      printf, 60, ' 1) '+x_comp3+'   0   #  position x, y'
      printf, 60, ' 2) '+y_comp3+'   0   #  position x, y'
      printf, 60, ' 3) '+mag_comp3+'       '+string(x1+1)+'       # total magnitude   '  
      if comp3_type eq 'sersic' then begin
        printf, 60, ' 4) '+Re_comp3+'   0 band  #  R_e (half-light radius)   [pix]'
        printf, 60, ' 5) '+n_comp3+'             0 band  #  Sersic index n (de Vaucouleurs n=4) '
        printf, 60, ' 9) '+q_comp3+'        0 band  #  axis ratio (b/a)  '
        printf, 60, '10) '+pa_comp3+'   0 band  #  position angle (PA) [deg: Up=0, Left=90]'
      endif
      printf, 60, ' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
    endif    
    if n_comp eq 1001 or n_comp eq 1101 or n_comp eq 1111 or n_comp eq 1011 then begin
      printf, 60, ' 0) '+comp4_type+'                # object type'
      printf, 60, ' 1) '+x_comp4+'   0   #  position x, y'
      printf, 60, ' 2) '+y_comp4+'   0   #  position x, y'
      printf, 60, ' 3) '+mag_comp4+'       '+string(x1+1)+'       # total magnitude   '  
      if comp4_type eq 'sersic' then begin
        printf, 60, ' 4) '+Re_comp4+'   0 band  #  R_e (half-light radius)   [pix]'
        printf, 60, ' 5) '+n_comp4+'             0 band  #  Sersic index n (de Vaucouleurs n=4) '
        printf, 60, ' 9) '+q_comp4+'        0 band  #  axis ratio (b/a)  '
        printf, 60, '10) '+pa_comp4+'   0 band  #  position angle (PA) [deg: Up=0, Left=90]'
      endif
      printf, 60, ' Z) 0                  #  Skip this model in output image?  (yes=1, no=0)'
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
      printf, 60, ' '
    endif
    
    if file_test(root+stars_file) eq 1 then begin
      readcol,root+stars_file,format='f,f,f',x_star,y_star,mag_star,comment='#',/SILENT
      for j=0,n_elements(x_star)-1,1 do begin
        x_pos=string(x_star[j],format='(F07.2)')
        y_pos=string(y_star[j],format='(F07.2)')
        mag_pos=string(mag_star[j],format='(F05.2)')
        
        for n=1,no_bins-1,1 do begin
          x_pos+=','+string(x_star[j],format='(F07.2)')
          y_pos+=','+string(y_star[j],format='(F07.2)')
          mag_pos+=','+string(mag_star[j],format='(F05.2)')
        endfor
        
        printf, 60, ' # Object number:  '+string(j)
        printf, 60, ' 0) psf                 #  object type'
          printf, 60, ' 1) '+x_pos+'   '+string(no_bins)+' band  #  position x'
          printf, 60, ' 2) '+y_pos+'   '+string(no_bins)+' band  #  position y'
          printf, 60, ' 3) '+mag_pos+'   '+string(no_bins)+' band  #  Integrated magnitude'
        printf, 60, ' '
        printf, 60, ' '
        printf, 60, ' '
      endfor
    endif    
    printf, 60, '================================================================================'
    
  
    close,60
    
    printf,70,galfitm+' galfitm_'+string(loop,format='(I4.4)')+'.feedme'
    printf,70,'mv imgblock_'+string(loop,format='(I4.4)')+'.fits imgblock_'+string(loop,format='(I4.4)')+'_fit.fits'
    
  endfor


  close,70


endif

delvarx, res
delvarx,sky_grad
delvarx,x_D
delvarx,y_D
delvarx,mag_D
delvarx,Re_D
delvarx,n_D
delvarx,PA_D
delvarx,q_D
delvarx,x_B
delvarx,y_B
delvarx,mag_B
delvarx,Re_B
delvarx,n_B
delvarx,PA_B
delvarx,q_B
delvarx,x_comp3
delvarx,y_comp3
delvarx,mag_comp3
delvarx,Re_comp3
delvarx,n_comp3
delvarx,PA_comp3
delvarx,q_comp3
delvarx,x_pos
delvarx,y_pos
delvarx,mag_pos


end

