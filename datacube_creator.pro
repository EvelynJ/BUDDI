


pro datacube_creator,setup,info,wavelength,original_datacube,bestfit_datacube,$
  residual_datacube,disk_datacube,residual_sky_datacube,bulge_datacube,$
  comp3_datacube,comp4_datacube,MANGA=manga,CALIFA=califa,KEEP_CUBES=keep_cubes


  root=setup.root
  decomp=setup.decomp
  decomp_dir=setup.decomp_dir
  kinematics=setup.kinematics
  galaxy_ref=setup.galaxy_ref
  file=setup.file
  slices_dir=setup.slices_dir
  n_comp=setup.n_comp
  comp3_type=setup.comp3_type
  comp4_type=setup.comp4_type
  no_slices=setup.no_slices
  

first_image=info[0]
final_image=info[1]
no_bins=info[2]
images_per_bin=info[3]
start_wavelength=info[4]
end_wavelength=info[5]
;no_images=no_slices

total_images=final_image-first_image+1
x1=(total_images mod no_slices)     ;total number of images in the last feedme file


;fits_read,root+kinematics+file+'_counts.fits',temp_input,h
fits_read,root+decomp+galaxy_ref+'_smoothed_FLUX.fits',temp_input,h
side1=sxpar(h,'NAXIS1')
side2=sxpar(h,'NAXIS2')
images=sxpar(h,'NAXIS3')

result = FILE_TEST(root+decomp+decomp_dir, /DIRECTORY) 
if result eq 0 then file_mkdir,root+decomp+decomp_dir

bulge_datacube=fltarr(side1,side2,total_images)
disk_datacube=fltarr(side1,side2,total_images)
comp3_datacube=fltarr(side1,side2,total_images)
comp4_datacube=fltarr(side1,side2,total_images)
original_datacube=fltarr(side1,side2,total_images)
bestfit_datacube=fltarr(side1,side2,total_images)
residual_datacube=fltarr(side1,side2,total_images)    ;residuals after subtracting model from original
residual_sky_datacube=fltarr(side1,side2,total_images)  ;residual sky level fitted as component 1

galfit_or_galfitm='galfitm'



if galfit_or_galfitm eq 'galfitm' then begin
  
  
  j=0
  result_subcomps = file_search(root+decomp+slices_dir+'subcomps*.fits',COUNT=nfiles_subcomps)
  result_feedme = file_search(root+decomp+slices_dir+'galfitm_*.feedme',COUNT=nfiles_feedme)
;  result = file_search(root+decomp+slices_dir+'subcomps*.fits',COUNT=nfiles1)
  
  nfiles=nfiles_feedme
  for n=0,nfiles-1,1 do begin
      ;For bulge and disc, read in the fits file, the try to select every 
      ;third image (change when I include the PSF fit) after the first
      ;50 to go into the bulge or disc arrays. 
      tempy=file_search(root+decomp+slices_dir+'subcomps_'+string(n,format='(i4.4)')+'.fits',COUNT=nfiles_subcomps)
      result=file_test(root+decomp+slices_dir+'subcomps_'+string(n,format='(i4.4)')+'.fits')
      
      if result eq 1 then begin
        fits_open,root+decomp+slices_dir+'subcomps_'+string(n,format='(I4.4)')+'.fits',subcomps
        fits_open,root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'_fit.fits',imgblock
        if n ne nfiles-1 then no_images=no_slices else no_images=x1
        if no_images gt 1 then begin
          for m=0,no_images-1,1 do begin
            fits_read,subcomps,disk_in,header_in,EXTNAME='COMPONENT_2_sersic _'+string(m,format='(I3.3)')
            disk_datacube[*,*,j]=disk_in


            fits_read,subcomps,sky_in,header_in,EXTNAME='COMPONENT_1_sky _'+string(m,format='(I3.3)')
            residual_sky_datacube[*,*,j]=sky_in

            if n_comp ge 1100 then begin
              fits_read,subcomps,bulge_in,header_in,EXTNAME='COMPONENT_3_sersic _'+string(m,format='(I3.3)')
              bulge_datacube[*,*,j]=bulge_in 

            endif
          
          
            if  n_comp eq 1110 or n_comp eq 1111 then begin
              if comp3_type eq 'sersic' then begin
;                if n_comp eq 1010 or n_comp eq 1011 then fits_read,subcomps,comp3_in,header_in,EXTNAME='COMPONENT_3_sersic _'+string(m,format='(I3.3)')
                if n_comp eq 1110 or n_comp eq 1111 then fits_read,subcomps,comp3_in,header_in,EXTNAME='COMPONENT_4_sersic _'+string(m,format='(I3.3)')
              endif else begin
;                if n_comp eq 1010 or n_comp eq 1011 then fits_read,subcomps,comp3_in,header_in,EXTNAME='COMPONENT_3_psf _'+string(m,format='(I3.3)')
                if n_comp eq 1110 or n_comp eq 1111 then fits_read,subcomps,comp3_in,header_in,EXTNAME='COMPONENT_4_psf _'+string(m,format='(I3.3)')
              endelse
              comp3_datacube[*,*,j]=comp3_in

            endif
          
            if n_comp eq 1111 then begin
              if comp4_type eq 'sersic' then begin
;                if n_comp eq 1001 then fits_read,subcomps,comp4_in,header_in,EXTNAME='COMPONENT_3_sersic _'+string(m,format='(I3.3)')
;                if n_comp eq 1101 or n_comp eq 1011 then fits_read,subcomps,comp4_in,header_in,EXTNAME='COMPONENT_4_sersic _'+string(m,format='(I3.3)')
                if n_comp eq 1111 then fits_read,subcomps,comp4_in,header_in,EXTNAME='COMPONENT_5_sersic _'+string(m,format='(I3.3)')
              endif else begin
;                if n_comp eq 1001 then fits_read,subcomps,comp4_in,header_in,EXTNAME='COMPONENT_3_psf _'+string(m,format='(I3.3)')
;                if n_comp eq 1101 or n_comp eq 1011 then fits_read,subcomps,comp4_in,header_in,EXTNAME='COMPONENT_4_psf _'+string(m,format='(I3.3)')
                if n_comp eq 1111 then fits_read,subcomps,comp4_in,header_in,EXTNAME='COMPONENT_5_psf _'+string(m,format='(I3.3)')
              endelse  
              comp4_datacube[*,*,j]=comp4_in

            endif
              
              
               
            fits_read,imgblock,original_in,header_in,EXTNAME='INPUT_'+string(m,format='(I3.3)')
            fits_read,imgblock,bestfit_in,header_in,EXTNAME='MODEL_'+string(m,format='(I3.3)')
            fits_read,imgblock,residuals_in,header_in,EXTNAME='RESIDUAL_'+string(m,format='(I3.3)')
            
            original_datacube[*,*,j]=original_in;[*,*,a+(n*tot_images)-1:z+(n*tot_images)-1]
            bestfit_datacube[*,*,j]=bestfit_in;[*,*,a+(n*tot_images)-1:z+(n*tot_images)-1]
            residual_datacube[*,*,j]=residuals_in;[*,*,a+(n*tot_images)-1:z+(n*tot_images)-1]
            

            
            
            
;            if n ne 0 and m eq 0 then begin
;;              original_datacube[*,*,j-1]=0.5*(original_datacube[*,*,j-2]+original_datacube[*,*,j])
;              bestfit_datacube[*,*,j-1]=0.5*(bestfit_datacube[*,*,j-2]+bestfit_datacube[*,*,j])
;              residual_datacube[*,*,j-1]= 0.5*(residual_datacube[*,*,j-2]+residual_datacube[*,*,j])
;              disk_datacube[*,*,j-1]= 0.5*(disk_datacube[*,*,j-2]+disk_datacube[*,*,j])
;              if n_comp ge 1100 then bulge_datacube[*,*,j-1]= 0.5*(bulge_datacube[*,*,j-2]+bulge_datacube[*,*,j])
;              if n_comp eq 1110 or n_comp eq 1010 or n_comp eq 1111 then comp3_datacube[*,*,j-1]= 0.5*(comp3_datacube[*,*,j-2]+comp3_datacube[*,*,j])
;              if n_comp eq 1111 then comp4_datacube[*,*,j-1]= 0.5*(comp4_datacube[*,*,j-2]+comp4_datacube[*,*,j])
;            endif
            
            j+=1
          endfor
        endif
        fits_close,imgblock
        fits_close,subcomps
        
      endif else begin
        if n ne nfiles-1 then no_images=no_slices else no_images=x1
        
        for m=0,no_images-1,1 do begin
          disk_datacube[*,*,j]=-99
          bulge_datacube[*,*,j]=-99
          comp3_datacube[*,*,j]=-99
          comp4_datacube[*,*,j]=-99
          residual_sky_datacube[*,*,j]=-99
          
          original_datacube[*,*,j]=-99
          bestfit_datacube[*,*,j]=-99
          residual_datacube[*,*,j]=-99
          j+=1
        endfor
        
      endelse 
  endfor
  
  h_temp=headfits(root+file+'.fits')
  h_tempy = headfits(root+decomp+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits')
;  fits_read,root+decomp+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits',tempycrap,h
  wavelength0=sxpar(h_tempy,'WAVELENG') 
  step=setup.step;sxpar(h_temp,'CD3_3')               
  ;print,'***',wavelength0,step
  
  tempf=mrdfits(root+file+'.fits',1,h_flux)
  delvarx,tempf,h_tempy
  if keyword_set(manga) then begin
    ;read in header for original image. Use the CRVAL3 and CD3_3 parameters to work out the step size in log units
    sxaddpar,h_temp,'Wave0',(wavelength0)
    sxaddpar,h_temp,'CRVAL3',alog10(wavelength0)
    sxaddpar,h_temp,'CD3_3',step
    sxaddpar,h_flux,'Wave0',(wavelength0)
    sxaddpar,h_flux,'CRVAL3',alog10(wavelength0)
    sxaddpar,h_flux,'CD3_3',step
  endif else if keyword_set(califa) then begin
    sxaddpar,h_temp,'Wave0',alog(wavelength0)
    sxaddpar,h_temp,'CRVAL3',(wavelength0)
    sxaddpar,h_temp,'CD3_3',step
    sxaddpar,h_flux,'WAVE0',alog(wavelength0)
    sxaddpar,h_flux,'CRVAL3',(wavelength0)
    sxaddpar,h_flux,'CD3_3',step
  endif
  s=size(bestfit_datacube)
  sxaddpar,h_flux,'NAXIS3',s[3]
  
  if keyword_set(keep_cubes) then begin
    
    s=size(bestfit_datacube)
    for j=9,s[3]-5,10 do bestfit_datacube[*,*,j]=0.5*(bestfit_datacube[*,*,j-1]+bestfit_datacube[*,*,j+1])
    for j=9,s[3]-5,10 do residual_datacube[*,*,j]=0.5*(residual_datacube[*,*,j-1]+residual_datacube[*,*,j+1])
    s=size(disk_datacube)
    for j=9,s[3]-5,10 do disk_datacube[*,*,j]=0.5*(disk_datacube[*,*,j-1]+disk_datacube[*,*,j+1])
    for j=9,s[3]-5,10 do residual_sky_datacube[*,*,j]=0.5*(residual_sky_datacube[*,*,j-1]+residual_sky_datacube[*,*,j+1])

    fits_write,root+decomp+decomp_dir+'original_cube.fits',original_datacube,h_flux,extname='FLUX'
    modfits,root+decomp+decomp_dir+'original_cube.fits',0,h_temp
    ;modfits,root+decomp+decomp_dir+'original_cube.fits',1,h_flux,extname='FLUX'
    fits_write,root+decomp+decomp_dir+'bestfit_cube.fits',bestfit_datacube,h_flux,extname='FLUX'
    modfits,root+decomp+decomp_dir+'bestfit_cube.fits',0,h_temp
    ;modfits,root+decomp+decomp_dir+'bestfit_cube.fits',1,h_flux,extname='FLUX'
    fits_write,root+decomp+decomp_dir+'residuals_cube.fits',residual_datacube,h_flux,extname='FLUX'
    modfits,root+decomp+decomp_dir+'residuals_cube.fits',0,h_temp
    ;modfits,root+decomp+decomp_dir+'residuals_cube.fits',1,h_flux,extname='FLUX'
  
    fits_write,root+decomp+decomp_dir+'component1_cube.fits',disk_datacube,h_flux,extname='FLUX'
    modfits,root+decomp+decomp_dir+'component1_cube.fits',0,h_temp
    ;modfits,root+decomp+decomp_dir+'component1_cube.fits',1,h_flux,extname='FLUX'
    fits_write,root+decomp+decomp_dir+'residual_sky_cube.fits',residual_sky_datacube,h_flux,extname='FLUX'
    modfits,root+decomp+decomp_dir+'residual_sky_cube.fits',0,h_temp
    ;modfits,root+decomp+decomp_dir+'residual_sky_cube.fits',1,h_flux,extname='FLUX'
    if n_comp ge 1100 then begin
      for j=9,s[3]-5,10 do bulge_datacube[*,*,j]=0.5*(bulge_datacube[*,*,j-1]+bulge_datacube[*,*,j+1])
      fits_write,root+decomp+decomp_dir+'component2_cube.fits',bulge_datacube,h_flux,extname='FLUX'
      modfits,root+decomp+decomp_dir+'component2_cube.fits',0,h_temp
      ;modfits,root+decomp+decomp_dir+'component2_cube.fits',1,h_flux,extname='FLUX'
    endif
    if  n_comp eq 1110 or n_comp eq 1111 then begin
      for j=9,s[3]-5,10 do comp3_datacube[*,*,j]=0.5*(comp3_datacube[*,*,j-1]+comp3_datacube[*,*,j+1])
      fits_write,root+decomp+decomp_dir+'component3_cube.fits',comp3_datacube,h_flux,extname='FLUX'
      modfits,root+decomp+decomp_dir+'component3_cube.fits',0,h_temp
      ;modfits,root+decomp+decomp_dir+'component3_cube.fits',1,h_flux,extname='FLUX'
    endif
    if n_comp eq 1111 then begin
      for j=9,s[3]-5,10 do comp4_datacube[*,*,j]=0.5*(comp4_datacube[*,*,j-1]+comp4_datacube[*,*,j+1])
      fits_write,root+decomp+decomp_dir+'component4_cube.fits',comp4_datacube,h_flux,extname='FLUX'
      modfits,root+decomp+decomp_dir+'component4_cube.fits',0,h_temp
      ;modfits,root+decomp+decomp_dir+'component4_cube.fits',1,h_flux,extname='FLUX'
    endif
  endif else begin
    fits_write,root+decomp+decomp_dir+'component1_cube.fits',disk_datacube,h_flux,extname='FLUX'
    modfits,root+decomp+decomp_dir+'component1_cube.fits',0,h_temp
    ;modfits,root+decomp+decomp_dir+'component1_cube.fits',1,h_flux,extname='FLUX'
  endelse
  
  
  delvarx,original_in,bestfit_in,residuals_in,comp4_in
  delvarx,comp3_in,bulge_in,disk_in

;======================================================================
endif else if galfit_or_galfitm eq 'galfit' then begin
  
  
  j=0
;  result = file_search(root+decomp+slices_dir+'galfit*.feedme',COUNT=nfiles)
  result = file_search(root+decomp+slices_dir+'imgblock*_fit.fits',COUNT=nfiles)
  
  
  
  for n=0,nfiles-1,1 do begin
      ;For bulge and disc, read in the fits file, the try to select every 
      ;third image (change when I include the PSF fit) after the first
      ;50 to go into the bulge or disc arrays. 
      if n eq 0 then a=0 else a=1
      if n eq nfiles-1 then z=x1 else z=no_images
      
      fits_open,root+decomp+slices_dir+'imgblock_'+string(n,format='(I4.4)')+'.fits',imgblock
      print,'loop, run ',n
      
      for m=a,z,1 do begin
  ;    ;For best fit and residuals, need to read in batches of 50, after 
  ;    ;the first 50, which are blank.
  ;    
        fits_read,imgblock,original_in,header_in,EXTNAME='INPUT_'+string(m,format='(I3.3)')
        fits_read,imgblock,bestfit_in,header_in,EXTNAME='MODEL_'+string(m,format='(I3.3)')
        fits_read,imgblock,residuals_in,header_in,EXTNAME='RESIDUAL_'+string(m,format='(I3.3)')
        
        original_datacube[*,*,j]=original_in;[*,*,a+(n*tot_images)-1:z+(n*tot_images)-1]
        bestfit_datacube[*,*,j]=bestfit_in;[*,*,a+(n*tot_images)-1:z+(n*tot_images)-1]
        residual_datacube[*,*,j]=residuals_in;[*,*,a+(n*tot_images)-1:z+(n*tot_images)-1]
        j+=1
      endfor
      fits_close,imgblock
      
  endfor
  
  fits_Read,root+decomp+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits',tempy,h_temp
;  h_temp = headfits(root+decomp+slices_dir+'image_'+string(first_image,format='(I4.4)')+'.fits')
  wavelength0=sxpar(h_temp,'WAVELENG')
  sxaddpar,h,'CRVAL3',alog10(wavelength0)
  
  fits_write,root+decomp+decomp_dir+'original.fits',original_datacube,h
  fits_write,root+decomp+decomp_dir+'bestfit.fits',bestfit_datacube,h
  fits_write,root+decomp+decomp_dir+'residuals.fits',residual_datacube,h
  


  j=0
  result = file_search(root+decomp+slices_dir+'models/subcomps*.feedme',COUNT=nfiles)
    
    for n=0,nfiles-1,1 do begin
      ;For bulge and disc, read in the fits file, the try to select every 
      ;third image (change when I include the PSF fit) after the first
      ;50 to go into the bulge or disc arrays. 
      if n eq 0 then a=0 else a=1
      if n eq nfiles-1 then z=x1 else z=no_images
      
      fits_open,root+decomp+slices_dir+'models/subcomps_'+string(n+first_image,format='(I4.4)')+'.fits',sub_comps
      print,'loop2, run ',n
      
      
      fits_read,sub_comps,disk_in,header_in,EXTEN_NO=2;EXTNAME='COMPONENT_2_sersic _'+string(m,format='(I3.3)')
      disk_datacube[*,*,n]=disk_in
      if n_comp eq 110 or n_comp eq 111 then begin
        fits_read,sub_comps,bulge_in,header_in,EXTEN_NO=3;EXTNAME='COMPONENT_3_sersic _'+string(m,format='(I3.3)')
        bulge_datacube[*,*,n]=bulge_in
      endif
        
      fits_close,sub_comps
      
  endfor  
  
  
  fits_write,root+decomp+decomp_dir+'disk.fits',disk_datacube,h
  if n_comp eq 110 or n_comp eq 111 then fits_write,root+decomp+decomp_dir+'bulge.fits',bulge_datacube,h





endif else begin
  
  error_message,'datacube_creator.pro ==> Please choose Galfit or Galfitm'
  
endelse


end

