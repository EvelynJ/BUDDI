@/home/ejohnsto/GitHub/BUDDI-GC/sex2ds9reg_buddi.pro
@/home/ejohnsto/GitHub/BUDDI-GC/add_1_component.pro


PRO detect_and_add, galfit_fits_output, new_output_filename, name_add,sextractor_executable,sextractor_setup,ds9_executable, nostop=nostop
; .run /home/ejohnsto/GitHub/BUDDI-GC/detect_and_add.pro
; detect_and_add, '/home/ejohnsto/NGC1407_GC/IFU_decomp/binned_images/imgblock.fits', 'Hallelujah.fits', 'eureka'
  zeropoint = 12
  ;sextractor_executable = '/usr/bin/sex'
  sextractor_parameters = '/raid/ejohnston/GitHub/BUDDI/buddi_gc_sex_params'
  ;sextractor_setup = '/home/ejohnsto/NGC1407_GC/sex_standard_setup.sex'
  ;ds9_executable = '/usr/local/bin/ds9'
  obj_standard_mag = 25.
  
; this code
;   runs SExtractor on a co-added image
;   add a new component for each detected object
;      should enable both SERSIC, GAUSSIAN and PSF (others? classified how?)


; derive and change folder
  spawn, 'pwd', infolder
  IF strpos(galfit_fits_output,'/') NE -1 THEN BEGIN
     workfolder = strmid(galfit_fits_output,0,strpos(galfit_fits_output,'/',/REVERSE_SEARCH))
     CD, workfolder
  ENDIF


; find all residual images and add them up to get an image for detection
; derive number of images/headers
  fits_info, galfit_fits_output, N_ext=layer_number, /silent
;  read galfit output file to derive header number, etc
;  rdfits_struct, galfit_fits_output, galfit_output_header,/header_only
;  rdfits_struct, galfit_fits_output, galfit_output_all
; identify residual images
  res_names = ' '
  FOR hn=0,layer_number DO BEGIN
     head=' '
     head = headfits(galfit_fits_output, exten = hn,/silent)
     IF strmid(sxpar(head, 'EXTNAME'),0,8) EQ 'RESIDUAL' THEN res_names = [res_names, sxpar(head, 'EXTNAME')]
  ENDFOR
  res_names = res_names[1:n_elements(res_names)-1]
;  print, 'found '+strtrim(n_elements(res_names),2)+' residuals'

; now read them into a structure
; read first image
  res_image = mrdfits(galfit_fits_output, res_names[0],/silent)
  res_image_str = fltarr(n_elements(res_image[*,0]),n_elements(res_image[0,*]),n_elements(res_names))
  FOR rn=0,n_elements(res_names)-1 DO BEGIN
     res_image_new = mrdfits(galfit_fits_output, res_names[rn],/silent)
     res_image_str[*,*,rn] = res_image_new 
     delvarx, res_image_new
  ENDFOR
  
; now add up images while ignoring first and last residual (fake images)
  n_add = 1
  res_image = res_image_str[*,*,1]
  FOR rn=2,n_elements(res_names)-2 DO BEGIN
     res_image += res_image_str[*,*,rn]
     n_add += 1
  ENDFOR
  writefits, 'average_residual.fits', res_image/n_add

  delvarx, res_image_str, res_image
  
; run SExtractor on this image
  sexcommand1 = sextractor_executable+' average_residual.fits -c '+sextractor_setup+ $
               ' -CATALOG_NAME detected_objects -CATALOG_TYPE ASCII' + $
               ' -PARAMETERS_NAME '+sextractor_parameters + $
               ' -MAG_ZEROPOINT '+strtrim(zeropoint,2)+ $
               ' -CHECKIMAGE_TYPE segmentation,apertures -CHECKIMAGE_NAME detected_objects_seg.fits,detected_objects_aper.fits'
  spawn, sexcommand1
  sexcommand2 = sextractor_executable+' average_residual.fits -c '+sextractor_setup+ $
               ' -CATALOG_NAME detected_objects.fits -CATALOG_TYPE FITS_1.0' + $
               ' -PARAMETERS_NAME '+sextractor_parameters + $
               ' -MAG_ZEROPOINT '+strtrim(zeropoint,2)+ $
               ' -CHECKIMAGE_TYPE segmentation,apertures -CHECKIMAGE_NAME detected_objects_seg.fits,detected_objects_aper.fits'
  spawn, sexcommand2  

; print object list as .reg file
  sex2ds9reg_buddi, 'detected_objects.fits','detected_objects.reg', 8

; open image and region file in ds9
  spawn, ds9_executable+' average_residual.fits -scale log -scale minmax -regions load detected_objects.reg &'

; find the latest galfit restart file
; find all matching galfit.??.band files and select newest
  spawn, 'ls '+strmid(galfit_fits_output,0,strpos(galfit_fits_output,'.fits'))+'.galfit.*', list

;; throw away all the files ending in 'band' or 'output'
  list = list[where(strmid(list,4,/reverse_offset) NE '.band')]
  list = list[where(strmid(list,6,/reverse_offset) NE '_output')]
  list = list[where(strmid(list,5,/reverse_offset,4) EQ 'fit.')]

  list2 = list
; isolate counting number
  FOR i=0,n_elements(list)-1 DO list2[i] = strmid(list[i],strpos(list[i],'.',/reverse_search)+1)

  list2 = fix(list2)   
; select latest file (file with highest number)
  wh = where(list2 EQ max(list2))
  obj = list[wh]
  new_name = obj+'_'+name_add
    
; print out instructions
  print, ' '
  print, ' A fairly sensitive SExtractor setup has been run. Please feel free to adapt this setting to detect the GCs'
  print, '  and run this (part of the) code again. The setup used is saved in '+sextractor_setup
  print, ' '
  print, 'Because detecting faint objects, especially in a non-ideal residual image, is terribly difficult'
  print, 'in a fully automated way, you NEED TO DO A FEW THINGS MANUALLY!'
  print, ' '
  print, 'In the ds9, that just opened (including regions to mark the detections), please to the following:'
  print, '- mark all false detections (you can mark several regions by clicking SHIFT while marking them)'
  print, '- remove those'
  print, '- make sure that all the detections you want to keep are nicely centered on the target.'
  print, '  This is important as in the fit, these positions will be constrained to wihtin 3 pixels!'
  print, '- you can also add missed detections by simply placing them over the missed objects.'
  print, '  For this you can either create new regions (best: circles) or copy existing ones.'
  print, ' '
  print, '- save the remaining GC detections (save regions) to a file called "detected_objects_psf.reg" in the same folder'
  print, '  This HAS to be done in the format "xy" amd coordinate system "image" '
  print, '  These objects will be inserted as PSFs into the fit (e.g. unresolved objects)'
  print, ' '
  print, '- You can do the same with objects that whould be inserted as SERSIC profiles.'
  print, '  These should be saved in a file calles detected_objects_sersic.reg'
  print, '  If these objects have been detected, the code will use the parameters found'
  print, '  In case they have not been detected (within 2 pixels), the fits will be starting at something small and round'
  print, '  If these do not converge, you might have to adapt their starting parameters in the start file.'
  print, ' '
  print, '- For bright objects (e.g. neigbouring galaxies, stars), please check the GALFIT start file that they are not already included!'
  print, ' '
  print, '- IN case you need to re-load the automatic regions, they are saved in "detected_objects.reg"'
  print, ' '
  print, 'All new objects will be added to '+strtrim(obj,2)
  print, ' and saved in '+new_name
  print, ' This file can then be restarted in GALFIT to fit these objects'
  print, ' '
  print, 'Once all this is done, type ".c"'
  print, ' '

; STOP to let the user to his stuff
  IF NOT keyword_set(nostop) THEN stop
  
  spawn, 'cp '+obj+' galfit_file_progress'
  det = mrdfits('detected_objects.fits',1,/silent)
  
; create new constraints file to be used
; go through setup file one line at a time and change line with output file format to include all information
  openr, filer, 'galfit_file_progress', /get_lun
  line = ''
  WHILE ~ EOF(filer) DO BEGIN
; Read a line of text:
     readf, filer, line
; get name of constraints file
     IF strpos(strtrim(line, 2), 'G) ') EQ 0 THEN BEGIN
        start = strtrim(strmid(line,0,strpos(line,')')+2),2)
        content = strtrim(strmid(line,strpos(line,')')+2, strpos(line,'#')-strpos(line,')')-2),2)
        comment = strtrim(strmid(line,strpos(line,'#')),2)
        old_constraints_file = content
     ENDIF
     print, line
  ENDWHILE
  close, filer
  FREE_LUN, filer
  new_constraints_name = old_constraints_file+'_'+name_add
  spawn, 'cp '+old_constraints_file+' '+new_constraints_name

; read in list for SERSIC
  IF file_test('detected_objects_sersic.reg') THEN BEGIN
     readcol, 'detected_objects_sersic.reg', xcen, ycen, format='F,F',/silent
     print, 'adding '+strtrim(n_elements(xcen),2)+' sersic objects'
     
; now add all new objects
; correlating to automatically detected objects'
     FOR ser=0,n_elements(xcen)-1 DO BEGIN
; getting parameters for detected objects
        srccor, det.x_image, det.y_image, xcen[ser], ycen[ser], 2, d, m, option = 1, /silent
        IF d[0] NE -1 THEN BEGIN
           xcengal = det[d].x_image
           ycengal = det[d].y_image
           maggal = det[d].mag_best
           regal = det[d].flux_radius
           argal = 1-det[d].ellipticity
           PAgal = det[d].theta_image-90.
        ENDIF ELSE BEGIN
; making up parameters for non-detected objects
           xcengal = xcen[ser]
           ycengal = ycen[ser]
           maggal = obj_standard_mag
           regal = 5.
           argal = 1.
           PAgal = 0.
        ENDELSE
        IF PAgal GT 180 THEN PAgal -= 180
        IF PAgal GT 90 THEN PAgal -= 180
        
; run add_1_component for all objects
        add_1_component, 'galfit_file_progress', 'ser'+strtrim(ser,2), new_output_filename, new_constraints_name, xcengal, ycengal, maggal, regal, 2.5, argal, PAgal, 'sersic'  
        spawn, 'mv '+'galfit_file_progress_ser'+strtrim(ser,2)+' galfit_file_progress'
     ENDFOR
  ENDIF
  
; read in list for PSF
  IF file_test('detected_objects_psf.reg') THEN BEGIN
     readcol, 'detected_objects_psf.reg', xcen, ycen, format='F,F',/silent
     print, 'adding '+strtrim(n_elements(xcen),2)+' PSF objects'

; correlating to automatically detected objects'
     FOR po=0,n_elements(xcen)-1 DO BEGIN
; getting parameters for detected objects
        srccor, det.x_image, det.y_image, xcen[po], ycen[po], 2, d, m, option = 1, /silent
        IF d[0] NE -1 THEN BEGIN
           xcengal = det[d].x_image
           ycengal = det[d].y_image
           maggal = det[d].mag_best
        ENDIF ELSE BEGIN
; making up parameters for non-detected objects
           xcengal = xcen[po]
           ycengal = ycen[po]
           maggal = obj_standard_mag
        ENDELSE
        
; run add_1_component for all objects
        add_1_component, 'galfit_file_progress', 'psf'+strtrim(po,2), new_output_filename, new_constraints_name, xcengal, ycengal, maggal, 0., 0., 0., 0., 'psf'  
        spawn, 'mv '+'galfit_file_progress_psf'+strtrim(po,2)+' galfit_file_progress'
     ENDFOR
  ENDIF
  spawn, 'mv galfit_file_progress '+new_name

  CD, infolder

END
