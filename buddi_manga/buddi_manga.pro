;
; Preliminary code to try to automate the running of BUDDI on all MaNGA galaxies.
; Version 0.0 :  February 2020, EJo, Santiago, Chile
; Version 1.0 :  May 2020, EJo & BHa. Added IDL bridge, then removed it as it didn't work
;
; To run in terminal:
;     buddi_manga,1
; Note: the number after calling buddi_manga is the session number. This helps the 
; user to keep track of when a fit fails if they are running the code on multiple 
; terminal windows. The session number appears in the Ǵalaxies_in_progress' file, 
; so when a session fails, one can check which galaxy was the last one used for 
; that session number
;

pro BUDDI_MaNGA,session,download_only=download_only
  print, systime(0)
;###########################################
;set up directories for fits and results
;###########################################
  
;####### USER INPUT ######
output_folder = '/home/kjegathe/BUDDI_out/DR17_fits'
static_files = '/home/kjegathe/BUDDI/static_files/'
stellar_lib = '/home/kjegathe/BUDDI/IDL_BUDDI/miles_models'         ;path to stellar library
galfit_exe = '/usr/bin/galfitm' ;path to galfit and version
how_many_gals = 0    ; define how many objects you want to look through (in catalogue, not actual galaxies) 
                                ; set == 0 to switch off this feature
                                ; and do ALL galaxies 
;####### END USER INPUT ######


  
; add / to every path and file! to avoid obvious and fatal errors
  output_folder = output_folder+'/'
  static_files = static_files+'/'
  stellar_lib = stellar_lib+'/'
  
; define working and final sub-directories
  root_fits=output_folder+'/fits/'
  root_output=output_folder+'/output/'
; make those files automatically
  IF NOT file_test(root_fits) THEN spawn, 'mkdir -p '+root_fits
  IF NOT file_test(root_output) THEN spawn, 'mkdir -p '+root_output
  
  CD, root_fits, CURRENT=start_dir
  
;###########################################
; * Read in a list of MaNGA galaxies
;###########################################
;; DR15
;  Manga_list=mrdfits(static_files+'dapall-v2_4_3-2.2.1.fits',1)
;  Pymorph_list=mrdfits(static_files+'manga-pymorph-DR15.fits',2) ;extensions 1,2,3 are g,r,i band results
;  Pymorph_SPA=mrdfits(static_files+'manga-pymorph-DR15-SPA.fits',1) 

; DR17
  Manga_list = mrdfits(static_files+'dapall-v3_1_1-3.1.0.fits',1)

;  WHICH LAYER????
  Pymorph_list = mrdfits(static_files+'manga-pymorph-DR17.fits',2) ;extensions 1,2,3 are g,r,i band results, DR17 NOW INCLUDES SPA!!
  gz3d_list = mrdfits(static_files+'gz3d_metadata.fits',1)
    
; cut off empty spaces in mangaIDs!!
  manga_list.mangaid = strtrim(manga_list.mangaid,2)
  gz3d_list.mangaid = strtrim(gz3d_list.mangaid,2)
  pymorph_list.manga_id = strtrim(pymorph_list.manga_id,2)

  totalgals = n_elements(manga_list.plate)
; select only objects on large IFUs
  manga_list = manga_list[where(manga_list.ifudesign GE 9000)]
  largegals = n_elements(manga_list.plate)

; select UNIQUE objects!
;  manga_list = manga_list[uniq(manga_list.mangaid,sort(manga_list.mangaid))]
  uniquegals = n_elements(manga_list.plate)
; this loses one of the pymorph fits, if it has been fit twice!!
  Pymorph_list = pymorph_list[uniq(pymorph_list.manga_id,sort(pymorph_list.manga_id))]  

;;; DR15  
;; SPA_PYM=Pymorph_SPA.SPA_R
;; DR17 (done directly below)
;  SPA_PYM = Pymorph_list.SPA

;###########################################
;  * Go to the next galaxy on the list
;###########################################
  
; make sure the files below exist, and create them if not
  exist=file_test(root_fits+'galaxies_in_progress.txt')
  spawn,'touch '+root_fits+'galaxies_in_progress.txt'
  spawn,'touch '+root_fits+'galaxies_completed_SS.txt'
  spawn,'touch '+root_fits+'galaxies_completed_SE.txt'
  spawn,'touch '+root_fits+'galaxies_failed_SS.txt'
  spawn,'touch '+root_fits+'galaxies_failed_SE.txt'
  spawn,'touch '+root_fits+'galaxies_not_in_gz3d.txt'
  spawn,'touch '+root_fits+'galaxies_not_in_pymorph.txt'
  
;if the files above did not exist, open the newly created files and add a line at the top to prevent readcol crashing later on
  if exist eq 0 then begin
     openw,1,root_fits+'galaxies_in_progress.txt'
     openw,2,root_fits+'galaxies_completed_SS.txt'
     openw,3,root_fits+'galaxies_completed_SE.txt'
     openw,4,root_fits+'galaxies_failed_SS.txt'
     openw,5,root_fits+'galaxies_failed_SE.txt'
     openw,6,root_fits+'galaxies_not_in_gz3d.txt'
     openw,7,root_fits+'galaxies_not_in_pymorph.txt'
     printf,1,'0000-00000 0'
     printf,2,'0000-00000 0'
     printf,3,'0000-00000 0'
     printf,4,'0000-00000 0'
     printf,5,'0000-00000 0'
     printf,6,'0000-00000 0'
     printf,7,'0000-00000 0'
     close,1
     close,2
     close,3
     close,4
     close,5
     close,6
     close,7
  endif
    
  cntpymorph = 0
  cntnotpymorph = 0
  cntpymorphssonly = 0
  cntpymorphseonly = 0  
  cntpymorphnone = 0  
  cntpymorphboth = 0
  cntnotingz3d = 0
  cntingz3d = 0
  cntbuddi_done_ss = 0
  cntbuddi_done_se = 0
  IF how_many_gals EQ 0 THEN how_many_gals = n_elements(Manga_list.plate)
; REMOVE LIMIT!!!!!
  FOR loop=0, how_many_gals-1 DO BEGIN 
     print, ' '
     print, strtrim(loop,2)+'    '
;     loop-=1 ;subtract a value so that the next line works without skipping a galaxy
;     repeat loop+=1 until Manga_list[loop].ifudesign ge 9000
     
     print,'***Now running for galaxy '+strtrim(loop,2)+' out of '+strtrim(n_elements(Manga_list.plate),2)
     
;;identify galaxy in the Pymorph catalog using Manga_list[loop].plateifu
;     element = where(strtrim(Pymorph_list.plateifu,2) eq strtrim(Manga_list[loop].plateifu,2), cntwh)
;identify galaxy in the Pymorph catalog using Manga_ID
     element = where(strtrim(Pymorph_list.manga_ID,2) eq strtrim(Manga_list[loop].mangaid,2), cntwh)
     print,'Manga IFU plate: '+strtrim(Manga_list[loop].plateifu,2)
     print,'Manga ID: '+strtrim(manga_list[loop].mangaid,2)
    
;read in the list of galaxies that have been tried to identify the next in the sequence   
     readcol,root_fits+'galaxies_in_progress.txt',format='a',in_prog,comment='#',/silent
     readcol,root_fits+'galaxies_completed_SS.txt',format='a',completed_SS,comment='#',/silent
     readcol,root_fits+'galaxies_completed_SE.txt',format='a',completed_SE,comment='#',/silent
     readcol,root_fits+'galaxies_not_in_pymorph.txt',format='a',not_in_pymorph,comment='#',/silent
     
; stop if more than one object was found (because it indicates an error)
     IF cntwh GT 1 THEN BEGIN
        print, 'more than one objects found in pymorph catalogue, choosing first entry!'
        element=element[0]
        stop
     ENDIF
     IF element EQ -1 THEN BEGIN
        cntnotpymorph += 1
        print,'galaxy is not in PYMORPH catalog. Moving on'
; avoid reporting each object mulitple times
        pymorph_test = where(not_in_pymorph eq Manga_list[loop].plateifu,cnt)
        IF cnt EQ 0  THEN BEGIN
           openw,1,root_fits+'galaxies_not_in_pymorph.txt',/APPEND
           printf,1,Manga_list[loop].plateifu,'    ',session
           close,1
        ENDIF
        CONTINUE
     ENDIF
     
     run_test = where(in_prog eq Manga_list[loop].plateifu,cnt)
     if cnt ge 1 then begin
;if code finds that galaxy has been started
        print,'galaxy is running in another window, or has already been done. Moving on'
        continue
     endif else begin
;if code finds that galaxy has not been started
        print,'galaxy is untouched, starting loop'
        
;adding the next galaxy to the in progress file to prevent another session repeating that fit
        IF NOT keyword_set(download_only) THEN BEGIN
           openw,1,root_fits+'galaxies_in_progress.txt',/APPEND
           printf,1,Manga_list[loop].plateifu,'    ',session
           close,1
        ENDIF

        print, systime(0)
        
;if the fit was successful with PyMorph, run the buddi fit
        Failed_SS=Pymorph_list.FLAG_FAILED_S
        Failed_SE=Pymorph_list.FLAG_FAILED_SE
        SPA=pymorph_list[element].spa
        
; ADDITIONAL SELECTION OF BAD FITS
;identify bad fits from the PyMorph catalog (i.e. with Re < 1 pixel) and 
;flag them so the fits don't run
        Re_arcsec_ss=Pymorph_list.A_hl_S
        Re_arcsec_se_bulge=Pymorph_list.A_hl_SE_BULGE
        Re_arcsec_se_disk=Pymorph_list.A_hl_SE_DISK
        if (Re_arcsec_ss[element]/0.5) lt 1 then Failed_SS[element]=1
        if (Re_arcsec_se_disk[element]/0.5) lt 1 or (Re_arcsec_se_bulge[element]/0.5) lt 1 then Failed_SE[element]=1
                
; SS success AND SE failed
        if Failed_SS[element] eq 0 and Failed_SE[element] eq 1 then cntpymorphssonly += 1
; SS failed AND SE success
        if Failed_SS[element] eq 1 and Failed_SE[element] eq 0 then cntpymorphseonly += 1
; both failed
        if Failed_SS[element] eq 1 and Failed_SE[element] eq 1 then cntpymorphnone += 1
; for completeness: both success
        if Failed_SS[element] eq 0 and Failed_SE[element] eq 0 then cntpymorphboth += 1

        working_dir=root_fits+'/'+Manga_list[loop].plateifu+'/'
; SS failed OR SE success -> at least one valid model!
        if Failed_SS[element] eq 0 or Failed_SE[element] eq 0  then begin
           cntpymorph += 1
;  * Download the log-binned cube from MaNGA
           print,'galaxy = '+Manga_list[loop].plateifu
           print, 'MangaID = '+MANGA_list[loop].mangaid
           cube='manga-'+Manga_list[loop].plateifu+'-LOGCUBE'
           IF NOT file_test(root_fits+'/'+Manga_list[loop].plateifu+'/',/DIRECTORY) then begin 
              spawn,'mkdir '+ root_fits+'/'+Manga_list[loop].plateifu+'/'
           endif
 
; first get GalaxyZoo:3D file
           gzwh = where(gz3d_list.MANGAID EQ MANGA_list[loop].mangaid, cnt_gzwh)
           IF cnt_gzwh EQ 0 THEN BEGIN 
              print, 'not in GZ3D catalogue. Manga ID: '+strtrim(MANGA_list[loop].mangaid,2)
              spawn, 'touch '+working_dir+'/gz3d_mask_does_not_exist'
              openw,1,root_fits+'galaxies_not_in_gz3d.txt',/APPEND
              printf,1,Manga_list[loop].plateifu,'    ',session
              close,1
              ;print, 'skipping object'
              cntnotingz3d += 1
              ;continue
           ENDIF
           IF cnt_gzwh GT 1 THEN BEGIN
              print, 'more than one entry in GZ3D catalogue'
              stop
           ENDIF
           IF cnt_gzwh EQ 1 THEN BEGIN
              IF NOT file_test(working_dir+'/gz3d_'+strtrim(gz3d_list[gzwh].file_name,2)) THEN BEGIN
                 spawn,'rsync -avz --progress --no-motd rsync://data.sdss.org/dr17/env/MANGA_MORPHOLOGY/galaxyzoo3d/v4_0_0/gz3d_'+strtrim(gz3d_list[gzwh].file_name,2)+'.gz '+working_dir
                 spawn,'gunzip '+working_dir+'/gz3d_'+strtrim(gz3d_list[gzwh].file_name,2)+'.gz'
              ENDIF
              cntingz3d += 1
           ENDIF

;  * only download and unzip the datacube if it's not already downloaded. Note- slow on Orion
           if file_test(working_dir+cube+'.fits*') eq 0 then $
;; DR16! 
;              spawn,'rsync -avz --progress --no-motd rsync://data.sdss.org/dr16/manga/spectro/redux/v2_4_3/'+string(Manga_list[loop].plate,format='(i4)')+'/stack/'+cube+'.fits.gz '+working_dir
; DR17
              spawn,'rsync -avz --progress --no-motd rsync://data.sdss.org/dr17/manga/spectro/redux/v3_1_1/'+strtrim(Manga_list[loop].plate,2)+'/stack/'+cube+'.fits.gz '+working_dir
           if file_test(working_dir+cube+'.fits') eq 0 then $
              spawn,'gunzip '+working_dir+cube+'.fits.gz'
           
           IF keyword_set(download_only) THEN continue
           
;  * Create a PSF datacube
           if NOT file_test(working_dir+cube+'_PSF.fits') then PSF_manga,working_dir,cube,static_files      
           
;  * run the buddi_prep code to convert the units and create the bad pixel cube
           fits_read,working_dir+cube+'.fits',wave,h_wave,extname='WAVE'
           wave_log=alog10(wave)
           delvarx, wave
           if NOT file_test(working_dir+cube+'_FLUX.fits') then BUDDI_manga_prep,working_dir,cube,/JY,/BADPIX
           
;  * Create 2 BUDDI start files, one for a single Sersic fit and one for a Sersic+exponential fit
           create_buddi_start_file,working_dir,cube,Manga_list[loop].plateifu,SPA,Manga_list,Pymorph_list,loop,element,wave_log,galfit_exe,stellar_lib,/ONECOMP
           create_buddi_start_file,working_dir,cube,Manga_list[loop].plateifu,SPA,Manga_list,Pymorph_list,loop,element,wave_log,galfit_exe,stellar_lib,/TWOCOMP
           
           
;###################################
; 1 component fit
;###################################
           readcol,root_fits+'galaxies_failed_SS.txt',format='a',failed_fit_SS,comment='#',/silent
           readcol,root_fits+'galaxies_failed_SE.txt',format='a',failed_fit_SE,comment='#',/silent
           
;only run loop if galaxy hasn already failed and doesnt have a failed flag in Pymorph catalog
           temp=where(failed_fit_SS eq Manga_list[loop].plateifu)
           temp2=where(completed_SS eq Manga_list[loop].plateifu)
           
           if Failed_SS[element]  ne 1  and temp[0] eq -1 and temp2[0] eq -1 then begin
;here, we run buddi
              cntbuddi_done_ss += 1 
              print,'***Starting SS fit, galaxy '+Manga_list[loop].plateifu
              buddi,working_dir+'BUDDI_onecomp_'+Manga_list[loop].plateifu+'.txt',/AUTO
              openw,1,root_fits+'galaxies_completed_SS.txt',/APPEND
              printf,1,Manga_list[loop].plateifu
              close,1
              
;test if the fit was successful (i.e. the fit was completed)
              success=file_search(working_dir+'/IFU_decomp_1comp/image_slices/subcomps*', count=nfiles_subcomp)
              success2=file_search(working_dir+'/IFU_decomp_1comp/image_slices/galfitm*.feedme', count=nfiles_feedme)
              
;consider the fit successful if the final number of subcomps files with within 100 of the number 
;of input/feedme files (i.e. less than 100 fits failed)
              if nfiles_subcomp ge nfiles_feedme-100 and nfiles_feedme ne 0 then begin
                 buddi_manga_final_flux,working_dir+'BUDDI_onecomp_'+Manga_list[loop].plateifu+'.txt'
                 openw,1,root_fits+'galaxies_successful_SS.txt',/APPEND
                 printf,1,Manga_list[loop].plateifu
                 close,1
                 buddi_plot_results, root_fits+'/'+Manga_list[loop].plateifu+'/',Manga_list[loop].plateifu,/onecomp
;clean up the directories to save on storage
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/image_*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/badpix/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/PSF/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/sigma/*.fits'
              endif 
              
;if the fit failed, record that in the failed list
              if NOT file_test(root_fits+Manga_list[loop].plateifu+'/IFU_decomp_1comp/decomposed_data/',/DIRECTORY) and Failed_SS[element] eq 0 then begin
                 openw,1,root_fits+'galaxies_failed_SS.txt',/APPEND
                 printf,1,Manga_list[loop].plateifu+'    failed_BUDDI_fit'
                 close,1        
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/image_*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/badpix/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/PSF/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/sigma/*.fits'
              endif
           endif else if Failed_SS[element]  eq 1 then begin
              openw,1,root_fits+'galaxies_failed_SS.txt',/APPEND
              printf,1,Manga_list[loop].plateifu+'    no_input_params'
              close,1
              Failed_SS[element]=1
              spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/image_*.fits'
              spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/badpix/*.fits'
              spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/PSF/*.fits'
              spawn,'rm '+working_dir+'/IFU_decomp_1comp/image_slices/sigma/*.fits'
           endif
           
;###################################
; 2 component fit
;###################################
           
           temp=where(failed_fit_SE eq Manga_list[loop].plateifu)
           temp2=where(completed_SE eq Manga_list[loop].plateifu)
           
           if Failed_SE[element] eq 0 and temp[0] eq -1 and temp2[0] eq -1 then begin
              cntbuddi_done_se += 1 
              print,'***Starting SE fit, galaxy '+Manga_list[loop].plateifu
              
;here, we run buddi
              buddi,working_dir+'BUDDI_twocomp_'+Manga_list[loop].plateifu+'.txt',/AUTO
              openw,1,root_fits+'galaxies_completed_SE.txt',/APPEND
              printf,1,Manga_list[loop].plateifu
              close,1
              
;test if the fit was successful
              success=file_search(working_dir+'/IFU_decomp_2comp/image_slices/subcomps*', count=nfiles_subcomp)
              success2=file_search(working_dir+'/IFU_decomp_2comp/image_slices/galfitm*.feedme', count=nfiles_feedme)
              
              if nfiles_subcomp ge nfiles_feedme-10 and nfiles_feedme ne 0 then begin
                 buddi_manga_final_flux,working_dir+'BUDDI_twocomp_'+Manga_list[loop].plateifu+'.txt'
                 openw,1,root_fits+'galaxies_successful_SE.txt',/APPEND
                 printf,1,Manga_list[loop].plateifu
                 close,1
                 buddi_plot_results, root_fits+'/'+Manga_list[loop].plateifu+'/',Manga_list[loop].plateifu,/twocomp
                 
;clean up the directories to save on memory
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/image_*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/badpix/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/PSF/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/sigma/*.fits'
                 
                 
              endif 
              
              if NOT file_test(root_fits+Manga_list[loop].plateifu+'/IFU_decomp_2comp/decomposed_data/',/DIRECTORY) and Failed_SE[element] eq 0 then begin
                 openw,1,root_fits+'galaxies_failed_SE.txt',/APPEND
                 printf,1,Manga_list[loop].plateifu+'    failed_BUDDI_fit'
                 close,1
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/image_*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/badpix/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/PSF/*.fits'
                 spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/sigma/*.fits'
              endif
              
           endif else if Failed_SE[element] eq 1 then begin
              openw,1,root_fits+'galaxies_failed_SE.txt',/APPEND
              printf,1,Manga_list[loop].plateifu+'    no_input_params'
              close,1
              Failed_SE[element]=1
              spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/image_*.fits'
              spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/badpix/*.fits'
              spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/PSF/*.fits'
              spawn,'rm '+working_dir+'/IFU_decomp_2comp/image_slices/sigma/*.fits'
           endif
           
           
           
;  * Clean up by deleting all the downloaded and intermediate files
;      spawn,'rm -r '+ root_fits+'/'+Manga_list[loop].plateifu+'/'
        endif else if Failed_SS[element] ne 0 and Failed_SE[element] ne 0  then begin
;if no input paramaters in the pyMorph catalog for either fit, make a note of this issue
           openw,1,root_fits+'galaxies_failed_SE.txt',/APPEND
           printf,1,Manga_list[loop].plateifu+'    no_input_params'
           close,1
           openw,1,root_fits+'galaxies_failed_SS.txt',/APPEND
           printf,1,Manga_list[loop].plateifu+'    no_input_params'
           close,1
           cntnotpymorph += 1
        endif
;delete the original fits.gz file to save space
; do not delete first 2, otherwise it will have to be downloaded again
;    spawn,'rm '+working_dir+'/'+cube+'.fits.gz'
;    spawn,'rm '+working_dir+'/'+cube+'_*.fits'
        spawn,'rm '+working_dir+'/IFU_decomp_2comp/*.fits'
        spawn,'rm '+working_dir+'/IFU_decomp_1comp/*.fits'
        
;  * Move onto the next galaxy
        CD, root_fits
     endelse
  endfor
  print, 'no more galaxies to do'
  CD,start_dir
  
  print, '  '
  print, '  '
  print, 'total number of galaxies in MANGA catalogue: '+strtrim(totalgals,2)
  print, 'total number of galaxies in MANGA catalogue in large IFUs: '+strtrim(largegals,2)
  print, 'unique galaxies in MANGA catalogue: '+strtrim(uniquegals,2)
  print, '  '
  print, 'not in Pymorph catalogue at all: '+strtrim(cntnotpymorph,2)  
  print, 'no Pymorph SS fit: '+strtrim(cntpymorphseonly,2)  
  print, 'no Pymorph SE fit: '+strtrim(cntpymorphssonly,2)  
  print, 'no Pymorph fit at all: '+strtrim(cntpymorphnone,2)  
  print, 'both Pymorph fits: '+strtrim(cntpymorphboth,2)
  print, '  '
  print, 'with Pymorph results (searched in GZ3D): '+strtrim(cntpymorph,2)  
  print, 'no GZ3D masks: '+strtrim(cntnotingz3d,2)  
  print, 'with GZ3D masks (downloaded): '+strtrim(cntingz3d,2)  
  print, 'BUDDI attempted SS (new, this session): '+strtrim(cntbuddi_done_ss,2)  
  print, 'BUDDI attempted SE (new, this session): '+strtrim(cntbuddi_done_se,2)  

end
