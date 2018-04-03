PRO add_1_component, galfit_sf, name_add, new_output_filename, new_constraints_name, x_pos, y_pos, mag, re, n, ar, PA, type
;  this code
;    adds one component (SERSIC, GAUSS or PSF) to a galfit start file.
;    use constraints for previously established parameters??
; !!   currently holds all values constant with wavelength

; currently works only with PSF, GAUSS and SERSIC profiles. More would
; need to be added. Technically, GC are KING profiles, but those are
; harder to use!

; in case of GAUSS, the re mean FWHM!!
  
; derive and change folder
  spawn, 'pwd', infolder
  IF strpos(galfit_sf,'/') NE -1 THEN BEGIN
     workfolder = strmid(galfit_sf,0,strpos(galfit_sf,'/',/REVERSE_SEARCH))
     CD, workfolder
     obj = strmid(galfit_sf,strpos(galfit_sf,'/',/REVERSE_SEARCH)+1)
  ENDIF

; go through setup file one line at a time and change line with output file format to include all information
  openr, filer, galfit_sf, /get_lun
  line = ''

  openw, filew, galfit_sf+'_'+name_add, /get_lun
  comp_count = 0
  
  WHILE ~ EOF(filer) DO BEGIN
; Read a line of text:
     readf, filer, line
; in first setup line, count number of images
     IF strpos(strtrim(line, 2), 'A) ') EQ 0 THEN BEGIN
        start = strtrim(strmid(line,0,strpos(line,')')+2),2)
        content = strtrim(strmid(line,strpos(line,')')+2, strpos(line,'#')-strpos(line,')')-2),2)
        comment = strtrim(strmid(line,strpos(line,'#')),2)
        content_elements = strsplit(content,',',/extract)
        nband = n_elements(content_elements)
        if nband eq 1 then bandstr = ' '
        if nband gt 1 then bandstr = 'band'
     ENDIF

; change output image name
     IF strpos(strtrim(line, 2), 'B) ') EQ 0 THEN BEGIN
        start = strtrim(strmid(line,0,strpos(line,')')+2),2)
        content = strtrim(strmid(line,strpos(line,')')+2, strpos(line,'#')-strpos(line,')')-2),2)
        comment = strtrim(strmid(line,strpos(line,'#')),2)
        printf, filew, start+' '+new_output_filename+'            '+comment
     ENDIF ELSE IF strpos(strtrim(line, 2), 'G) ') EQ 0 THEN BEGIN
        start = strtrim(strmid(line,0,strpos(line,')')+2),2)
        content = strtrim(strmid(line,strpos(line,')')+2, strpos(line,'#')-strpos(line,')')-2),2)
        comment = strtrim(strmid(line,strpos(line,'#')),2)
        printf, filew, start+' '+new_constraints_name+'            '+comment
     ENDIF ELSE printf, filew, line

     IF strpos(strtrim(line, 2), '0) ') EQ 0 THEN BEGIN
        start = strtrim(strmid(line,0,strpos(line,')')+2),2)
        content = strtrim(strmid(line,strpos(line,')')+2, strpos(line,'#')-strpos(line,')')-2),2)
        comment = strtrim(strmid(line,strpos(line,'#')),2)
        comp_count += 1
     ENDIF
   
  ENDWHILE
  close, filer
  FREE_LUN, filer

; at this point, the new file should be an identical copy to the old one
; now add the component

; now add the wanted target!
; X Position
  x_po = string(strarr(nband)+x_pos,format = '('+string(nband)+'(A,","))')
  x_po = strtrim(strcompress(x_po, /remove_all), 2)
  x_po = ' 1) '+strmid(x_po, 0, strlen(x_po)-1)+'  1  '+bandstr+'   # position x     [pixel]'

; Y Position
  y_po = string(strarr(nband)+y_pos,format = '('+string(nband)+'(A,","))')
  y_po = strtrim(strcompress(y_po, /remove_all), 2)
  y_po = ' 2) '+strmid(y_po, 0, strlen(y_po)-1)+'  1  '+bandstr+'   # position y     [pixel]'

; MAG 
  mag_po = string(strarr(nband)+mag,format = '('+string(nband)+'(A,","))')
  mag_po = strtrim(strcompress(mag_po, /remove_all), 2)
  mag_po = ' 3) '+strmid(mag_po, 0, strlen(mag_po)-1)+'  '+strtrim(nband-2,2)+'  '+bandstr+'   # magnitude'

; RE 
  re_po = string(strarr(nband)+re,format = '('+string(nband)+'(A,","))')
  re_po = strtrim(strcompress(re_po, /remove_all), 2)
  re_po = ' 4) '+strmid(re_po, 0, strlen(re_po)-1)+'  1  '+bandstr+'   # halflight_radius   [pixel]'

; SERSIC INDEX 
  n_po = string(strarr(nband)+n,format = '('+string(nband)+'(A,","))')
  n_po = strtrim(strcompress(n_po, /remove_all), 2)
  n_po = ' 5) '+strmid(n_po, 0, strlen(n_po)-1)+'  1  '+bandstr+'   # sersic index'

; AR 
  q_po = string(strarr(nband)+ar,format = '('+string(nband)+'(A,","))')
  q_po = strtrim(strcompress(q_po, /remove_all), 2)
  q_po = ' 9) '+strmid(q_po, 0, strlen(q_po)-1)+'  1  '+bandstr+'   # axis ratio'

; PA 
  pa_po = string(strarr(nband)+pa,format = '('+string(nband)+'(A,","))')
  pa_po = strtrim(strcompress(pa_po, /remove_all), 2)
  pa_po = '10) '+strmid(pa_po, 0, strlen(pa_po)-1)+'  1  '+bandstr+'   # position angle'
   

; create new constraints file
  openr, filerc, new_constraints_name, /get_lun
  openw, filewc, new_constraints_name+'_progress', /get_lun
; copy each line from old file
  WHILE ~ EOF(filerc) DO BEGIN
; Read a line of text:
     readf, filerc, line
     printf, filewc, line
  ENDWHILE
  close, filerc
  FREE_LUN, filerc
  
  IF strtrim(type) EQ 'sersic' THEN BEGIN
     printf, filew
     printf, filew, ' 0) sersic             # Object type --- SERSIC'
     printf, filew, x_po
     printf, filew, y_po
     printf, filew, mag_po
     printf, filew, re_po
     printf, filew, n_po
     printf, filew, q_po
     printf, filew, pa_po
     printf, filew, ' Z) 0                  # output image (see above)'
; add constraints
     printf, filewc, strtrim(comp_count,2)+'     x      -3    3'
     printf, filewc, strtrim(comp_count,2)+'     y      -3    3'
     printf, filewc, strtrim(comp_count,2)+'     n       0.2 to 8'
     printf, filewc, strtrim(comp_count,2)+'     pa     -180 to 180'
  ENDIF
  
  IF strtrim(type) EQ 'psf' THEN BEGIN
     printf, filew
     printf, filew, ' 0) psf             # Object type --- PSF'
     printf, filew, x_po
     printf, filew, y_po
     printf, filew, mag_po
     printf, filew, ' Z) 0                  # output image (see above)'
; add constraints
     printf, filewc, strtrim(comp_count,2)+'     x      -3    3'
     printf, filewc, strtrim(comp_count,2)+'     y      -3    3'
  ENDIF
  
  IF strtrim(type) EQ 'gaussian' THEN BEGIN
     printf, filew
     printf, filew, ' 0) gaussian             # Object type --- Gaussian'
     printf, filew, x_po
     printf, filew, y_po
     printf, filew, mag_po
     printf, filew, re_po
     printf, filew, q_po
     printf, filew, pa_po
     printf, filew, ' Z) 0                  # output image (see above)'
; add constraints
     printf, filewc, strtrim(comp_count,2)+'     x      -3    3'
     printf, filewc, strtrim(comp_count,2)+'     y      -3    3'
  ENDIF
  
;  IF strtrim(type) EQ 'king' THEN BEGIN
;   printf, filew
;   printf, filew, ' 0) king             # Object type --- KING'
;   printf, filew, x_po
;   printf, filew, y_po
;   printf, filew, mag_po
;
; 3) 14.9805    1       #  mu(0) 
; 4) 10.1328    1       #   Rc 
; 5) 51.0968    1       #   Rt 
; 6) 2.0485     1       # alpha 
; 9) 0.9918     1       # axis ratio (b/a)  
;10) 20.7684    1       # position angle (PA) 
;   printf, filew, ' Z) 0                  # output image (see above)'
;  ENDIF

  close, filew
  FREE_LUN, filew
  close, filewc
  FREE_LUN, filewc

; copy work constraints file to final name
  spawn, 'mv '+new_constraints_name+'_progress '+new_constraints_name

  
  CD, infolder

END
