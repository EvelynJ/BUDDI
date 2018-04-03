;requires pros from galapagos

PRO sex2ds9reg_buddi, sexcat, regfile, radius
;add_column: 2d string array: [[label0, value0],[label1, value1],...]
;e.g.:       add_col = ['FILE', '" "']
;num does not work together with add_column - number has preference
;add_column does not work together with tag - add_column has preference

  tab = mrdfits(sexcat,1)
  print, 'catalogue contains a total of '+strtrim(n_elements(tab.alpha_J2000),2)+' objects '
  
  openw, 1, regfile
  printf, 1, '# Region file format: DS9 version 3.0'
  printf, 1, '# Filename: '+sexcat
  
  IF n_elements(color) GT 0 THEN col = color ELSE col = 'green'
  printf, 1, 'global color='+col+' font="helvetica 10 normal" ' + $
          'select=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
  printf, 1, 'image'
  
  n = n_elements(tab)
    
  FOR i=0ul, n-1 DO $
     printf, 1, 'circle('+strtrim(tab[i].x_image, 2)+','+ $
             strtrim(tab[i].y_image, 2)+','+strtrim(radius, 2)+'p)'
  
  close, 1
END
