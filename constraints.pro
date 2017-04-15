

pro constraints,dir,binned_dir,slices_dir,median_dir,n_comp,constraint,SINGLE=single,DOUBLE=double

if keyword_set(single) then begin
  for run=1,2,1 do begin
    if run eq 1 then openw,1,dir+median_dir+'galfitm.constraints'
    if run eq 2 then openw,1,dir+binned_dir+'galfitm.constraints'


 
    printf,1,'   2           q          0.1 to 1    # Soft constraint: Constrains the '
    printf,1,'                                      # sersic index n to within values '
    printf,1,'                                      # from 0.7 to 5.'
    
    printf,1,'   2             n        0.2 to 8    # Soft constraint: Constrains the '
    printf,1,'                                       # position angle to within values '
    printf,1,'                                       # from -180 to 180.'
        
    printf,1,'   2              pa        -180 to 180    # Soft constraint: Constrains the '
    printf,1,'                                       # position angle to within values '
    printf,1,'                                       # from -180 to 180.'

    close,1
  endfor 
endif

if keyword_set(double) then begin
    openw,1,dir+median_dir+'galfitm.constraints'
    openw,2,dir+binned_dir+'galfitm.constraints'
  
  
    for n=1,2,1 do begin
      
      if constraint eq 'y' then begin
        printf,n,'# Component/    parameter   constraint Comment'
        printf,n,'# operation  (see below)   range'
        printf,n,'  2_3        x          offset # Hard constraint: Constrains the'
        printf,n,'                                 # to have RELATIVE positions defined'
        printf,n,'                                 # by the initial parameter file.'
        printf,n,'   '
        printf,n,'   '
        
        printf,n,'# Component/    parameter   constraint Comment'
        printf,n,'# operation  (see below)   range'
        printf,n,'  2_3        y          offset # Hard constraint: Constrains the'
        printf,n,'                                 # to have RELATIVE positions defined'
        printf,n,'                                 # by the initial parameter file.'
        printf,n,'   '
        printf,n,'   '
      endif
      
    printf,n,'   2              pa        -180 to 180    # Soft constraint: Constrains the '
    printf,n,'                                       # position angle to within values '
    printf,n,'                                       # from -180 to 180.'
    
    if n_comp ge 1100 then begin
      printf,n,'   3              pa        -180 to 180    # Soft constraint: Constrains the '
      printf,n,'                                       # position angle to within values '
      printf,n,'                                       # from -180 to 180.'
    endif 
    
    printf,n,'   2              n        0.2 to 8    # Soft constraint: Constrains the '
    printf,n,'                                       # position angle to within values '
    printf,n,'                                       # from -180 to 180.'

    if n_comp ge 1100 then begin
      printf,n,'   3              n        0.2 to 8    # Soft constraint: Constrains the '
      printf,n,'                                       # position angle to within values '
      printf,n,'                                       # from -180 to 180.'
    endif
;      
;    printf,n,'   2              mag        5 to 20    # Soft constraint: Constrains the '
;    printf,n,'                                       # position angle to within values '
;    printf,n,'                                       # from -180 to 180.'
;      
;    printf,n,'   3              mag        5 to 20    # Soft constraint: Constrains the '
;    printf,n,'                                       # position angle to within values '
;    printf,n,'                                       # from -180 to 180.'
      
;      printf,n,'   2           q          0.1 to 1    # Soft constraint: Constrains the '
;      printf,n,'                                      # sersic index n to within values '
;      printf,n,'                                      # from 0.7 to 5.'
;      
;      printf,n,'   3           q          0.1 to 1    # Soft constraint: Constrains the '
;      printf,n,'                                      # sersic index n to within values '
;      printf,n,'                                      # from 0.7 to 5.'
      
      if n_comp ge 1110 and constraint eq 'y' then begin
        printf,n,'# Component/    parameter   constraint Comment'
        printf,n,'# operation  (see below)   range'
        printf,n,'  2_4        x          offset # Hard constraint: Constrains the'
        printf,n,'                                 # to have RELATIVE positions defined'
        printf,n,'                                 # by the initial parameter file.'
        printf,n,' '
        printf,n,' '
        printf,n,'# Component/    parameter   constraint Comment'
        printf,n,'# operation  (see below)   range'
        printf,n,'  2_4        y          offset # Hard constraint: Constrains the'
        printf,n,'                                 # to have RELATIVE positions defined'
        printf,n,'                                 # by the initial parameter file.'
      endif        

      if n_comp eq 1100 and constraint eq 'y' then begin
        printf,n,'# Component/    parameter   constraint Comment'
        printf,n,'# operation  (see below)   range'
        printf,n,'  2_3        x          offset # Hard constraint: Constrains the'
        printf,n,'                                 # to have RELATIVE positions defined'
        printf,n,'                                 # by the initial parameter file.'
        printf,n,' '
        printf,n,' '
        printf,n,'# Component/    parameter   constraint Comment'
        printf,n,'# operation  (see below)   range'
        printf,n,'  2_3        y          offset # Hard constraint: Constrains the'
        printf,n,'                                 # to have RELATIVE positions defined'
        printf,n,'                                 # by the initial parameter file.'
      endif        
      close,n
    endfor
endif

end