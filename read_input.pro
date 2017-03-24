PRO read_input, setup_file, setup
;read in the input file for the IFU decomposition.
;
;Note- this code is based on the read_setup code distributed as part
;of the galapagos distribution (https://github.com/MegaMorph/galapagos, 
;Boris Haeussler, BorisHaeussler.astro@gmail.com)
;
;example usage
;   > read_input, setup_file, setup
;
;example:
  if file_test(setup_file) eq 0 then begin
     print, 'input file does not exist'
     stop
  ENDIF
  
  ON_IOERROR, bad_input
  
  len_num = 4                   ;length of the numbering scheme, e.g. 4 for 'A00)'
  
  
  ; 
  ; *** create structure to contail all the input information***
  ; create_struct works in pair, where the first is the name of the 
  ; variable and is followed by its value. the structure is created 
  ; using create_struct here, and then the input file is read in, 
  ; and the values put into the correct places in the structure using 
  ; the CASE keyword.
  ; 
  setup = create_struct('root', '', $
                        'galaxy_ref', '', $
                        'file', '', $
                        'kinematics', '', $
                        'decomp', '', $
                        'median_dir', '', $
                        'binned_dir', '', $
                        'slices_dir', '', $
                        'decomp_dir', '', $
                        'psf_file', '', $
                        'stellib_dir','',$
                        'stars_file','',$
                        'sigma_cube','',$
                        'badpix_cube','',$
                        'badpix_file','',$
                        'galfitm', '', $
                        'x_centre', 0., $
                        'y_centre', 0., $
                        'cont_wavelength', 0., $
                        'cont_range', 0., $
                        'targetSN', 0., $
                        'Redshift', 0., $
                        'PA', 0., $
                        'central_wavelength', 0., $
                        'wave0', 0., $
                        'step', 0., $
                        'start_wavelength', 0., $
                        'end_wavelength', 0., $
                        'no_bins', 0., $
                        'no_slices', 0., $
                        'bin_data', '', $
                        'measure_kinematics', '', $
                        'the_bias_value_is_known', '', $
                        'plot_kinematics', '', $
                        'correct_kinematics', '', $
                        'bin_datacube', '', $
                        'decompose_median_image', '', $
                        'decompose_binned_images', '', $
                        'decompose_image_slices', '', $
                        'create_subcomps', '', $
                        'create_decomposed_cubes', '', $
                        'visualise_results', '', $
                        'n_comp', 0., $
                        'constraint', '', $
                        'magzpt', 0, $
                        'sky_input', 0., $
                        'disk_type', '', $
                        'disk_mag', 0., $
                        'disk_re', 0., $
                        'disk_n', 0., $
                        'disk_q', 0., $
                        'disk_pa', 0., $
                        'disk_re_polynomial', 0., $
                        'disk_mag_polynomial', 0., $
                        'disk_n_polynomial', 0., $
                        'bulge_type', '', $
                        'bulge_mag', 0., $
                        'bulge_re', 0., $
                        'bulge_n', 0., $
                        'bulge_q', 0., $
                        'bulge_pa', 0., $
                        'bulge_re_polynomial', 0., $
                        'bulge_mag_polynomial', 0., $
                        'bulge_n_polynomial', 0., $
                        'comp3_type', '', $
                        'comp3_mag', 0., $
                        'comp3_re', 0., $
                        'comp3_n', 0., $
                        'comp3_q', 0., $
                        'comp3_pa', 0., $
                        'comp3_re_polynomial', 0., $
                        'comp3_mag_polynomial', 0., $
                        'comp3_n_polynomial', 0., $
                        'comp4_type', '', $
                        'comp4_x', 0., $
                        'comp4_y', 0., $
                        'comp4_mag', 0., $
                        'comp4_re', 0., $
                        'comp4_n', 0., $
                        'comp4_q', 0., $
                        'comp4_pa', 0.)


  
; check format for backwards compatibility
  block_bd = 0
  line = ''
  openr, 1, setup_file
  WHILE NOT eof(1) DO BEGIN
     readf, 1, line       
;get rid of leading and trailing blanks
     line = strtrim(line, 2)
;comment or empty line encountered?
     IF strmid(line, 0, 1) EQ '#' OR strlen(line) EQ 0 THEN CONTINUE       
;comment at end of line?
     pos = strpos(line, '#')
     IF pos EQ -1 THEN pos = strlen(line)
     content = strtrim(strmid(line, len_num, pos-len_num), 2)
     IF strupcase(strmid(line, 0, len_num)) eq 'G00)' then block_bd = 1
  ENDWHILE
  close, 1
  
  line = ''
  openr, 1, setup_file
  WHILE NOT eof(1) DO BEGIN
     readf, 1, line
     
;get rid of leading and trailing blanks
     line = strtrim(line, 2)
;    print, linestop
     
     
;comment or empty line encountered?
     IF strmid(line, 0, 1) EQ '#' OR strlen(line) EQ 0 THEN CONTINUE
     
;comment at end of line?
     pos = strpos(line, '#')
     IF pos EQ -1 THEN pos = strlen(line)
     
     content = strtrim(strmid(line, len_num, pos-len_num), 2)
     
     CASE strupcase(strmid(line, 0, len_num)) OF
;CHANGE=====have to trigger bad input / proper filenames
        
        'A00)': setup.root = content
        'A02)': setup.galaxy_ref = content
        'A03)': setup.file = content
        'A04)': setup.kinematics = content
        'A05)': setup.decomp = content
        'A06)': setup.median_dir = content
        'A07)': setup.binned_dir = content
        'A08)': setup.slices_dir = content
        'A09)': setup.decomp_dir = content
        'A10)': setup.psf_file = content
        'A11)': setup.stellib_dir = content
        'A12)': setup.stars_file = content
        'A13)': setup.sigma_cube = content
        'A14)': setup.badpix_cube = content
        'A15)': setup.badpix_file = content
        
        'B00)': setup.galfitm = content
         
        'C00)': setup.x_centre = float(content)
        'C01)': setup.y_centre = float(content)
        'C02)': setup.cont_wavelength = float(content)
        'C03)': setup.cont_range = float(content)
        'C04)': setup.targetSN = float(content)
        'C05)': setup.Redshift = float(content)
        'C06)': setup.PA = float(content)
        'C07)': setup.central_wavelength = float(content)
        
        'D00)': setup.wave0 = float(content)
        'D01)': setup.step = float(content)
        'D03)': setup.start_wavelength = float(content)
        'D04)': setup.end_wavelength = float(content)
        'D05)': setup.no_bins = float(content)
        'D06)': setup.no_slices = float(content)
        
        'E00)': setup.bin_data = content
        'E01)': setup.measure_kinematics = content
        'E02)': setup.the_bias_value_is_known = content
        'E03)': setup.plot_kinematics = content
        'E04)': setup.correct_kinematics = content
        'E05)': setup.bin_datacube = content
        'E06)': setup.decompose_median_image = content
        'E07)': setup.decompose_binned_images = content
        'E08)': setup.decompose_image_slices = content
        'E09)': setup.create_subcomps = content
        'E10)': setup.create_decomposed_cubes = content
        'E11)': setup.visualise_results = content
        
;        'F00)': setup.mag1 = float(content)
;        'F01)': setup.axis_ratio1 = float(content)
;        'F02)': setup.PA1 = float(content)
        
        'F00)': setup.n_comp = float(content)
        'F01)': setup.constraint = content
        'F02)': setup.magzpt = content
        'F03)': setup.sky_input = content
        
        'F10)': setup.disk_type = content
        'F11)': setup.disk_mag = float(content)
        'F12)': setup.disk_re = float(content)
        'F13)': setup.disk_n = float(content)
        'F14)': setup.disk_q = float(content)
        'F15)': setup.disk_pa = float(content)
        'F16)': setup.disk_re_polynomial = float(content)
        'F17)': setup.disk_mag_polynomial = float(content)
        'F18)': setup.disk_n_polynomial = float(content)
        
;        if setup.n_comp gt 1 then begin
          'F20)': setup.bulge_type = content
          'F21)': setup.bulge_mag = float(content)
          'F22)': setup.bulge_re = float(content)
          'F23)': setup.bulge_n = float(content)
          'F24)': setup.bulge_q = float(content)
          'F25)': setup.bulge_pa = float(content)
          'F26)': setup.bulge_re_polynomial = float(content)
          'F27)': setup.bulge_mag_polynomial = float(content)
          'F28)': setup.bulge_n_polynomial = float(content)
        
;          if setup.n_comp gt 2 then begin
          'F30)': setup.comp3_type = content
          'F31)': setup.comp3_mag = float(content)
          'F32)': setup.comp3_re = float(content)
          'F33)': setup.comp3_n = float(content)
          'F34)': setup.comp3_q = float(content)
          'F35)': setup.comp3_pa = float(content)
          'F36)': setup.comp3_re_polynomial = float(content)
          'F37)': setup.comp3_mag_polynomial = float(content)
          'F38)': setup.comp3_n_polynomial = float(content)

            'F40)': setup.comp4_type = content
            'F41)': setup.comp4_x = float(content)
            'F42)': setup.comp4_y = float(content)
            'F43)': setup.comp4_mag = float(content)
            'F44)': setup.comp4_re = float(content)
            'F45)': setup.comp4_n = float(content)
            'F46)': setup.comp4_q = float(content)
            'F47)': setup.comp4_pa = float(content)
            
;          endif
;        endif
     ENDCASE
  ENDWHILE
  close, 1
  
  
  return
  
bad_input:
  message, 'Invalid Entry in '+setup_file
END