

pro create_BUDDI_start_file,dir,cube,plate_ifu,SPA,Manga_list,Pymorph_list,loop,element,wave_log,galfit_root,stellar_lib,ONECOMP=onecomp,TWOCOMP=twocomp

close,01
if keyword_set(TWOCOMP) then openw,01,dir+'BUDDI_twocomp_'+plate_ifu+'.txt' $
  else if keyword_set(ONECOMP) then openw,01,dir+'BUDDI_onecomp_'+plate_ifu+'.txt' $
  else openw,01,dir+'BUDDI_input_onecomp_'+plate_ifu+'.txt'
  
;openw,01,dir+'BUDDI_input_'+plate_ifu+'.txt'

printf,01,'########################################################################'
printf,01,'#     INPUT FILE FOR IFU_WRAPPER.PRO'
printf,01,'#'
printf,01,'#     comments and questions should be sent to Evelyn Johnston'
printf,01,'#     (evelynjohnston.astro@gmail.com)'
printf,01,'#     ======================================================'
printf,01,'#     Version history:'
printf,01,'#     V1.0M  -  first version of the wrapper script for the automatic '
printf,01,'#                of BUDDI-MaNGA project'
printf,01,'#                EJ, Santiago, Chile 20 Feb 2020'
printf,01,'#'
printf,01,'#'
printf,01,'########################################################################'
printf,01,' '
printf,01,' '
printf,01,' '
printf,01,' '

printf,01,'# Set up directory structure for all output files'
printf,01,'A00)   '+dir+'          # [string] main data directory (root)'
printf,01,'A02)   '+plate_ifu+'                # [string] galaxy reference (galaxy_ref)'
printf,01,'A03)   '+cube+'_FLUX        # [string] input datacube, without .fits (file)'
printf,01,'A04)   kinematics/                     # [string] kinematics directory (kinematics)'
if keyword_set(TWOCOMP) then decomp_dir='IFU_decomp_2comp/' $
else if keyword_set(ONECOMP) then decomp_dir='IFU_decomp_1comp/' $
else decomp_dir='IFU_decomp_1comp/'

printf,01,'A05)   '+decomp_dir+'                     # [string] directory for IFU decomposition results (decomp)'
printf,01,'A06)   median_image/                   # [string] directory for mean galaxy image (median_dir)'
printf,01,'A07)   binned_images/                  # [string] directory for binned images (binned_dir)'
printf,01,'A08)   image_slices/                   # [string] directory for image slices (slices_dir)'
printf,01,'A09)   decomposed_data/                # [string] directory for image slices (decomp_dir)'
printf,01,'A10)   '+cube+'_PSF.fits                        # [string] input PSF datacube (psf_file)'
printf,01,'A11)   '+stellar_lib+' # [string] directory for stellar library (stellib_dir)'
printf,01,'A12)   stars.txt               # [string] positions of stars to fit  (stars_file)'
printf,01,'A13)   '+cube+'_SIGMA  # [string] sigma datacube  (sigma_cube)'
printf,01,'A14)   '+cube+'_BADPIX               # [string] bad pixel datacube  (badpix_cube)'
printf,01,'A15)   badpix.txt                # [string] text file of bad pixels (ascii)  (badpix_file)'
printf,01,' '
printf,01,' '

printf,01,'# Define software versions'
printf,01,'B00)   '+galfit_root+'   # [string] directory and executable for GalfitM (galfitm)'
printf,01,' '
printf,01,' '

h=headfits(dir+cube+'_FLUX.fits',exten=1)
x0=sxpar(h,'NAXIS1')/2.0
y0=sxpar(h,'NAXIS1')/2.0
pix_scale=abs(sxpar(h,'CD1_1')*3600)  ;convert pixel scale from degrees to arcsec

Z=Manga_list.Z
PA=Manga_list.NSA_SERSIC_PHI
if PA[loop] lt -180 or PA[loop] gt 180 then begin
  PA_temp=Pymorph_list.PA_S
  PA[loop]=-((90-PA_temp[element])-SPA[element])
  if PA[loop] ge 0 then PA[loop]-=90 else PA[loop]+=90
endif

printf,01,'# Provide basic information for the datacube'
printf,01,'C00)   '+strtrim(round(x0),2)+'           # [integer] x position of centre of galaxy (x_centre)'
printf,01,'C01)   '+strtrim(round(y0),2)+'           # [integer] y position of centre of galaxy (y_centre)'
printf,01,'C02)   6500         # [float]   central wavelength of continuum band for measuring S/N (cont_wavelength)'
printf,01,'C03)   50           # [float]   wavelength range for measuring S/N in continuum (cont_range)'
printf,01,'C04)   50           # [float]   target S/N value for binning (targetSN)'
printf,01,'C05)   '+strtrim(Z[loop],2)+'    # [float]   Redshift from NED (Redshift)'
printf,01,'C06)   '+strtrim(PA[loop],2)+'          # [float]   kinematic Position Angle of galaxy (from NED) (PA)'
printf,01,'C07)   5500         # [float]   central wavelength of spectrum for kinematics corrections (central_wavelength)'
printf,01,'C08)   2.6          # [float]   spectral resolution in AA (FWHM)'
printf,01,' '
printf,01,' '


printf,01,'# Set wavelength range for decomposition'
printf,01,'D00)   '+strtrim(sxpar(h,'CRVAL3'),2)+'      # [float]   wavelength of first pixel- log10 units (wave0)'
printf,01,'D01)   '+strtrim(sxpar(h,'CD3_3'),2)+'  # [float]   step size (step)'
printf,01,'D03)   '+strtrim(10^wave_log[0],2)+'         # [float]   start wavelength for decomposition (start_wavelength)'
printf,01,'D04)   '+strtrim(10^wave_log[-1],2)+'         # [float]   end wavelength for decomposition (end_wavelength)'
printf,01,'D05)   12           # [float]   number of bins in wavelength direction to get GalfitM fit polynomials (no_bins)'
printf,01,'D06)   10           # [float]   number of image slices to include in each set of fits (no_slices)'
printf,01,' '
printf,01,' '


printf,01,'# Determine which parts of the code to run'
if keyword_set(TWOCOMP) and file_test(dir+'kinematics/',/DIRECTORY) then begin
  printf,01,'E00)  n            # [y/n] voronoi_bin_data'
  printf,01,'E01)  n            # [y/n] measure_kinematics'
  printf,01,'E03)  n            # [y/n] plot_kinematics'
endif else if keyword_set(TWOCOMP) and ~file_test(dir+'kinematics/',/DIRECTORY) then begin
  printf,01,'E00)  y            # [y/n] voronoi_bin_data'
  printf,01,'E01)  y            # [y/n] measure_kinematics'
  printf,01,'E03)  y            # [y/n] plot_kinematics'
endif else if keyword_set(ONECOMP) then begin
  printf,01,'E00)  y            # [y/n] voronoi_bin_data'
  printf,01,'E01)  y            # [y/n] measure_kinematics'
  printf,01,'E03)  y            # [y/n] plot_kinematics'
endif
printf,01,'E04)  y            # [y/n] correct_kinematics'
printf,01,'E05)  y            # [y/n] bin_datacube'
printf,01,'E06)  y            # [y/n] decompose_median_image'
printf,01,'E07)  y            # [y/n] decompose_binned_images'
printf,01,'E08)  n            # [y/n] add GCs'
printf,01,'E09)  y            # [y/n] decompose_image_slices'
printf,01,'E10)  y            # [y/n] create_subcomps'
printf,01,'E12)  y            # [y/n] visualise_results'
printf,01,' '
printf,01,' '

if keyword_set(TWOCOMP) then ncomp='1100' $
  else if keyword_set(ONECOMP) then ncomp='1000' $
  else ncomp='1000'
  
printf,01,'#Initial estimates for single Sersic fit'
printf,01,'F00)   '+ncomp+'      # Number of components to fit'
printf,01,'F01)   y        # Constrain components 1, 2 and 3 to have the same centre?'
printf,01,'F02)   8.9      # magnitude zeropoint'
printf,01,'F03)   1e-9     # sky background'
printf,01,'F04)   n        # boxy/disky bulge (valid values: b (boxy), d (disky), n (none)'
printf,01,' '
printf,01,' '

;Absmag_g=Simard_SS.ggMag
;appmag_g=Simard_SS.gg2d
;DM=Absmag_g[loop]-appmag_g[loop]  ;calculate conversion between absolute and apparent magnitudes



if ncomp eq 1000 then begin
;  mag=Simard_SS.gg2d
;  Re_kpc=Simard_SS.Rhlg
;  Scale=Simard_SS.Scale
;  n=Simard_SS.ng
;  e=Simard_SS.e
;  PA=Simard_SS.phi
  mag=Pymorph_list.M_S
  Re_arcsec=Pymorph_list.A_hl_S
  n=Pymorph_list.N_S
  BA=Pymorph_list.BA_S
  PA_temp=Pymorph_list.PA_S
  PA=-((90-PA_temp[element])-SPA[element])
  if PA ge 0 then PA-=90 else PA+=90
  
  printf,01,'F10)   sersic    # Type of profile for disk (ALWAYS sersic)'
  printf,01,'F11)    '+strtrim(mag[element],2)+'      # magnitude estimate for disk'
;  printf,01,'F12)    '+strtrim(((Re_kpc[element]/Scale[loop])/pix_scale)*1.678,2)+'      # Re for disk'
  printf,01,'F12)    '+strtrim(Re_arcsec[element]/pix_scale,2)+'      # Re for disk'
  printf,01,'F13)    '+strtrim(n[element],2)+'       # n for disk'
  printf,01,'F14)    '+strtrim(BA[element],2)+'       # axis ratio for disk'
  printf,01,'F15)    '+strtrim(PA,2)+'       # position angle for disk'
  printf,01,'F16)    1        # polynomial order for variation in Re (negative for free)'
  printf,01,'F17)    1       # polynomial order for variation in n (negative for free)'
  printf,01,'F18)    1        # polynomial order for variation in q (negative for free)'
  printf,01,'F19)    1        # polynomial order for variation in pa (negative for free)'
  printf,01,' '
  printf,01,' '
  
  
  printf,01,'F20)   sersic    # Type of profile for bulge  (ALWAYS sersic)'
  printf,01,'F21)   15.0      # magnitude estimate for bulge'
  printf,01,'F22)   18.0       # Re for bulge'
  printf,01,'F23)   1       # n for bulge'
  printf,01,'F24)   0.8       # axis ratio for bulge'
  printf,01,'F25)   -10       # position angle for bulge'
  printf,01,'F26)    1        # polynomial order for variation in Re (negative for free)'
  printf,01,'F27)    1       # polynomial order for variation in n (negative for free)'
  printf,01,'F28)    1        # polynomial order for variation in q (negative for free)'
  printf,01,'F29)    1        # polynomial order for variation in pa (negative for free)'
  printf,01,' '
  printf,01,' '
endif

if ncomp eq 1100 then begin
;  mag=Simard_BD.gbMag
;  Re_kpc=Simard_BD.Rd
;  PA=Simard_BD.phid
;  e=Simard_SS.e
;  Scale=Simard_BD.Scale
  mag=Pymorph_list.M_SE_DISK
  Re_arcsec=Pymorph_list.A_hl_SE_DISK
  n=Pymorph_list.N_SE_DISK
  BA=Pymorph_list.BA_SE_DISK
  PA_temp=Pymorph_list.PA_SE_DISK
  PA=(90-PA_temp[element])-SPA[element]
  
  ;account for 'flipped' galaxies in the Pymorph catalog
  if n[element] ne 1.0 then begin
    mag=Pymorph_list.M_SE_BULGE
    Re_arcsec=Pymorph_list.A_hl_SE_BULGE
    n=Pymorph_list.N_SE_BULGE
    BA=Pymorph_list.BA_SE_BULGE
    PA_temp=Pymorph_list.PA_SE_BULGE
    PA=(90-PA_temp[element])-SPA[element]
  endif
  printf,01,'F10)   sersic    # Type of profile for disc  (ALWAYS sersic)'
;  printf,01,'F11)    '+strtrim(mag[element]-DM,2)+'      # magnitude estimate for disc'
  printf,01,'F11)    '+strtrim(mag[element],2)+'      # magnitude estimate for disc'
;  printf,01,'F12)    '+strtrim(((Re_kpc[element]/Scale[loop])/pix_scale)*1.678,2)+'       # Re for disc'  ;convert scale length to effective radius for exponential disc case
  printf,01,'F12)    '+strtrim(Re_arcsec[element]/pix_scale,2)+'       # Re for disc'  ;convert scale length to effective radius for exponential disc case
  printf,01,'F13)   1       # n for disc'
  printf,01,'F14)    '+strtrim(BA[element],2)+'        # axis ratio for disc'
  printf,01,'F15)    '+strtrim(PA,2)+'       # position angle for disc'
  printf,01,'F16)    1        # polynomial order for variation in Re (negative for free)'
  printf,01,'F17)    0       # polynomial order for variation in n (negative for free)'
  printf,01,'F18)    1        # polynomial order for variation in q (negative for free)'
  printf,01,'F19)    1        # polynomial order for variation in pa (negative for free)'
  printf,01,' '
  printf,01,' '


;  mag=Simard_BD.gbMag
;  Re_kpc=Simard_BD.Re
;  n=Simard_BD.nb
;  e=Simard_BD.e
;  PA=Simard_BD.phib
  mag=Pymorph_list.M_SE_BULGE
  Re_arcsec=Pymorph_list.A_hl_SE_BULGE
  n=Pymorph_list.N_SE_BULGE
  BA=Pymorph_list.BA_SE_BULGE
  PA_temp=Pymorph_list.PA_SE_BULGE
  PA=(90-PA_temp[element])-SPA[element]

  ;account for 'flipped' galaxies in the Pymorph catalog
  if n[element] eq 1.0 then begin
    mag=Pymorph_list.M_SE_DISK
    Re_arcsec=Pymorph_list.A_hl_SE_DISK
    n=Pymorph_list.N_SE_DISK
    BA=Pymorph_list.BA_SE_DISK
    PA_temp=Pymorph_list.PA_SE_DISK
    PA=(90-PA_temp[element])-SPA[element]
  endif
  
  printf,01,'F20)   sersic    # Type of profile for bulge (ALWAYS sersic)'
;  printf,01,'F21)    '+strtrim(mag[element]-DM,2)+'      # magnitude estimate for bulge'
  printf,01,'F21)    '+strtrim(mag[element],2)+'      # magnitude estimate for bulge'
  printf,01,'F22)    '+strtrim(Re_arcsec[element]/pix_scale,2)+'      # Re for bulge'
  printf,01,'F23)    '+strtrim(n[element],2)+'       # n for bulge'
  printf,01,'F24)    '+strtrim(BA[element],2)+'       # axis ratio for bulge'
  printf,01,'F25)    '+strtrim(PA,2)+'       # position angle for bulge'
  printf,01,'F26)    1        # polynomial order for variation in Re (negative for free)'
  printf,01,'F27)    1       # polynomial order for variation in n (negative for free)'
  printf,01,'F28)    1        # polynomial order for variation in q (negative for free)'
  printf,01,'F29)    1        # polynomial order for variation in pa (negative for free)'
  printf,01,' '
  printf,01,' '
  
  
endif


printf,01,' '
printf,01,'###Ignore the following lines. They are necessary for BUDDI but wont affect the fits'
printf,01,'F30)   sersic       # Type of profile for 3rd component (sersic or psf)'
printf,01,'F31)   18.0      # magnitude estimate for 3rd component'
printf,01,'F32)   13.0       # Re for 3rd component'
printf,01,'F33)   1.0       # n for 3rd component'
printf,01,'F34)   0.4       # axis ratio for 3rd component'
printf,01,'F35)   -10       # position angle for 3rd component'
printf,01,'F36)    1        # polynomial order for variation in Re (negative for free)'
printf,01,'F37)    1       # polynomial order for variation in n (negative for free)'
printf,01,'F38)    1        # polynomial order for variation in q (negative for free)'
printf,01,'F39)    1        # polynomial order for variation in pa (negative for free)'
printf,01,' '
printf,01,' '

printf,01,'F40)   sersic       # Type of profile for 3rd component (sersic or psf)'
printf,01,'F41)   18.0      # magnitude estimate for 3rd component'
printf,01,'F42)   3.0       # Re for 3rd component'
printf,01,'F43)   0.5       # n for 3rd component'
printf,01,'F44)   0.9       # axis ratio for 3rd component'
printf,01,'F45)   -10       # position angle for 3rd component'
printf,01,'F46)    1        # polynomial order for variation in Re (negative for free)'
printf,01,'F47)    1       # polynomial order for variation in n (negative for free)'
printf,01,'F48)    1        # polynomial order for variation in q (negative for free)'
printf,01,'F49)    1        # polynomial order for variation in pa (negative for free)'
printf,01,' '
printf,01,' '

printf,01,'#GC input information'
printf,01,'G00)   /usr/bin/sex    #sextractor_executable'
printf,01,'G01)   sex_standard_setup.sex    #sextractor_setup'
printf,01,'G02)   /usr/local/bin/ds9    #ds9_executable'
printf,01,' '
printf,01,' '




close,01
end
