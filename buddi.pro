; first attempt at wrapper script for IFU decomposition
; 
; This script should be started in the main directory for 
; the galaxy to be decomposed, and will start by reading 
; the input file. In the input file, you can determine which 
; parts of the code to run. It has been written such that 
; any part can be run independently, assumign that all 
; previous parts have completed successfully.
; 
; v0.1  May 2015: Code ready for ecomposing galaxies. Development 
;                 paused to focus on tests of the idea
; v0.2  Jan 2016: updated to give the option of measuring kinematics 
;                 with or without gas emission included. it is faster 
;                 without gas emission, but if there is significant 
;                 emission present then this must be measured and removed.
; v1.0  Nov 2016: starting to get BUDDI ready for public release
; 


;===============================================================

function S_N_calculator,x,y,cont_array,root,galaxy_ref,x_centre,y_centre,xy,limit
;*** Calculate S/N for each spaxel. Using the cont_wavelength and cont_range, 
;*** the Signal is calculated as the mean flux value within the defined wavelength 
;*** range, while the nise is the standard deviation in that same range. This 
;*** calculation assumes that the SD is dominated by nise and contains n 
;*** significant spectral features



S_N_array=fltarr(x,y)
close,01
openw,01,root+galaxy_ref+'_S_N_array.txt'
printf,01,'##########################################################################'
printf,01,'#          Xpix           Ypix               Signal              noise               S/N '
printf,01,'##########################################################################'

close,11
openw,11,root+galaxy_ref+'_pixel_array_full.txt'
printf,11,'##########################################################################'
printf,11,'#          Xpix           Ypix     '
printf,11,'##########################################################################'

close,21
openw,21,root+galaxy_ref+'_S_N_array_full.txt'
printf,21,'##########################################################################'
printf,21,'#          Xpix           Ypix               Signal              noise'
printf,21,'##########################################################################'

off=4
if xy[0] ne -100 and xy[1] ne -100 then background_source='y' else background_source='n'
;print,xy,background_source

for column=0,x-1,1 do begin
  for row=0,y-1,1 do begin
    Signal=mean(cont_array[column,row,*]);/100
    noise=stddev(cont_array[column,row,*]);/100
    ;if directory eq 'Commissioning_12701/' and row-y_centre ge 17 then Signal=0
    S_N=Signal/noise
    S_N_array[column,row]=S_N
    if column ge xy[0]-off-1 and column le xy[0]+off-1 and $
      row ge xy[1]-off-1 and row le xy[1]+off-1 and $
      background_source eq 'y' then begin
          Signal=0.0001
          noise=1
          ;print,column-x_centre,row-y_centre,Signal,noise
    endif
;    if Signal eq 0 then begin
;      Signal=0.0001
;      noise=1
;    endif    
    S_N=Signal/noise
    S_N_array[column,row]=S_N
    
    
    if S_N gt limit then printf,01,column-x_centre,row-y_centre,Signal,noise,Signal/noise,format='(2f16.5,3f20.4)'
    if S_N gt 0 then printf,11,column-x_centre,row-y_centre,format='(2f16.5)'
    printf,21,column-x_centre,row-y_centre,Signal,noise,format='(2f16.5,2f20.4)'
    
  endfor
endfor



close,01
close,11
close,21
return,S_N_array
end


;===============================================================
;===============================================================
pro bin2d_display_pixels, x, y, counts, pixelSize
COMPILE_OPT IDL2, HIDDEN

; Plots colored pixels with the same pixels size

PLOT, [MIN(x)-pixelSize,MAX(x)+pixelSize], [MIN(y)-pixelSize,MAX(y)+pixelSize], $
    /NODATA, /XSTYLE, /YSTYLE, XTITLE='arcsec', YTITLE='arcsec', /ISO, color=cgcolor('black')
x1 = [-0.5, -0.5, +0.5, +0.5, -0.5] * pixelSize
y1 = [+0.5, -0.5, -0.5, +0.5, +0.5] * pixelSize
color = bytscl(counts)
FOR j=0, N_ELEMENTS(counts)-1 DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]

END


;##############################################################
;##############################################################
;##############################################################
pro BUDDI,input_file

;***read in necessary information
;input_file='IFU_wrapper_input.txt'
read_input, input_file, setup



;*** Set up directory structure for all output files
root=setup.root
galaxy_ref=setup.galaxy_ref
file=setup.file
kinematics=setup.kinematics
decomp=setup.decomp
median_dir=setup.median_dir
binned_dir=setup.binned_dir
slices_dir=setup.slices_dir
psf_file=setup.psf_file
stellib_dir=setup.stellib_dir


;*** Define software versions
galfitm=setup.galfitm

;*** Provide basic information for the datacube
x_centre=fix(setup.x_centre-1)             ;x position of centre of galaxy, -1 to convert to position in array
y_centre=fix(setup.y_centre-1)             ;y position of centre of galaxy, -1 to convert to position in array
cont_wavelength=setup.cont_wavelength        ;central wavelngth of continuum band for measuring S/N
cont_range=setup.cont_range             ; wavelength range for measuring S/N in continuum
targetSN=setup.targetSN               ;target S/N value for binning
Redshift=setup.Redshift               ;value from NED after doing an position search
PA=setup.PA                     ;kinematic Position Angle of galaxy (from NED)
central_wavelength=setup.central_wavelength     ;central wavelength of spectrum for kinematics corrections
badpix_file=setup.badpix_file

;*** Set wavelength range for decomposition
start_wavelength=setup.start_wavelength       ;start wavelength for decomposition
end_wavelength=setup.end_wavelength         ;end wavelength for decomposition
no_bins=setup.no_bins
no_slices=setup.no_slices

n_comp=setup.n_comp
constraints=setup.constraint

if n_comp lt 1000 then message,'You must have at least one component (F1*)'

;for consistency, change bulge and disk types to all small letters
if setup.disk_type eq 'Sersic' then setup.disk_type='sersic' 
if setup.bulge_type eq 'Sersic' then setup.bulge_type='sersic' 
if setup.comp3_type eq 'Sersic' then setup.comp3_type='sersic' 
if setup.comp4_type eq 'Sersic' then setup.comp4_type='sersic' 
if setup.disk_type eq 'PSF' then setup.disk_type='psf' 
if setup.bulge_type eq 'PSF' then setup.bulge_type='psf' 
if setup.comp3_type eq 'PSF' then setup.comp3_type='psf' 
if setup.comp4_type eq 'PSF' then setup.comp4_type='psf' 

;;*** Initial estimates for Galfit single Sersic fit
if setup.disk_type eq 'sersic' or setup.disk_type eq 'Sersic' then disk_type=0 else disk_type=1
estimates_disk=[disk_type,setup.disk_mag,setup.disk_re,setup.disk_n,setup.disk_q,setup.disk_pa]
disk_re_polynomial=setup.disk_re_polynomial
disk_mag_polynomial=setup.disk_mag_polynomial
disk_n_polynomial=setup.disk_n_polynomial
magzpt=setup.magzpt

if n_comp ge 1100 or n_comp eq 1010 or n_comp eq 1011 then begin
  if setup.bulge_type eq 'sersic' or setup.bulge_type eq 'Sersic' then bulge_type=0 else bulge_type=1
  estimates_bulge=[bulge_type,setup.bulge_mag,setup.bulge_re,setup.bulge_n,setup.bulge_q,setup.bulge_pa] 
  bulge_re_polynomial=setup.bulge_re_polynomial
  bulge_mag_polynomial=setup.bulge_mag_polynomial
  bulge_n_polynomial=setup.bulge_n_polynomial
endif else begin 
  estimates_bulge=0
  bulge_re_polynomial=99
  bulge_mag_polynomial=99
  bulge_n_polynomial=99
endelse

if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
  if setup.comp3_type eq 'sersic' or setup.comp3_type eq 'Sersic' then comp3_type=0 else comp3_type=1
  estimates_comp3=[comp3_type,setup.comp3_mag,setup.comp3_re,setup.comp3_n,setup.comp3_q,setup.comp3_pa] 
  comp3_re_polynomial=setup.comp3_re_polynomial
  comp3_mag_polynomial=setup.comp3_mag_polynomial
  comp3_n_polynomial=setup.comp3_n_polynomial
endif else begin 
  estimates_comp3=0
  comp3_re_polynomial=99
  comp3_mag_polynomial=99
  comp3_n_polynomial=99
endelse



if setup.comp4_type eq 'sersic' or setup.comp4_type eq 'Sersic' then comp4_type=0 else  comp4_type=1
if n_comp eq 1111 or n_comp eq 1101 or n_comp eq 1001 or n_comp eq 1011 then estimates_comp4=[comp4_type,setup.comp4_x,setup.comp4_y,setup.comp4_mag,setup.comp4_re,setup.comp4_n,setup.comp4_q,setup.comp4_pa] $
else estimates_comp4=0 

bias=0.0
c = 299792.458d ; speed of light in km/s

;calculate wavelength solution. Datacube should already be log10 binned, 
;with the first wavelength value and the step size listed int he input file
;

fits_read,root+file+'.fits',input_IFU,header_IFU

wavelength=fltarr(sxpar(header_IFU,'NAXIS3'))
wavelength0=setup.wave0
step=setup.step
for m=0,sxpar(header_IFU,'NAXIS3')-1,1 do wavelength[m]=wavelength0+(m*step)



if setup.bin_data eq 'y' then begin
  ;========================================================
  ; *** Read in basic information about IFU data cube
  fits_read,root+file+'.fits',input_IFU,header_IFU
  x=sxpar(header_IFU,'NAXIS1')
  y=sxpar(header_IFU,'NAXIS2')
  z=sxpar(header_IFU,'NAXIS3')
  wavelength0=setup.wave0;sxpar(header_IFU,'CRVAL3')
  step=setup.step;sxpar(header_IFU,'CD3_3')
  
  
  ;========================================================
  ; *** identify the elements in the array corresponding to the region over 
  ; *** which to measure the S/N of each spaxel
  
  cont_wavelength_log=alog10(cont_wavelength)
  cont_low_wavelength_log=alog10(cont_wavelength-cont_range/2)
  cont_high_wavelength_log=alog10(cont_wavelength+cont_range/2)
   
  sample=where(wavelength ge cont_low_wavelength_log and wavelength le cont_high_wavelength_log)
  
  
  




  
  
  result = FILE_TEST(root+kinematics,/DIRECTORY) 
  if result eq 0 then spawn,'mkdir '+root+kinematics
  
  datacube=input_IFU
  
  ;========================================================
  ; *** Calculate S/N for each element in the array
  cont_array=datacube[*,*,sample]
  limit=5           ;only include pixels with S/N above this value in the binning
  
  ;mask regions covered by backgorund/foreground object
  if setup.n_comp eq 1001 or setup.n_comp eq 1101 or setup.n_comp eq 1111 or setup.n_comp eq 1011 then xy=[setup.comp4_x,setup.comp4_y] else xy=[-100,-100]
  
  S_N_array=S_N_calculator(x,y,cont_array,root,galaxy_ref,x_centre,y_centre,xy,limit)
  
  
  ;========================================================
  ; *** Apply voronoi binning to the image
  readcol,root+galaxy_ref+'_S_N_array.txt',format='F,F,F,F',Xpix,Ypix,signal,noise,comment='#',/SILENT
  
  ; Load a colortable and open a graphic window
;  set_plot,'x'
;  Device, Decomposed=0
;  loadct, 5
;  ;r = GET_SCREEN_SIZE()      ;ubuntu laptop
;  device,get_screen_size=r    ;macbook
;  window, xsize=r[0]*0.4, ysize=r[1]*0.8
  
  ; Perform the actual computation. The vectors
  ; (binNum, xnde, ynde, xBar, yBar, sn, nPixels, scale)
  ; are all generated in *output*
  ;signal/10 for victors datacubes
  voronoi_2d_binning, xpix, ypix, signal, noise, targetSN, $
      binNum, xnde, ynde, xBar, yBar, sn, nPixels, scale, root,galaxy_ref, /QUIET;,/PLOT
  
  ; Save to a text file the initial coordinates of each pixel together
  ; with the corresponding bin number computed by this procedure.
  ; binNum uniquely specifies the bins and for this reason it is the only
  ; number required for any subsequent calculation on the bins.
  astrolib
  forprint, xpix, ypix, binNum, TEXTOUT=root+galaxy_ref+'_voronoi_2d_binning_output.txt', $
      COMMENT='          X"              Y"           BIN_NUM'
  
  ; Print out fits file with binned spectra
  n_bins=max(binNum)+1
  binned_spec=fltarr(z,n_bins)   ;2D array for binned spectra
  
   
  
  close,03
  openw,03,root+galaxy_ref+'_bin_centers.txt'
  printf,03,'#################################################################'
  printf,03,'#       Bin       Xcentre        Ycentre '
  printf,03,'#################################################################'
  
  for bin=0,n_bins-1, 1 do begin
  ;  print,bin
  ;  print,'xBar',xBar
  ;  print,'yBar',yBar
    printf,03,bin,xBar[bin],yBar[bin],format='(i9,2f15.3)'
    pixels=where(binNum eq bin)       ;identify elements of array to be used for each bin
    spec=fltarr(z)
    for n=0,n_elements(pixels)-1,1 do begin
      m=pixels[n]
      x_new=xpix[m]+x_centre
      y_new=ypix[m]+y_centre
      spec[*]+=datacube[x_new,y_new,*]    ;coadd relevant spectra
    endfor
  
    binned_spec[*,bin]=spec[*]
    
  endfor
  close,03
  
  
  
  sxaddpar,header_IFU,'CRVAL1',sxpar(header_IFU,'CRVAL3')
  sxaddpar,header_IFU,'CD1_1',sxpar(header_IFU,'CD3_3')
  sxaddpar,header_IFU,'CDELT1',sxpar(header_IFU,'CDELT3')
  sxaddpar,header_IFU,'CRVAL2',n_elements(binned_spec[*,0])
  sxaddpar,header_IFU,'CD2_2',1
  sxaddpar,header_IFU,'CDELT2',1
  sxdelpar,header_IFU,'CRVAL3' 
  sxdelpar,header_IFU,'CDELT3' 
  sxdelpar,header_IFU,'CD3_3'
   
  fits_write,root+kinematics+galaxy_ref+'_binned_spectra.fits',binned_spec,header_IFU
  

endif else begin
  fits_read,root+kinematics+galaxy_ref+'_binned_spectra.fits',binned_spec,header_IFU
  z=sxpar(header_IFU,'NAXIS1')
  n_bins=sxpar(header_IFU,'NAXIS2')
  
 
endelse








  
;========================================================
; *** Measure the kinematics of the binned spectra using PPXF

; Read in galaxy spectrum.
; The spectrum should already be log rebinned.
;Start by measuring the kinematics of the central bin

if setup.measure_kinematics eq 'y' then begin

  a=(Redshift+1)*(Redshift+1)
  Velocity=c*((a-1)/(a+1))
  
  output=root+kinematics+galaxy_ref
  result = FILE_TEST(output+'_kinematics.txt') 
  if result eq 1 then begin
    spawn,'mv '+output+'_kinematics.txt '+output+'_old_kinematics.txt'
    spawn,'mv '+output+'_kinematics_gas.txt '+output+'_old_kinematics_gas.txt'
    spawn,'mv '+output+'_emission_info.txt '+output+'_old_emission_info.txt'
  endif

  
  ;read in first spectrum to get velocisty scale
  flux=binned_spec[*,0]

  
; Read the list of filenames from the Single Stellar Population library
;
  stellar_lib = file_search(stellib_dir+'*.fits',COUNT=nfiles_stellib)
  FWHM_tem = 2.51 ; MIUSCAT spectra have a resolution FWHM of 2.51A.
     

; Only use the wavelength range in common between galaxy and stellar library.
; Check user defined wavelength range, and if larger than that of stellar 
; templates, resample for kinematics measurements
  if nfiles_stellib ne 0 then begin      
    nfiles=nfiles_stellib
    stellib_header=headfits(stellar_lib[0])
    stellib_wave1=sxpar(stellib_header,'CRVAL1')
    stellib_step=sxpar(stellib_header,'CDELT1')
    stellib_length=sxpar(stellib_header,'NAXIS1')
    stellib_wave2=stellib_wave1+stellib_step*(stellib_length-1)
  endif else begin
    stellar_lib = file_search(stellib_dir+'*V',COUNT=nfiles_stellib_txt)
    nfiles=nfiles_stellib_txt
    readcol,stellar_lib[0],stellar_wave,stellar_flux,format='F,F'
    stellib_wave1=stellar_wave[0]
    stellib_step=stellar_wave[1]-stellar_wave[0]
    stellib_length=n_elements(stellar_wave)
    stellib_wave2=stellar_wave[stellib_length-1]      
  endelse
  
  
  if end_wavelength-start_wavelength le stellib_wave2-stellib_wave1 then $
      sample = where(wavelength gt alog10(start_wavelength) and wavelength lt alog10(end_wavelength)) $
      else sample = where(wavelength gt alog10(stellib_wave1+300) and wavelength lt alog10(stellib_wave2-300)) 

  ;sample = where(wavelength gt alog10(start_wavelength) and wavelength lt alog10(end_wavelength))
  galaxy = flux[sample]/median(flux[sample])  ; normalize spectrum to avoid numerical issues
  wave = 10^wavelength[sample]
  noise = galaxy*0 + 1           ; Assume constant noise per pixel here


; Convert velocity step below to km/s
;
    
  velScale = alog(wave[1]/wave[0])*c
  FWHM_gal = 2.52 ; MaNGA has an instrumental resolution FWHM of 2.52A.
;  FWHM_gal = FWHM_gal/(1+Redshift)   ; Adjust resolution in Angstrom


; Extract the wavelength range and logarithmically rebin one spectrum
; to the same velocity scale of the SAURON galaxy spectrum, to determine
; the size needed for the array which will contain the template spectra.
;

    if nfiles_stellib ne 0 then begin
      fits_read, stellar_lib[0], ssp, h2
      lamRange2 = sxpar(h2,'CRVAL1') + [0d,sxpar(h2,'CDELT1')*(sxpar(h2,'NAXIS1')-1d)]
    endif else begin
      readcol,stellar_lib[0],stellar_wave,ssp,format='F,F'
      lamRange2 = [stellar_wave[0],stellar_wave[stellib_length-1]]
    endelse
    log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
    stars_templates = dblarr(n_elements(sspNew),nfiles)
    

; Convolve the whole stellar library of spectral templates 
; with the quadratic difference between the galaxy and the 
; stellar instrumental resolution. Logarithmically rebin 
; and store each template as a column in the array TEMPLATES.
;
; Quadratic sigma difference in pixels stellar_lib --> SDSS
; The formula below is rigorously valid if the shapes of the 
; instrumental spectral profiles are well approximated by Gaussians. 
;

    if abs(FWHM_tem) gt abs(FWHM_gal) then stop
    FWHM_dif = SQRT(FWHM_gal^2 - FWHM_tem^2)
    sigma = FWHM_dif/2.355/stellib_step;sxpar(h2,'CDELT1') ; Sigma difference in pixels 
  

; IMPORTANT: To avoid spurious velocity offsets of the templates, the
; NPIXEL keyword in PSF_GAUSSIAN must be an odd integer as done below.
;
    lsf = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)
    for j=0,nfiles-1 do begin
        if nfiles_stellib ne 0 then fits_read, stellar_lib[j], ssp $
          else readcol,stellar_lib[j],ssp,format='X,F'
          

        ;fits_read, stellar_lib[j], ssp
        ssp = convol(ssp,lsf) ; Degrade template to SDSS resolution
        ; From IDL 8.1 one can use the following commented line instead of the 
        ; above one, and the line with PSF_GAUSSIAN is nt needed any more. 
        ; ssp = gauss_smooth(ssp,sigma)
        log_rebin, lamRange2, ssp, sspNew, logLam2, VELSCALE=velScale
        stars_templates[*,j] = sspNew/median(sspNew) ; nrmalizes templates 

    endfor




  for run=0,n_bins-1,1 do begin
    print,'*** now measuring the kinematics of binned spectrum number '+string(run,format='(i4.4)')+' out of '+string(n_bins-1,format='(i4.4)')
    flux=binned_spec[*,run]
  
    ;sample = where(wavelength gt alog10(start_wavelength) and wavelength lt alog10(end_wavelength))
    galaxy = flux[sample]/median(flux[sample])  ; normalize spectrum to avoid numerical issues
    wave = 10^wavelength[sample]
    noise = galaxy*0 + 1           ; Assume constant noise per pixel here
    mult_factor=median(flux[sample])

; Convert velocity step below to km/s
;
    
    velScale = alog(wave[1]/wave[0])*c
    FWHM_gal = 2.52 ; MaNGA has an instrumental resolution FWHM of 2.52A.

 
; If the galaxy is at a significant redshift (z > 0.03), one would need to apply 
; a large velocity shift in PPXF to match the template to the galaxy spectrum.
; This would require a large initial value for the velocity (V > 1e4 km/s) 
; in the input parameter START = [V,sig]. This can cause PPXF to stop! 
; The solution consists of bringing the galaxy spectrum roughly to the 
; rest-frame wavelength, before calling PPXF. In practice there is no 
; need to modify the spectrum in any way, given that a red shift 
; corresponds to a linear shift of the log-rebinned spectrum. 
; One just needs to compute the wavelength range in the rest-frame
; and adjust the instrumental resolution of the galaxy observations.
; This is done with the following four lines:
;
;  if Redshift ge 0.03 then begin
;    vc=Velocity/c
;    z = sqrt((1+vc)/(1-vc))-1 ; Initial estimate of the galaxy redshift
;    ;z = vc ; Initial estimate of the galaxy redshift
;    wave = wave/(1+z) ; Compute approximate restframe wavelength
;    FWHM_gal = FWHM_gal/(1+z)   ; Adjust resolution in Angstrom
;  endif



    
; The galaxy and the template spectra do nt have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV. This assume the redshift is negligible.
; In the case of a high-redshift galaxy one should de-redshift its 
; wavelength to the rest frame before using the line below (see above).
;


    dv = (logLam2[0] - alog(wave[0]))*c ; km/s
    goodPixels = ppxf_determine_goodPixels(alog(wave),lamRange2,Velocity)

; Here the actual fit starts. The best fit is plotted on the screen.
; Gas emission lines are excluded from the pPXF fit using the GOODPIXELS keyword.
;
    start = [Velocity, 2*velScale] ; (km/s), starting guess for [V,sigma]
;    start = [central_vel, central_sigma] ; (km/s), starting guess for [V,sigma]

    output=root+kinematics+galaxy_ref


;;;;Original verison of code, no gas measurements  

    ppxf_v479, stars_templates, galaxy, noise, velScale, start, sol,$
        GOODPIXELS=goodPixels, MOMENTS=2, DEGREE=4, $
        VSYST=dv, ERROR=error, BIAS=Bias, BESTFIT=bestfit;,MDEGREE=6, /PLOT

    close,50
    openw,50,output+'_kinematics.txt',/APPEND
      
      
    if run eq 0 then printf,50, '#', 'bin','V', 'sigma', 'h3', 'h4', 'h5', 'h6', FORMAT='(8A10)'
    printf,50, run,sol[0:5,0], FORMAT='(i5,f10.1,4f10.3,A10)'
    close,50

    set_plot,'ps'
    !p.multi=0
    device,file=output+'_kinematics_'+string(run,format='(i4.4)')+'.eps',/color,xoffset=0,yoffset=0
    mn = min(bestfit[goodPixels], MAX=mx)
    resid = mn + galaxy - bestfit
;    cgplot, galaxy, XTITLE='pixels', YTITLE='counts', /XSTYLE, /YNOZERO, $
;        YRANGE=[min(resid[goodPixels]),mx], YSTYLE=2, XRANGE=[-0.02,1.02]*s2[1]*s2[0]
; ;       YRANGE=[-2,mx], YSTYLE=2, XRANGE=[-0.02,1.02]*s2[1]*s2[0]
;    cgplot, bestfit, COLOR='red', THICK=2, /OVERPLOT
;    n = n_elements(goodPixels)
;    cgplot, goodPixels, replicate(mn,n*s2[0]), PSYM=3, COLOR='forest green', /OVERPLOT
;    cgplot, goodPixels, resid[goodPixels], PSYM=4, COLOR='forest green', SYMSIZE=0.3, /OVERPLOT
;    w = where((goodPixels[1:*] - goodPixels) gt 1, m)
;    for j=0,m-1 do begin
;        x = range(goodPixels[w[j]],goodPixels[w[j]+1])
;        cgplot, x, resid[x], COLOR='blue', /OVERPLOT
;    endfor
;    w = (m gt 0) ? [0,w,w+1,n-1] : [0,n-1]  ; Add first and last point
;    for j=0,n_elements(w)-1 do $
;        cgplot, goodPixels[w[[j,j]]], [mn,bestfit[goodPixels[w[j]]]], COLOR='forest green', /OVERPLOT


    cgplot, wave, galaxy, color='black', XTITLE='Observed Wavelength A', $
        YTITLE='Relative Flux', /XSTYLE, /YSTYLE, YRANGE=[-1, 1.2*max(galaxy)];, YRANGE=[-0.1, 2]
    cgoplot, wave, bestfit, color='orange'
    cgOPlot, wave, galaxy - bestfit, PSYM=4, COLOR='limegreen'

    device,/close




      
;;moments: 2= Vel and sigma only, 4= Vel, sigma, h3 and h4
;;degree: orer of polynomial to fit to continuum





;#########################################

;    print, 'Formal errors:    dV    dsigma       dh3       dh4'
;    print, error*sqrt(sol[6]), FORMAT='(10x,2f10.1,2f10.3)'

; If the galaxy is at significant redshift z and the wavelength has been
; de-redshifted with the three lines "z = 1.23..." near the beginning of 
; this procedure, the best-fitting redshift is nw given by the following 
; commented line (equation 2 of Cappellari et al. 2009, ApJ, 704, L34):
;    
;  print, 'Best-fitting redshift z:', (z + 1)*(1 + sol[0]/c) - 1
  vc=sol[0]/c
  print, 'Best-fitting redshift z:', sqrt((1+vc)/(1-vc))-1
  
    
  endfor
;close,50
;close,60
;close,70  
endif



;========================================================
; *** Use the kinematics measurements to plot the kinematics maps

if setup.plot_kinematics eq 'y' then begin
;read in kinematics
  readcol,root+kinematics+galaxy_ref+'_kinematics.txt',format='F,F,F,F,F,X,X',bin_n,vel_in,sigma_in,h3_in,h4_in,comment='#',/SILENT
  ;read in positions of each bin on the ccd image
  readcol,root+galaxy_ref+'_bin_centers.txt',format='F,F,F',bin_n_in,xbar,ybar,/silent
  
;read in pixel scale in x and y directions from file header
  fits_read,root+file+'.fits',input_IFU,header_IFU
  x_scale=abs(sxpar(header_IFU,'CD1_1')*3600)      ;*3600 to convert to arcsec
  y_scale=0.2;abs(sxpar(header_IFU,'CD2_2')*3600)
  

;calculate distance of each bin from the centre of the galaxy, taking 
;into account the inclination
  n_bins=max(bin_n_in)+1
  radius=fltarr(n_bins)
  if PA lt 0 then PA=180-PA
  if PA gt 0   and PA le 90  then case_x='case2'
  if PA gt 180 and PA le 270 then case_x='case2'
  if PA gt 90  and PA le 180 then case_x='case2'
  if PA gt 270 and PA le 360 then case_x='case2'

  if PA gt 180 then PA=PA-180
  ;if PA gt 270 then PA=PA-180
  
  for n=0,n_bins-1,1 do begin
    temp=sqrt((xbar[n]*x_scale)^2+(ybar[n]*y_scale)^2)
    ;distance(xbar[n],ybar[n],0,0) 
    ;if S/N of pixels aren't measured relative to centre of the galaxy, need 
    ;to change 0,0 to the x and y position of galaxy.
    ;need to determine positive or negative radii depending on PA position of each bin
    if case_x eq 'case1' then begin
      PA_rad=PA/360*(2*!pi)           ;convert PA into radians
      y_x=(xbar[n])*tan(PA_rad)
      if ybar[n] ge y_x then temp=temp
      if ybar[n] lt y_x then temp=-temp
    endif else begin
      alpha=(180-PA)/360*(2*!pi)           ;convert alpha into radians
      y_x=(-xbar[n])*tan(alpha)         
      ; negative sign added above to account for the fact that 
      ; in this orientation, when x is positive, the y value 
      ; of the minor axis is negative
      if ybar[n] ge y_x then temp=-temp
      if ybar[n] lt y_x then temp=temp
    endelse
    radius[n]=temp
  endfor
;  convert xbar and ybar back to pixels
;  xbar=xbar_temp
;  ybar=ybar_temp

;Make a new array for the kinematics, and rearrange in order of min to max velocity. 
; Then can remove outliers due to bad fits or foreground stars.
central_bin=where(radius eq 0.0)
if central_bin eq -1 then central_bin=where(abs(radius) eq min(abs(radius)))
vel_centre=vel_in[central_bin]
vel_centre= vel_centre[0]       ;defines vel_centre as a float, and not an array

kinematics_temp1=fltarr(8,n_bins)
kinematics_temp1[0,*]=bin_n
kinematics_temp1[1,*]=vel_in
kinematics_temp1[2,*]=sigma_in
kinematics_temp1[3,*]=h3_in
kinematics_temp1[4,*]=h4_in
kinematics_temp1[5,*]=radius

sort_by=1   ;use 1 to sort by line-of-sight velocity
kinematics_temp2=colsort(kinematics_temp1,sort_by)

kinematics_temp3=fltarr(8,n_bins)
new_central_bin=where(kinematics_temp2[1,*] eq vel_centre)

;identify any outlying points, and remove from fit
  ;outliers include those values where sigma=sigma_in or sigma_max (=1000)
badbins=fltarr(100)
badbins[*]=-1
m=0
l=0
for n=0,new_central_bin[0],1 do begin
  if (kinematics_temp2[1,n+10]-kinematics_temp2[1,n]) lt 150 then begin 
        kinematics_temp3[*,m]=kinematics_temp2[*,n]
        m+=1
  endif
  if (kinematics_temp2[1,n+10]-kinematics_temp2[1,n]) ge 150 then begin
        badbins[l]=kinematics_temp2[0,n]
        l+=1
  endif else if kinematics_temp2[2,n] eq 1000 AND badbins[l-1] ne kinematics_temp2[0,n] then begin
        badbins[l]=kinematics_temp2[0,n]
        l+=1
  endif else if kinematics_temp2[2,n] eq 80.0 AND badbins[l-1] ne kinematics_temp2[0,n] then begin
        badbins[l]=kinematics_temp2[0,n]
        l+=1
  endif
endfor

for n2=new_central_bin[0]+1,n_bins-1,1 do begin
  if (kinematics_temp2[1,n2]-kinematics_temp2[1,n2-10]) lt 150 then begin 
        kinematics_temp3[*,m]=kinematics_temp2[*,n2]
        m+=1
  endif 
  if (kinematics_temp2[1,n2]-kinematics_temp2[1,n2-10]) ge 150 then begin
        badbins[l]=kinematics_temp2[0,n2]
        l+=1
  endif else if kinematics_temp2[2,n2] eq 1000 AND badbins[l-1] ne kinematics_temp2[0,n2] then begin
        badbins[l]=kinematics_temp2[0,n2]
        l+=1
  endif else if kinematics_temp2[2,n2] eq 80.0 AND badbins[l-1] ne kinematics_temp2[0,n2] then begin
        badbins[l]=kinematics_temp2[0,n2]
        l+=1
  endif
    
endfor
non_zero=where(kinematics_temp3[1,*] gt 0,count)
kinematics_sorted=kinematics_temp3[*,0:count-1]

;Plot kinematics maps to a ps file
;  central_bin=where(radius eq 0.0)
;  vel_centre=vel_in[central_bin]
  max_vel=max(kinematics_sorted[1,*])
  min_vel=min(kinematics_sorted[1,*])
  max_sigma=max(kinematics_sorted[2,*])+50
  

  if max_vel-vel_centre gt vel_centre-min_vel then y_range=max_vel-vel_centre+100
  if max_vel-vel_centre le vel_centre-min_vel then y_range=vel_centre-min_vel+100

  set_plot,'ps'
  !p.multi=0
  loadct,0
  device,file=root+kinematics+galaxy_ref+'_kinematics.eps',xoffset=0,yoffset=0
  !p.multi=[0,2,2]
  !p.symsize=0.7
;  -130 for 
  plot,kinematics_sorted[5,*],kinematics_sorted[1,*],psym=sym(1),yrange=[vel_centre-y_range,vel_centre+y_range],ytitle='V!ILOS!N (km/s)',$
      xtitle='distance (arcsec)',xrange=[-25,25],xmargin=[10,3],ytickformat='(i5)',symsize=0.5,/ystyle,/xstyle
;  oplot,[-100,100],[velocity_NED,velocity_NED],linestyle=1
  xyouts,-24,(2*y_range)*0.9+vel_centre-y_range,galaxy_ref, charsize=0.7
  plot,kinematics_sorted[5,*],kinematics_sorted[2,*],psym=sym(1),yrange=[0,max_sigma],ytitle=greek('sigma')+' (km/s)',$
      xtitle='distance (arcsec)',xrange=[-25,25],xmargin=[10,3],/ystyle,/xstyle,symsize=0.5
  xyouts,-24,0.9*(max_sigma),galaxy_ref, charsize=0.7

  plot,kinematics_sorted[5,*],kinematics_sorted[3,*],psym=sym(1),xrange=[-25,25],yrange=[-0.3,0.3],ytitle='h!I3!N',$
      xtitle='distance (arcsec)',xmargin=[10,3],/ystyle,/xstyle,symsize=0.5
  xyouts,-24,0.9*0.6-0.3,galaxy_ref, charsize=0.7
  plot,kinematics_sorted[5,*],kinematics_sorted[4,*],psym=sym(1),xrange=[-25,25],yrange=[-0.3,0.3],ytitle='h!I4!N',$
      xtitle='distance (arcsec)',xmargin=[10,3],/ystyle,/xstyle,symsize=0.5
  xyouts,-24,0.9*0.6-0.3,galaxy_ref, charsize=0.7
  
  device,/close


; *** Plot kinematics maps
; need to first create new kinematics arrays for each pixel
  readcol,root+galaxy_ref+'_voronoi_2d_binning_output.txt',format='F,F,F',x_in,y_in,bin_new,comment='#',/SILENT
  vel_new=fltarr(n_elements(x_in))
  sigma_new=fltarr(n_elements(x_in))
  h3_new=fltarr(n_elements(x_in))
  h4_new=fltarr(n_elements(x_in))
  x_new=fltarr(n_elements(x_in))
  y_new=fltarr(n_elements(x_in))
  m=0


  for n=0,n_elements(x_in)-1,1 do begin 
    new_element=where(kinematics_sorted[0,*] eq bin_new[n],count)
    new_element=new_element[0]

    if new_element gt 0 then begin
      vel_new[m]=kinematics_sorted[1,new_element]-vel_centre
      sigma_new[m]=kinematics_sorted[2,new_element]
      h3_new[m]=kinematics_sorted[3,new_element]
      h4_new[m]=kinematics_sorted[4,new_element]
      x_new[m]=x_in[n]*x_scale
      y_new[m]=y_in[n]*y_scale
      ;if vel_new[m] gt 100 then print,n,m,vel_new[m]
      m+=1
    endif 
  endfor

  tempytemp=where(sigma_new ne 0,count)

  set_plot,'ps'
  cgloadct,33 
  device,file=root+kinematics+galaxy_ref+'_kinematics_maps.eps',/color,xoffset=0,yoffset=0
    !p.multi=[0,2,2]
    pixelSize=x_scale;1
    bin2d_display_pixels, x_new[0:count-1], y_new[0:count-1], vel_new[0:count-1], pixelSize
     cgCOLORBAR, NCOLORS=251,/right,/vertical,charsize=0.8,minrange=min(vel_new),maxrange=max(vel_new),title='V!ILOS!N (km/s)',position=[0.44,0.61,0.46,0.945],format='(i5)'
    bin2d_display_pixels, x_new[0:count-1], y_new[0:count-1], sigma_new[0:count-1], pixelSize
     cgCOLORBAR, NCOLORS=251,/right,/vertical,charsize=0.8,minrange=min(sigma_new[0:count-1]),maxrange=max(sigma_new[0:count-1]),title=greek('sigma')+' (km/s)',position=[0.97,0.61,0.99,0.945],format='(i5)'
    bin2d_display_pixels, x_new[0:count-1], y_new[0:count-1], h3_new[0:count-1], pixelSize
     ;cgCOLORBAR, NCOLORS=251,/right,/vertical,charsize=0.8,minrange=min(h3_new),maxrange=max(h3_new),title='h!I3!N',position=[0.37,0.11,0.39,0.445];,format='(f3.1)'
    bin2d_display_pixels, x_new[0:count-1], y_new[0:count-1], h4_new[0:count-1], pixelSize
     ;cgCOLORBAR, NCOLORS=251,/right,/vertical,charsize=0.8,minrange=min(h4_new),maxrange=max(h4_new),title='h!I4!N',position=[0.87,0.11,0.89,0.445];,format='(f3.1)'
    !p.multi=0

  device,/close

endif



;*** Apply kinematics corrections
; Need to identify correct bin for pixels with S/N below the limit 
; in order to apply correct kinematics corrections. To do this, I 
; should read in the centres of each bin, and calculate the distance 
; of each unbinned pixel from the bin centres to identify the closest 
; bin.
if setup.correct_kinematics eq 'y' then begin
  
  badbins=-1
  number=where(badbins ne -1,count)
  bad_bins=badbins[0:count-1]
  badpixels=fltarr(2,1000)
  badpixels[*,*]=-1

    
  ;read in positions of each bin on the ccd image
  readcol,root+galaxy_ref+'_bin_centers.txt',format='F,F,F',$
    bin_n_in,xbar,ybar,/silent

  ;read in list of binned pixels in full size array
  readcol,root+galaxy_ref+'_S_N_array_full.txt',format='F,F,F,F',$
    xpix_total,ypix_total,signal,noise,/silent

  ;read in bins for each pixel
  readcol,root+galaxy_ref+'_voronoi_2d_binning_output.txt',$
    format='F,F,F',xtemp,ytemp,binarray,/silent
  
  ;identify those pixels with S/N less than 5 and set to 0, and those 
  ;contained within the bad bins. Then work out 
  ;closest bin
  number=count
  for n=0,number-1,1 do begin
    temp_bin=where(binarray eq bad_bins[n],count)
    for m=0,count-1,1 do begin
      signal[temp_bin[m]]=1e-5     ;set bad bins to almost zero signal, but not exactly zero for later loops
      noise[temp_bin[m]]=1e-4
    endfor
  endfor
  
  x_final=fltarr(n_elements(xpix_total))
  y_final=fltarr(n_elements(xpix_total))
  bins_final=fltarr(n_elements(xpix_total))
  readcol,root+kinematics+galaxy_ref+'_kinematics.txt',$
      format='F,F,F,X,X,X,X',bin_kinematics,velocity_bin,sigma_bin,comment='#',/SILENT
  
  fits_read,root+file+'.fits',input_IFU,header_IFU
  x=sxpar(header_IFU,'NAXIS1')
  y=sxpar(header_IFU,'NAXIS2')
  z=sxpar(header_IFU,'NAXIS3')
;  wavelength0=sxpar(header_IFU,'CRVAL3')
;  step=sxpar(header_IFU,'CD3_3')
;  wavelength=fltarr(z)
;  for m=0,z-1,1 do wavelength[m]=wavelength0+(m*step)
  corrected_IFU=fltarr(x,y,z)
  intermediate_IFU=fltarr(x,y,z)

;  identify binned emission spectrum for each spaxel, and subtract
  
  tempy=where(xbar eq 0.0 and ybar eq 0.0)
  central_bin=bin_n_in[tempy]
  vel0=velocity_bin[central_bin]
  sigma0=0;sigma_bin[central_bin];max(sigma_bin)
  for j=0,n_elements(bin_n_in)-1,1 do begin
    radius=sqrt(xbar[j]^2+ybar[j]^2)
    if radius le 5 and sigma_bin[j] gt sigma0 then sigma0=sigma_bin[j]
  endfor
  
  
  vel_correction=fltarr(n_elements(bin_kinematics))
  sigma_correction=fltarr(n_elements(bin_kinematics))
  
  openw,22,root+kinematics+galaxy_ref+'_kinematics0.txt'
  printf,22,' Velocity= '+string(vel0)
  printf,22,' Sigma=    '+string(sigma0)
  close,22
  
  velScale = alog(10^wavelength[1]/10^wavelength[0])*3e5
  for l=0,n_elements(bin_kinematics)-1,1 do begin
    ;calculate velocity corrections
    wavelengthshift_linear=(velocity_bin[l]-vel0)/3e5*central_wavelength
    wavelengthshift_log=alog10(wavelengthshift_linear+central_wavelength)-alog10(central_wavelength)
    vel_correction[l]=wavelengthshift_log/step    ;wavelength correction necessary for each bin
    if vel_correction[l] ge 0 then vel_correction[l]+=0.5 else vel_correction[l]-=0.5
    vel_correction[l]=long(vel_correction[l])
    
    ;calculate sigma corrections
    if sigma_bin[l] le sigma0 then FWHM_dif = SQRT((2.355*sigma0/velScale)^2 - (2.355*sigma_bin[l]/velScale)^2)  $;divide sigma by velScale to get value in pixels
      else FWHM_dif=0.
    sigma_correction[l] = FWHM_dif/2.355*step ; Sigma difference in AA 
  endfor
  
  temparray=fltarr(3,n_elements(xpix_total))
  for m=0,n_elements(xpix_total)-1,1 do begin
    if signal[m]/noise[m] ge 5 then begin   ;identify pixels with known bins
      x_final[m]=xpix_total[m]
      y_final[m]=ypix_total[m]
      xxxxx=where(xtemp eq xpix_total[m])
      yyyyy=where(ytemp eq ypix_total[m])
      m_new=a_and_b(xxxxx,yyyyy)
      bins_final[m]=binarray[m_new[0]]
      ;subtract off emission spectrum if requested
      

      ;apply kinematics corrections here
      intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=shift(input_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*],-vel_correction[bins_final[m]])
      temparray[0,m]=x_final[m]+x_centre
      temparray[1,m]=y_final[m]+y_centre
      temparray[2,m]=-vel_correction[bins_final[m]]
      if sigma_correction[bins_final[m]] gt 0 then begin
        input_spec=intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]
        corrected_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=gauss_smooth(input_spec[*],sigma_correction[bins_final[m]]) 
;        corrected_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=gaussfold(wavelength,intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*],sigma_correction[bins_final[m]]) $
      endif else corrected_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]

      
      
    endif else if signal[m]/noise[m] lt 5 and signal[m]/noise[m] gt 0 then begin
      x_final[m]=xpix_total[m]
      y_final[m]=ypix_total[m]
      distance_to_bin=1000
      for n=0,n_elements(bin_n_in)-1,1 do begin
        deltax=x_final[m]-xbar[n]
        deltay=y_final[m]-ybar[n]
        distance_temp=sqrt(deltax*deltax+deltay*deltay)
        if distance_temp le distance_to_bin then begin
          distance_to_bin=distance_temp
          bins_final[m]=bin_n_in[n]
        endif
      endfor
      temparray[0,m]=x_final[m]+x_centre
      temparray[1,m]=y_final[m]+y_centre
      temparray[2,m]=-vel_correction[bins_final[m]]
      intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=shift(input_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*],-vel_correction[bins_final[m]])
      
      if sigma_correction[bins_final[m]] ne 0 then begin
        input_spec=intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]
        corrected_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=gauss_smooth(input_spec[*],sigma_correction[bins_final[m]]) 
;        corrected_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=gaussfold(wavelength,intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*],sigma_correction[bins_final[m]]) $
      endif else corrected_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]=intermediate_IFU[x_final[m]+x_centre,y_final[m]+y_centre,*]


     endif else begin        ;account for zero-value pixels
      x_final[m]=xpix_total[m]
      y_final[m]=ypix_total[m]
      bins_final[m]=-10
      temparray[0,m]=x_final[m]+x_centre
      temparray[1,m]=y_final[m]+y_centre
      temparray[2,m]=0
    endelse
  endfor
  
  
  result = FILE_TEST(root+decomp, /DIRECTORY) 
  if result eq 0 then file_mkdir,root+decomp
  fits_write, root+decomp+galaxy_ref+'_smoothed_kinematics.fits',corrected_IFU, header_IFU


  set_plot,'ps'
  cgloadct,33 
  device,file=root+kinematics+galaxy_ref+'_kinematics_corrections.eps',/color,xoffset=0,yoffset=0
    !p.multi=0;[0,2,2]
    pixelSize=1
    plot,temparray[0,*],temparray[1,*],/nodata,xtitle='x-axis',ytitle='y-axis',color=cgcolor('black'),xrange=[-2,max(temparray[0,*])+2],yrange=[-2,max(temparray[1,*])+2],/xstyle,/ystyle
    for n=0,n_elements(xpix_total)-1,1 do begin
      if signal[n]/noise[n] ge 5 then xyouts,temparray[0,n],temparray[1,n],string(temparray[2,n],format='(I2)'),charsize=0.4,color=cgcolor('red')
      if signal[n]/noise[n] lt 5 then xyouts,temparray[0,n],temparray[1,n],string(temparray[2,n],format='(I2)'),charsize=0.4,color=cgcolor('goldenrod')
      if signal[n] eq 0 then xyouts,temparray[0,n],temparray[1,n],string(temparray[2,n],format='(I2)'),charsize=0.4,color=cgcolor('black')
;      xyouts,temparray[0,*],temparray[1,*],string(temparray[2,*],format='(I2)'),charsize=0.4,color=colour
    endfor
  device,/close

endif





;*** Having now prepared the data cube, it can now be decomposed.

;1a. Slice the data cube up into individual images between the required wavelengths
;1b. Bin the datacube into higher quality images over the desired wavelength range
output=root

decompose='n'
if setup.bin_datacube            eq 'y' then decompose='y'  
if setup.decompose_median_image  eq 'y' then decompose='y' 
if setup.decompose_binned_images eq 'y' then decompose='y' 
if setup.decompose_image_slices  eq 'y' then decompose='y' 
if setup.create_subcomps         eq 'y' then decompose='y'  
if setup.create_decomposed_cubes eq 'y' then decompose='y'  
if setup.visualise_results       eq 'y' then decompose='y' 


if decompose eq 'y' then begin
        result = FILE_TEST(root+decomp, /DIRECTORY) 
        if result eq 0 then file_mkdir,root+decomp
        fits_read,root+decomp+galaxy_ref+'_smoothed_kinematics.fits',corrected_IFU, header_IFU
endif





if setup.bin_datacube eq 'y' then begin
    output=root
    bin_datacube,corrected_IFU, no_bins, root,decomp, binned_dir,slices_dir,median_dir,galaxy_ref,$
      file,start_wavelength, end_wavelength,  wavelength, binned_wavelengths,x_centre,y_centre,$
      psf_file,/galaxy,/PSF,/MANGA
endif

;2a. Download and decompose an sdss image of the galaxy
;   Need decomp parameters as a starting point
;   Take the median value for each pixel in the x and y plane 
;   within the wavelength range stipulated above

if decompose eq 'y' then readcol,root+decomp+slices_dir+'info.txt',format='X,F',info,/silent



if setup.decompose_median_image eq 'y' then begin
  ;first slice, last slice, number of bins, number of images ber bin, first wavelength, last wavelength
  x=sxpar(header_IFU,'NAXIS1')
  y=sxpar(header_IFU,'NAXIS2')
  print,'starting decomposition'
  median_image=fltarr(x,y)
  
  for j=0,x-1,1 do begin
    for k=0,y-1,1 do begin
      if corrected_IFU[j,k,100] eq 0 then mean=0 $
          else RESISTANT_Mean,corrected_IFU[j,k,info[0]+500:info[1]-500], 3, mean
      median_image[j,k]=mean
      
    endfor
  endfor
  
  result = FILE_TEST(root+decomp+median_dir, /DIRECTORY) 
  if result eq 0 then file_mkdir,root+decomp+median_dir
  ;h = headfits(root+galaxy_ref+'-LOGCUBE.fits');,exten='FLUX')

  fits_read,root+file+'.fits',temp_input,h
  
  
;  mkhdr, hdr0,median_image;,/IMAGE
;  
;  sxaddpar,hdr0,'EXPTIME',sxpar(h,'EXPTIME')/sxpar(h,'NEXP')
;  sxaddpar,hdr0,'NCOMBINE',sxpar(h,'NEXP')
  fits_write, root+decomp+median_dir+galaxy_ref+'_median_image.fits',median_image;,hdr0
  
  
  mkhdr, hdr0,median_image
  

  sxaddpar,hdr0,'SIMPLE','T',before='BITPIX'
  sxaddpar,hdr0,'EXPTIME',sxpar(h,'EXPTIME')/sxpar(h,'NEXP')
  sxaddpar,hdr0,'NCOMBINE',sxpar(h,'NEXP')

  
  modfits,root+decomp+median_dir+galaxy_ref+'_median_image.fits',0,hdr0
  

  
  
  
  ;fits_write, root+decomp+median_dir+galaxy_ref+'_median_image.fits',median_image,h
  ;modfits,root+decomp+median_dir+galaxy_ref+'_median_image.fits',median_image
  scale=[abs(sxpar(h,'CD1_1')*3600),abs(sxpar(h,'CD2_2')*3600)]
  
  
  ;identify bad pixels in the MaNGA data cube. Initally stick to 0-value pixels
  badpixelmask, root,decomp, galaxy_ref, badpix, binned_dir, median_dir, slices_dir, badpix_file
  
  ;write a constraints file, initially constraining the centres of the bulge 
  ;and disc fits to be together
  
  ;constraints,root+decomp,binned_dir,slices_dir,median_dir,/single
  if n_comp gt 1000 then constraints,root+decomp,binned_dir,slices_dir,median_dir,n_comp,constraints,/double


  
  ;run Galfitm for single band image to get best fit paramaters as 
  ;initial parameters for multiwaveband fits
  
  
  ;run initially for a single Sersic fit
  output=root+decomp
  
    crash_test1= file_search(output+median_dir+'imgblock_single.galfit.*',COUNT=nfiles1)
    if nfiles1 gt 0 then spawn,'rm '+output+median_dir+'imgblock_single.galfit.*'
    print,'*Now running Galfitm on single image with a single Sersic fit*'
    galfitm_single_band,output,median_dir,slices_dir,galaxy_ref,info,x,y,x_centre,$
                        y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,$
                        estimates_comp4,n_comp,disk_n_polynomial,bulge_n_polynomial,/median,/single
    CD,root+decomp+median_dir
    if n_comp eq 1000 or n_comp eq 1001 then spawn,galfitm+' galfitm_single.feedme' 

  
  ;rerun with a double Sersic profile
  if n_comp ge 1100 or n_comp eq 1011 or n_comp eq 1010 then begin
    print,'*Now running Galfitm on single image with a double Sersic fit*'
    crash_test1= file_search(output+median_dir+'imgblock_double.galfit.*',COUNT=nfiles1)
    if nfiles1 gt 0 then spawn,'rm '+output+median_dir+'imgblock_double.galfit.*'
    galfitm_single_band,output,median_dir,slices_dir,galaxy_ref,info,x,y,x_centre,$
                        y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,$
                        estimates_comp4,n_comp,disk_n_polynomial,bulge_n_polynomial,/median,/double
    CD,root+decomp+median_dir
    spawn,galfitm+' galfitm_double.feedme' 
  endif
  CD,root
  
endif

;3a. Read in the imgblock for the median image to get initial parameters.
;3b. Prepare galfitm readme file and run galfitm with all free parameters






if setup.decompose_binned_images eq 'y' then begin
  output=root+decomp
  
  fits_read,root+file+'.fits',temp_input,h
  scale=[abs(sxpar(h,'CD1_1')*3600),abs(sxpar(h,'CD2_2')*3600)]
  
  ;test to determine the polynomials for each parameter.
  ;run once with all parameters free, then ask the user
  ;if they are happy with those in the input file. If yes, 
  ;repeat with those values and ask again. If no, ask them 
  ;update the input file, then read in the new values
  
  for rep=1,99,1 do begin
    comp3_poly=[comp3_re_polynomial,comp3_mag_polynomial,comp3_n_polynomial]
      
    ;test if a free fit has already been carried out
    result=file_test(root+decomp+binned_dir+'imgblock_free.fits')
    if result eq 1 and rep eq 1 then rep=2
    
    if rep eq 1 then $
      galfitm_multiband,output,median_dir,binned_dir,slices_dir,galaxy_ref,info,x_centre,$
          y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,estimates_comp4,n_comp,no_slices,disk_re_polynomial, $
          disk_mag_polynomial,disk_n_polynomial,bulge_re_polynomial,bulge_mag_polynomial,bulge_n_polynomial,comp3_poly,$
          galfitm,rep,/binned,/header;file
    
    
    if rep ne 1 then begin
      ;if doing constrained fits, allow user to determine whether to use 
      ;results from median fit or free fit as starting parameters
      
      decision='n'
      print,'For the starting parameters, do you want to use: '
      print,'  a) results from the input file'
      print,'  b) results from the free binned fit'
      print,'  c) results from the median fit (not yet incorporated)'
      print,'  q) quit'
      while decision ne 'a' and decision ne 'b' and decision ne 'c' and decision ne 'q' do $
        READ, decision, PROMPT='Please select "a", "b", "c" or "q": '
      
      if decision eq 'q' then rep=999
      if decision eq 'a' then begin
        galfitm_multiband,output,median_dir,binned_dir,slices_dir,galaxy_ref,info,x_centre,$
          y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,estimates_comp4,n_comp,no_slices,disk_re_polynomial, $
          disk_mag_polynomial,disk_n_polynomial,bulge_re_polynomial,bulge_mag_polynomial,bulge_n_polynomial,comp3_poly,$
          galfitm,rep,/binned,/file
      endif else if decision eq 'b' then begin
        read_input, input_file, setup
        disk_re_polynomial=setup.disk_re_polynomial
        disk_mag_polynomial=setup.disk_mag_polynomial
        disk_n_polynomial=setup.disk_n_polynomial
        if n_comp ge 1100 or n_comp eq 1010 or n_comp eq 1011 then begin
          bulge_re_polynomial=setup.bulge_re_polynomial
          bulge_mag_polynomial=setup.bulge_mag_polynomial
          bulge_n_polynomial=setup.bulge_n_polynomial
        endif else begin 
          bulge_re_polynomial=99
          bulge_mag_polynomial=99
          bulge_n_polynomial=99
        endelse

        
        galfitm_multiband,output,binned_dir,binned_dir,slices_dir,galaxy_ref,info,x_centre,$
          y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,estimates_comp4,n_comp,no_slices,disk_re_polynomial, $
          disk_mag_polynomial,disk_n_polynomial,bulge_re_polynomial,bulge_mag_polynomial,bulge_n_polynomial,comp3_poly,$
          galfitm,rep,/binned,/header
        
      endif else if decision eq 'c' then begin
        galfitm_multiband,output,median_dir,binned_dir,slices_dir,galaxy_ref,info,x_centre,$
          y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,estimates_comp4,n_comp,no_slices,disk_re_polynomial, $
          disk_mag_polynomial,disk_n_polynomial,bulge_re_polynomial,bulge_mag_polynomial,bulge_n_polynomial,comp3_poly,$
          galfitm,1,/binned,/file
      endif

    endif
    
    if rep ne 999 then begin  
    
      CD,root+decomp+binned_dir
      spawn,'pwd'
      print,'*Now running Galfitm on multi-band images*'
      
      spawn,galfitm+' galfitm.feedme'
      
      ;save free fits
      if rep eq 1 then begin
        spawn,'mv imgblock.fits imgblock_free.fits'
        spawn,'mv imgblock.galfit.01 imgblock_free.galfit.01'
        spawn,'mv imgblock.galfit.01.band imgblock_free.galfit.01.band'
      endif
      
       
      
      
      ;plot summary of results to screen
      nband=info[2]
      if rep gt 1 then begin
        
        if setup.disk_type eq 'Sersic' then setup.disk_type='sersic' 
        if setup.bulge_type eq 'Sersic' then setup.bulge_type='sersic' 
        if setup.comp3_type eq 'Sersic' then setup.comp3_type='sersic' 
        if setup.comp4_type eq 'Sersic' then setup.comp4_type='sersic' 
        if setup.disk_type eq 'PSF' then setup.disk_type='psf' 
        if setup.bulge_type eq 'PSF' then setup.bulge_type='psf' 
        if setup.comp3_type eq 'PSF' then setup.comp3_type='psf' 
        if setup.comp4_type eq 'PSF' then setup.comp4_type='psf' 
        if n_comp eq 1000 then res=read_sersic_results_2comp(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        else if n_comp eq 1100 then res=read_sersic_results_2comp(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1101 and setup.comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1101 and setup.comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1001 and setup.comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        else if n_comp eq 1001 and setup.comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        
        else if n_comp eq 1010  and setup.comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        else if n_comp eq 1010  and setup.comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        else if n_comp eq 1110  and setup.comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1110  and setup.comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+binned_dir+'imgblock.fits', nband, bd=1) $
        else if n_comp eq 1011 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        else if n_comp eq 1011 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        else if n_comp eq 1011 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) $
        else if n_comp eq 1011 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+binned_dir+'imgblock.fits', nband, bd=0) 
      
      
        Re_disk=res.RE_GALFIT_BAND_D
        mag_disk=res.MAG_GALFIT_BAND_D
        n_disk=res.N_GALFIT_BAND_D
        pa_disk=res.PA_GALFIT_BAND_D
        q_disk=res.Q_GALFIT_BAND_D
        if n_comp ge 1100 then begin
          Re_bulge=res.RE_GALFIT_BAND_B
          mag_bulge=res.MAG_GALFIT_BAND_B
          n_bulge=res.N_GALFIT_BAND_B
          pa_bulge=res.PA_GALFIT_BAND_B
          q_bulge=res.Q_GALFIT_BAND_B
        endif
        if n_comp gt 1100 and n_comp ne 1101 then begin
          Re_comp3=res.RE_GALFIT_BAND_COMP3
          mag_comp3=res.MAG_GALFIT_BAND_COMP3
          n_comp3=res.N_GALFIT_BAND_COMP3
          pa_comp3=res.PA_GALFIT_BAND_COMP3
          q_comp3=res.Q_GALFIT_BAND_COMP3
        endif
        sky=res.sky_galfit_band
      endif
      
;      if rep eq 1 then begin
        if n_comp eq 1000 then res=read_sersic_results_2comp(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        else if n_comp eq 1100 then res=read_sersic_results_2comp(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1101 and setup.comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1101 and setup.comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1001 and setup.comp4_type eq 'psf' then res=read_sersic_results_3psf(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        else if n_comp eq 1001 and setup.comp4_type eq 'sersic' then res=read_sersic_results_3sersic(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        
        else if n_comp eq 1010  and setup.comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        else if n_comp eq 1010  and setup.comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        else if n_comp eq 1110  and setup.comp3_type eq 'psf' then res=read_sersic_results_2comp_p(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1110  and setup.comp3_type eq 'sersic' then res=read_sersic_results_2comp_s(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1111 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=1) $
        else if n_comp eq 1011 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'psf' then res=read_sersic_results_3psf_p(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        else if n_comp eq 1011 and setup.comp4_type eq 'psf' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3psf_s(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        else if n_comp eq 1011 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'psf' then res=read_sersic_results_3sersic_p(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) $
        else if n_comp eq 1011 and setup.comp4_type eq 'sersic' and setup.comp3_type eq 'sersic' then res=read_sersic_results_3sersic_s(root+decomp+binned_dir+'imgblock_free.fits', nband, bd=0) 
        
        
        Re_disk1=res.RE_GALFIT_BAND_D
        mag_disk1=res.MAG_GALFIT_BAND_D
        n_disk1=res.N_GALFIT_BAND_D
        pa_disk1=res.PA_GALFIT_BAND_D
        q_disk1=res.Q_GALFIT_BAND_D
        if n_comp ge 1100 then begin
          Re_bulge1=res.RE_GALFIT_BAND_B
          mag_bulge1=res.MAG_GALFIT_BAND_B
          n_bulge1=res.N_GALFIT_BAND_B
          pa_bulge1=res.PA_GALFIT_BAND_B
          q_bulge1=res.Q_GALFIT_BAND_B
        endif
        if n_comp gt 1100 and n_comp ne 1101 then begin
          Re_comp3_1=res.RE_GALFIT_BAND_COMP3
          mag_comp3_1=res.MAG_GALFIT_BAND_COMP3
          n_comp3_1=res.N_GALFIT_BAND_COMP3
          pa_comp3_1=res.PA_GALFIT_BAND_COMP3
          q_comp3_1=res.Q_GALFIT_BAND_COMP3
        endif
        sky1=res.sky_galfit_band

          
      set_plot,'ps'
      device,file=root+decomp+binned_dir+'summary_plots_'+string(rep,format='(I2.2)')+'.eps',/landscape;xoffset=0,yoffset=0,xsize=11,ysize=8,/inches,/color;,/landscape
      !P.thick=2
      !p.charthick=2
      !p.charsize=1
      !p.multi=0;[0,2,3] 
      
      if n_comp eq 1110 or n_comp eq 1101 then number=3  $
        else if n_comp eq 1100 then number=2 $
        else if n_comp eq 1000 then number=1
      multiplot, [number,5],  mXtitle='Wavelength ('+cgSymbol("angstrom")+')'
      symbolsize=0.5
    
      x1=3500
      x2=10500
      if n_comp ge 1100 then begin
        mag_max=max([mag_bulge1[1:-2],mag_disk1[1:-2]])+1.5
        mag_min=min([mag_bulge1[1:-2],mag_disk1[1:-2]])-1.5
        Re_max=max([Re_bulge1[1:-2],Re_disk1[1:-2]])+4
        n_max=max([n_bulge1[1:-2],n_disk1[1:-2]])+1
        PA_min=min([pa_bulge1[1:-2],pa_disk1[1:-2]])-10
        PA_max=max([pa_bulge1[1:-2],pa_disk1[1:-2]])+10
      endif else begin
        mag_max=max(mag_disk1[1:-2])+1.5
        mag_min=min(mag_disk1[1:-2])-1.5
        Re_max=max(Re_disk1[1:-2])+4
        n_max=max(n_disk1[1:-2])+1
        PA_min=min(pa_disk1[1:-2])-10
        PA_max=max(pa_disk1[1:-2])+10
      endelse
      plot,mag_disk1,yrange=[mag_max,mag_min],$
        ytitle='m!ITot!N',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,title='Disc'
      if rep gt 1 then oplot,mag_disk,linestyle=0
      multiplot
      if n_comp ge 1100 then begin
        plot,mag_bulge1,yrange=[mag_max,mag_min],$
          /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,$
          title='Bulge'
        if rep gt 1 then oplot,mag_bulge,linestyle=0
        multiplot
      endif
      if n_comp gt 1100 then begin
        plot,mag_comp3_1,yrange=[mag_max,mag_min],$
          /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,title='Component 3'
        if rep gt 1 and n_comp gt 1100 then oplot,mag_comp3,linestyle=0
        multiplot
      endif
      
      plot,Re_disk1,yrange=[0,Re_max],$
        ytitle='R!Ie!N!C(arcsec)',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,ytickinterval=4,yminor=4
      if rep gt 1 then oplot,Re_disk,linestyle=0
      multiplot
      if n_comp ge 1100 then begin
        plot,Re_bulge1,yrange=[0,Re_max],$
          /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,$
          ytickinterval=4,yminor=4
        if rep gt 1 then oplot,Re_bulge,linestyle=0
        multiplot
      endif
      if n_comp gt 1100 then begin
        plot,Re_comp3_1,yrange=[0,Re_max],$
          /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,ytickinterval=4,yminor=4
        if rep gt 1 then oplot,Re_comp3,linestyle=0
        multiplot
      endif
      
      plot,n_disk1,yrange=[0,n_max],$
        ytitle='n',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,ytickinterval=2
      if rep gt 1 then oplot,n_disk,linestyle=0
      multiplot
      if n_comp ge 1100 then begin
        plot,n_bulge1,yrange=[0,n_max],$
          /xstyle,/ystyle,psym=sym(1),$
          ytickinterval=2
         if rep gt 1 then oplot,n_bulge,linestyle=0
         multiplot
      endif
      if n_comp gt 1100 then begin
        plot,n_comp3_1,yrange=[0,n_max],$
          /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,ytickinterval=2
        if rep gt 1 then oplot,n_comp3,linestyle=0
        multiplot
      endif
      
      plot,pa_disk1,yrange=[PA_min,PA_max],$
        ytitle='PA!C(degrees)',/xstyle,/ystyle,ytickinterval=20,psym=sym(1),symsize=symbolsize
      if rep gt 1 then oplot,pa_disk,linestyle=0
      multiplot
      
      if n_comp ge 1100 then begin 
        plot,pa_bulge1,yrange=[PA_min,PA_max],$
         /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,$
         ytickinterval=20
         if rep gt 1 then oplot,pa_bulge,linestyle=0
         multiplot
       endif
       if n_comp gt 1100 then begin
         plot,pa_comp3_1,yrange=[PA_min,PA_max],$
           /xstyle,/ystyle,ytickinterval=20,psym=sym(1),symsize=symbolsize
         if rep gt 1 then oplot,pa_comp3,linestyle=0
         multiplot
       endif
       
       plot,q_disk1,yrange=[0,1],$
         ytitle='q',/xstyle,/ystyle,psym=sym(1),symsize=symbolsize,ytickinterval=0.5
       if rep gt 1 then oplot,q_disk,linestyle=0
       ;multiplot,/reset,/default
      if n_comp ge 1100 then begin
        multiplot
        plot,q_bulge1,yrange=[0,1],$
          /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,$
          ytickinterval=0.5
        if rep gt 1 then oplot,q_bulge,linestyle=0
        
      endif
      if n_comp gt 1100 then begin
        multiplot
        plot,q_comp3_1,yrange=[0,1],$
          /xstyle,/ystyle,psym=sym(1),symsize=symbolsize,ytickinterval=0.5
        if rep gt 1 then oplot,q_comp3,linestyle=0
      endif
      multiplot,/reset,/default
      device,/close      
      CD,root
      
      
      ;user inpiut to determine if polynomials are ok
      print,'******************************************************'
      print,'The latest fit has been plotted and saved as'
      print,root+decomp+binned_dir+'summary_plots_'+string(rep,format='(I2.2)')+'.eps'
      print,'******************************************************'
      
      decision='0'
      while decision ne 'y' and decision ne 'n' do $
        READ, decision, PROMPT='Are you happy with your current polynomial fits? [y/n]: '
      
      if decision eq 'y' and rep gt 1 then rep=1000 $
      else if decision eq 'n' then begin
        delay='n'
;        print,delay,'Please update the input file, and type "y" when done.  '
  
        ;;*** Initial estimates for Galfit single Sersic fit
        if setup.disk_type eq 'sersic' or setup.disk_type eq 'Sersic' then disk_type=0 else disk_type=1
        estimates_disk=[disk_type,setup.disk_mag,setup.disk_re,setup.disk_n,setup.disk_q,setup.disk_pa]
        disk_re_polynomial=setup.disk_re_polynomial
        disk_mag_polynomial=setup.disk_mag_polynomial
        disk_n_polynomial=setup.disk_n_polynomial
        if n_comp ge 1100 or n_comp eq 1010 or n_comp eq 1011 then begin
          if setup.bulge_type eq 'sersic' or setup.bulge_type eq 'Sersic' then bulge_type=0 else bulge_type=1
          estimates_bulge=[bulge_type,setup.bulge_mag,setup.bulge_re,setup.bulge_n,setup.bulge_q,setup.bulge_pa] 
          bulge_re_polynomial=setup.bulge_re_polynomial
          bulge_mag_polynomial=setup.bulge_mag_polynomial
          bulge_n_polynomial=setup.bulge_n_polynomial
        endif else begin 
          estimates_bulge=0
          bulge_re_polynomial=99
          bulge_mag_polynomial=99
          bulge_n_polynomial=99
        endelse
        
        if n_comp eq 1010 or n_comp eq 1011 or n_comp eq 1110 or n_comp eq 1111 then begin
          if setup.comp3_type eq 'sersic' or setup.comp3_type eq 'Sersic' then comp3_type=0 else comp3_type=1
          estimates_comp3=[comp3_type,setup.comp3_mag,setup.comp3_re,setup.comp3_n,setup.comp3_q,setup.comp3_pa] 
          comp3_re_polynomial=setup.comp3_re_polynomial
          comp3_mag_polynomial=setup.comp3_mag_polynomial
          comp3_n_polynomial=setup.comp3_n_polynomial
        endif else begin 
          estimates_comp3=0
          comp3_re_polynomial=99
          comp3_mag_polynomial=99
          comp3_n_polynomial=99
        endelse
  
        
      endif
    endif
  endfor
  
  
  
endif

;4a. Read in the imgblock for the binned images to get chebychev polynomials.
;4b. Prepare galfitm readme file and run galfitm with all free parameters

if setup.decompose_image_slices eq 'y' then begin
  output=root+decomp
  comp3_poly=[comp3_re_polynomial,comp3_mag_polynomial,comp3_n_polynomial]
  galfitm_multiband,output,median_dir,binned_dir,slices_dir,galaxy_ref,info,x_centre,$
    y_centre,scale,magzpt,estimates_bulge,estimates_disk,estimates_comp3,estimates_comp4,n_comp,no_slices,disk_re_polynomial, $
    disk_mag_polynomial,disk_n_polynomial,bulge_re_polynomial,bulge_mag_polynomial,bulge_n_polynomial,comp3_poly,galfitm,/slices

  result = FILE_TEST(root+decomp+slices_dir+'galfit.*') 
  if result eq 1 then spawn,'rm '+root+decomp+slices_dir+'galfit.*'
  
  
  CD,root+decomp+slices_dir;decomp+median_dir
  spawn,'pwd'
  print,'*Now running Galfitm on image_slices*'
  spawn,'rm *galfit.01*'
  spawn,'chmod a+x run_galfitm.sh'
  
  spawn,'./run_galfitm.sh'
  print,'***Note: PSF binning and implementation has not yet been coded***'
  
  CD,root
 
endif
;  print,'if you see this, then Im testing where the code has hung'


;5. Run subcomps program to obtain bulge and disc images
if setup.create_subcomps eq 'y' then begin
;  print,'if you see this, then the code has not hung'

  subcomps,root+decomp,slices_dir,decomp,info,no_slices,n_comp,setup.comp3_type,setup.comp4_type,wavelength,galfitm,/MANGA

  CD,root+decomp+slices_dir;+'models/';decomp+median_dir
  spawn,'pwd'
  print,'*Now running Galfitm -o3 to get subcomps*'
  spawn,'chmod a+x run_subcomps.sh'
  spawn,'./run_subcomps.sh'
  
  CD,root
endif

;6. create bulge and disc data cubes

;need to make sure the code can identify where no subcomps file exists, 
;say in the case that galfitm crashes for that feedme file

if setup.create_decomposed_cubes eq 'y' then datacube_creator,root,decomp,kinematics,galaxy_ref,file,slices_dir,info,n_comp,setup.comp3_type,setup.comp4_type,no_slices,wavelength,/MANGA


;7. Analysis of the datacubes.
;7a. Use datacubes to create !D bulge and disc spectra, and ps file to visualise results.
;

if setup.visualise_results eq 'y' then result_visualiser,root,decomp,galaxy_ref,slices_dir,info,x_centre,y_centre,start_wavelength,end_wavelength,wavelength,Redshift,n_comp,setup.comp3_type,setup.comp4_type,setup.comp4_x,setup.comp4_y,no_slices,/MANGA



print,'*** Code has finished running***'
end

