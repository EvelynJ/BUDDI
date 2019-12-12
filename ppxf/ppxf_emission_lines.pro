;#############################################################################
;
; Generates an array of Gaussian emission lines to be used as templates in PPXF.
; Additional lines can be easily added by editing this procedure.
;
; - logLam_temp is the natural log of the wavelength of the templates in Angstrom.
; logLam_temp should be the same as that of the stellar templates.
;
; - lamRange_gal is the estimated rest-frame fitted wavelength range
; Typically lamRange_gal = [min(wave), max(wave)]/(1 + z),
; where wave is the observed wavelength of the fitted galaxy pixels
; and z is an initial very rough estimate of the galaxy redshift.
;
; - FWHM_gal is the instrumantal FWHM of the galaxy spectrum under study in
; Angstrom. Here it is assumed constant. It could be a function of wavelength.
;
;- The [OI], [OIII] and [NII] doublets are fixed at theoretical flux ratio~3.
;
; V1.0.0: Michele Cappellari, Oxford, 7 January 2014
; V1.0.1: Also return line names and wavelengths. MC, Oxford, 16 November 2015
; V1.1.0: Treat line doublets as a single template, as in Python version.
;         Only returns lines included within the estimated fitted wavelength range. 
;         MC, Oxford, 25 January 2016
; V1.2.0: Perform integration over the pixels of the Gaussian line spread function
;         using the new function emline(). Thanks to Eric Emsellem for the suggestion.
;         MC, Oxford, 10 August 2016
;       
;#############################################################################

function emline, logLam_temp, line_wave, sigma

; Instrumental Gaussian line spread function integrated within the
; pixels borders. The function is normalized in such a way that
;
;     total(integ) = 1
;
; For sigma larger than one pixels, this function quickly converges
; to the normalized Gaussian function:
;
;     gauss = dLogLam * exp(-0.5*(x/xsig)**2) / (sqrt(2*!pi)*xsig)
;
; :param logLam_temp: alog(wavelength) in Angstrom
; :param line_wave: lines wavelength in Angstrom
; :param sigma: sigma in Angstrom
; :return: LSF computed for every logLam_temp

; Compute pixels borders for Gaussian integration
n = N_ELEMENTS(logLam_temp)
logLamBorders = (logLam_temp[1:n-1] + logLam_temp[0:n-2])/2
xsig = sigma/line_wave    ; sigma in logLambda units

; Perform pixel integration using
; equation (28) of Cappellari (2017, MNRAS, 466, 798)
x = logLamBorders - alog(line_wave)
integ = 0.5*cap_diff(erf(x/(sqrt(2)*xsig)))

return, [0, integ, 0]
end
;------------------------------------------------------------------------------

function ppxf_emission_lines, logLam_temp, lamRange_gal, FWHM_gal, $
    LINE_NAMES=line_names, LINE_WAVE=line_wave

sigma = FWHM_gal/2.355  ; Here sigma is assumed constant. But it can vary with wavelength

; Second dimension is the number of gas templates (lines + doublets)
emission_lines = dblarr(n_elements(logLam_temp), 12) 
k = 0  ; Indext of the column of the templates array

; Balmer Series: Hdelta Hgamma   Hbeta     Halpha
line_wave = [4101.76, 4340.47, 4861.33, 6562.80]   ; air wavelengths
line_names = ['Hdelta', 'Hgamma', 'Hbeta', 'Halpha']
for j=0,n_elements(line_wave)-1 do $
    emission_lines[*, k++] = emline(logLam_temp, line_wave[j], sigma)

;         -----[OII]-----    -----[SII]-----
lines = [3726.03, 3728.82, 6716.47, 6730.85]    ; air wavelengths
names = ['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731']
for j=0,n_elements(lines)-1 do $
    emission_lines[*, k++] = emline(logLam_temp, lines[j], sigma)
line_wave = [line_wave, lines]
line_names = [line_names, names]

;         -----[OIII]-----
lines = [4958.92, 5006.84]  ; air wavelengths
emission_lines[*, k++] = 0.33*emline(logLam_temp, lines[0], sigma) + emline(logLam_temp, lines[1], sigma)
line_wave = [line_wave, lines[1]]
line_names = [line_names, '[OIII]5007d']  ; single template for this doublet

;         -----[OI]-----
lines = [6300.30, 6363.67]  ; air wavelengths
emission_lines[*, k++] = emline(logLam_temp, lines[0], sigma) + 0.33*emline(logLam_temp, lines[1], sigma)
line_wave = [line_wave, lines[1]]
line_names = [line_names, '[OI]6300d']  ; single template for this doublet

;         -----[NII]-----
lines = [6548.03, 6583.41]  ; air wavelengths
emission_lines[*, k++] = 0.33*emline(logLam_temp, lines[0], sigma) + emline(logLam_temp, lines[1], sigma)
line_wave = [line_wave, lines[1]]
line_names = [line_names, '[NII]6583d']  ; single template for this doublet


;         -----[NeII]-----
lines = [3868.7]  ; air wavelengths
emission_lines[*, k++] = 0.33*emline(logLam_temp, lines[0], sigma)
line_wave = [line_wave, lines]
line_names = [line_names, '[NeII]']  ; single template for this doublet



; Only include lines falling within the estimated fitted wavelength range.
; This is important to avoid instabilities in the PPXF system solution
;
w = where((line_wave gt lamRange_gal[0]) and (line_wave lt lamRange_gal[1]))
emission_lines = emission_lines[*, w]
line_wave = line_wave[w]
line_names = line_names[w]

return, emission_lines
end
;------------------------------------------------------------------------------
