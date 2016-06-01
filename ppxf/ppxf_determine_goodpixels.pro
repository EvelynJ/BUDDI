;------------------------------------------------------------------------------
function ppxf_determine_goodPixels, logLam, lamRangeTemp, vel
;
; PPXF_DETERMINE_GOODPIXELS: Example routine to generate the vector of goodPixels 
;     to be used as input keyword for the routine PPXF. This is useful to mask 
;     gas emission lines or atmospheric absorptions. 
;     It can be trivially adapted to mask different lines.
; 
; INPUT PARAMETERS:
; - LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom 
;     of each pixel of the log rebinned *galaxy* spectrum.
; - LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the minimum and
;     maximum wavelength in Angstrom in the stellar *template* used in PPXF.
; - VEL: Estimate of the galaxy velocity in km/s.
; 
; V1.0: Michele Cappellari, Leiden, 9 September 2005
; V1.01: Made a separate routine and included additional common emission lines. 
;   MC, Oxford 12 January 2012
; V1.02: Included more lines. MC, Oxford, 7 Januray 2014

;        -----[OII]-----   Hdelta    Hgamma   Hbeta    -----[OIII]-----   [OI]    -----[NII]-----   Halpha   -----[SII]-----  
lines = [3726.03, 3728.82, 4101.76, 4340.47, 4861.33, 4958.92, 5006.84, 6300.30, 6548.03, 6583.41, 6562.80, 6716.47, 6730.85] 
dv = lines*0+800d ; width/2 of masked gas emission region in km/s
c = 299792.458d ; speed of light in km/s

flag = bytarr(n_elements(logLam))

for j=0,n_elements(lines)-1 do $
    flag or= logLam gt alog(lines[j]) + (vel - dv[j])/c $
         and logLam lt alog(lines[j]) + (vel + dv[j])/c

flag or= logLam lt alog(lamRangeTemp[0]) + (vel + 900d)/c ; Mask edges of
flag or= logLam gt alog(lamRangeTemp[1]) + (vel - 900d)/c ; stellar library

return, where(flag eq 0)
end
;------------------------------------------------------------------------------
