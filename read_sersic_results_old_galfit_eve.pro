FUNCTION read_sersic_results_old_galfit_eve, obj, bd=bd
IF file_test(obj) THEN BEGIN
    hd = headfits(obj, exten = 2)
    s0 = sxpar(hd, '1_SKY')
    sky = float(strmid(s0, 1, strpos(s0, ']')))
    
    mag0 = sxpar(hd, '2_MAG')
    mag0 = strmid(mag0, 1,strlen(mag0)-1)
    mag = float(mag0)
    re0 = sxpar(hd, '2_RE')
    re0 = strmid(re0, 1,strlen(re0)-1)
    re = float(re0)
    n0 = sxpar(hd, '2_N')
    n0 = strmid(n0, 1,strlen(n0)-1)
    n = float(n0)
    q0 = sxpar(hd, '2_AR')
    q0 = strmid(q0, 1,strlen(q0)-1)
    q = float(q0)
    pa0 = sxpar(hd, '2_PA')
    pa0 = strmid(pa0, 1,strlen(pa0)-1)
    pa = float(pa0)
    x0 = sxpar(hd, '2_XC')
    x0 = strmid(x0, 1,strlen(x0)-1)
    x = float(x0)
    y0 = sxpar(hd, '2_YC')
    y0 = strmid(y0, 1,strlen(y0)-1)
    y = float(y0)
    
    hd = headfits(obj, exten = 2)
    mag0_b = sxpar(hd, '3_MAG')
    mag0_b = strmid(mag0_b, 1,strlen(mag0_b)-1)
    mag_b = float(mag0_b)
    re0_b = sxpar(hd, '3_RE')
    re0_b = strmid(re0_b, 1,strlen(re0_b)-1)
    re_b = float(re0_b)
    n0_b = sxpar(hd, '3_N')
    n0_b = strmid(n0_b, 1,strlen(n0_b)-1)
    n_b = float(n0_b)
    q0_b = sxpar(hd, '3_AR')
    q0_b = strmid(q0_b, 1,strlen(q0_b)-1)
    q_b = float(q0_b)
    pa0_b = sxpar(hd, '3_PA')
    pa0_b = strmid(pa0_b, 1,strlen(pa0_b)-1)
    pa_b = float(pa0_b)
    x0_b = sxpar(hd, '3_XC')
    x0_b = strmid(x0_b, 1,strlen(x0_b)-1)
    x_b = float(x0_b)
    y0_b = sxpar(hd, '3_YC')
    y0_b = strmid(y0_b, 1,strlen(y0_b)-1)
    y_b = float(y0_b)
    
    psf0 = sxpar(hd, 'PSF') 
    psf= strtrim(psf0, 2)
     
; find number of neighbors
    comp=0
    repeat comp = comp +1 until sxpar(hd, 'COMP_'+strtrim(comp,2)) eq '0'
    neigh_galfit = comp-3
    flag_galfit = 2
    chisq_galfit = float(strmid(sxpar(hd, 'CHISQ'),2)) 
    ndof_galfit = fix(strmid(sxpar(hd, 'NDOF'),2))
    nfree_galfit = fix(strmid(sxpar(hd, 'NFREE'),2))
    nfix_galfit = fix(strmid(sxpar(hd, 'NFIX'),2))
    chi2nu_galfit = float(strmid(sxpar(hd, 'CHI2NU'),2))
ENDIF ELSE BEGIN
    mag = -999.
    magerr = 99999.
    re = -99.
    reerr = 99999.
    n = -99.
    nerr = 99999.
    q = -99.
    qerr = 99999.
    pa = 0.
    paerr = 99999.
    x = 0.
    xerr = 99999.
    y = 0.
    yerr = 99999.
    sky = -999.
    neigh_galfit = -99
    flag_galfit = 1
    chisq_galfit = -99.
    ndof_galfit = -99
    nfree_galfit = -99
    nfix_galfit = -99
    chi2nu_galfit = -99.
    psf='none'

    mag_b = -999.
    magerr_b = 99999.
    re_b = -99.
    reerr_b = 99999.
    n_b = -99.
    nerr_b = 99999.
    q_b = -99.
    qerr_b = 99999.
    pa_b = 0.
    paerr_b = 99999.
    x_b = 0.
    xerr_b = 99999.
    y_b = 0.
    yerr_b = 99999.
ENDELSE
    
    feedback = create_struct('mag_galfit_d', mag, $
                             're_galfit_d', re, $
                             'n_galfit_d', n,  $
                             'q_galfit_d', q,  $
                             'pa_galfit_d', pa, $
                             'x_galfit_d', x,  $
                             'y_galfit_d', y,  $
                             'mag_galfit_b', mag_b,  $
                             're_galfit_b', re_b,  $
                             'n_galfit_b', n_b,  $
                             'q_galfit_b', q_b,  $
                             'pa_galfit_b', pa_b,  $
                             'x_galfit_b', x_b,  $
                             'y_galfit_b', y_b,  $
                             'psf_galfit_bd', psf,  $
                             'mag_galfit_band_d', mag,  $
                             're_galfit_band_d', re,  $
                             'n_galfit_band_d', n,  $
                             'q_galfit_band_d', q,  $
                             'pa_galfit_band_d', pa,  $
                             'x_galfit_band_d', x,  $
                             'y_galfit_band_d', y, $
                             'mag_galfit_band_b', mag_b,  $
                             're_galfit_band_b', re_b,  $
                             'n_galfit_band_b', n_b,  $
                             'q_galfit_band_b', q_b, $
                             'pa_galfit_band_b', pa_b,  $
                             'x_galfit_band_b', x_b,  $
                             'y_galfit_band_b', y_b, $
                             'sky_galfit_band_bd', sky, $
                             'mag_galfit_cheb_d', mag,  $
                             're_galfit_cheb_d', re, $
                             'n_galfit_cheb_d', n,  $
                             'q_galfit_cheb_d', q,  $
                             'pa_galfit_cheb_d', pa,  $
                             'x_galfit_cheb_d', x,  $
                             'y_galfit_cheb_d', y,  $
                             'mag_galfit_cheb_b', mag_b,  $
                             're_galfit_cheb_b', re_b, $
                             'n_galfit_cheb_b', n_b, $
                             'q_galfit_cheb_b', q_b,  $
                             'pa_galfit_cheb_b', pa_b,  $
                             'x_galfit_cheb_b', x_b, $
                             'y_galfit_cheb_b', y_b,  $
                             'sky_galfit_cheb_bd', sky, $
                             'psf_galfit_band_bd', psf, $
                             'chisq_galfit_bd', chisq_galfit, $
                             'ndof_galfit_bd', ndof_galfit, $
                             'nfree_galfit_bd', nfree_galfit, $
                             'nfix_galfit_bd', nfix_galfit, $
                             'chi2nu_galfit_bd', chi2nu_galfit)
       
return, feedback
END