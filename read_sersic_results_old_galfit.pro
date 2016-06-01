FUNCTION read_sersic_results_old_galfit, obj, bd=bd
IF file_test(obj) THEN BEGIN
    hd = headfits(obj, exten = 2)
    mag0 = sxpar(hd, '2_MAG')
    mag = float(strmid(mag0, 0, strpos(mag0, '+/-')))
    magerr = float(strmid(mag0, strpos(mag0, '+/-')+3, strlen(mag0)))
    re0 = sxpar(hd, '2_RE')
    re = float(strmid(re0, 0, strpos(re0, '+/-')))
    reerr = float(strmid(re0, strpos(re0, '+/-')+3, strlen(re0)))
    n0 = sxpar(hd, '2_N')
    n = float(strmid(n0, 0, strpos(n0, '+/-')))
    nerr = float(strmid(n0, strpos(n0, '+/-')+3, strlen(n0)))
    q0 = sxpar(hd, '2_AR')
    q = float(strmid(q0, 0, strpos(q0, '+/-')))
    qerr = float(strmid(q0, strpos(q0, '+/-')+3, strlen(q0)))
    pa0 = sxpar(hd, '2_PA')
    pa = float(strmid(pa0, 0, strpos(pa0, '+/-')))
    paerr = float(strmid(pa0, strpos(pa0, '+/-')+3, strlen(pa0)))
    x0 = sxpar(hd, '2_XC')
    x = float(strmid(x0, 0, strpos(x0, '+/-')))
    xerr = float(strmid(x0, strpos(x0, '+/-')+3, strlen(x0)))
    y0 = sxpar(hd, '2_YC')
    y = float(strmid(y0, 0, strpos(y0, '+/-')))
    yerr = float(strmid(y0, strpos(y0, '+/-')+3, strlen(y0)))
    s0 = sxpar(hd, '1_SKY')
    sky = float(strmid(s0, 1, strpos(s0, ']')))
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
ENDELSE

if not keyword_set(bd) then BEGIN
    feedback = create_struct('mag_galfit', mag, 'magerr_galfit', magerr, $
                             're_galfit', re, 'reerr_galfit', reerr, $
                             'n_galfit', n, 'nerr_galfit', nerr, $
                             'q_galfit', q, 'qerr_galfit', qerr, $
                             'pa_galfit', pa, 'paerr_galfit', paerr, $
                             'x_galfit', x, 'xerr_galfit', xerr, $
                             'y_galfit', y, 'yerr_galfit', yerr, $
                             'psf_galfit', psf, 'sky_galfit', sky, $
                             'mag_galfit_band', mag, 'magerr_galfit_band', magerr, $
                             're_galfit_band', re, 'reerr_galfit_band', reerr, $
                             'n_galfit_band', n, 'nerr_galfit_band', nerr, $
                             'q_galfit_band', q, 'qerr_galfit_band', qerr, $
                             'pa_galfit_band', pa, 'paerr_galfit_band', paerr, $
                             'x_galfit_band', x, 'xerr_galfit_band', xerr, $
                             'y_galfit_band', y, 'yerr_galfit_band', yerr, $
                             'sky_galfit_band', sky, $
                             'mag_galfit_cheb', mag, 'magerr_galfit_cheb', magerr, $
                             're_galfit_cheb', re, 'reerr_galfit_cheb', reerr, $
                             'n_galfit_cheb', n, 'nerr_galfit_cheb', nerr, $
                             'q_galfit_cheb', q, 'qerr_galfit_cheb', qerr, $
                             'pa_galfit_cheb', pa, 'paerr_galfit_cheb', paerr, $
                             'x_galfit_cheb', x, 'xerr_galfit_cheb', xerr, $
                             'y_galfit_cheb', y, 'yerr_galfit_cheb', yerr, $
                             'sky_galfit_cheb', sky, $
                             'psf_galfit_band', psf, $
                             'chisq_galfit', chisq_galfit, $
                             'ndof_galfit', ndof_galfit, $
                             'nfree_galfit', nfree_galfit, $
                             'nfix_galfit', nfix_galfit, $
                             'chi2nu_galfit', chi2nu_galfit, $
                             'neigh_galfit', neigh_galfit, 'flag_galfit', flag_galfit, $
                             'X_GALFIT_DEG', -99, $
                             'Y_GALFIT_DEG', -99, $
                             'MAG_GALFIT_DEG', -99, $
                             'RE_GALFIT_DEG', -99, $
                             'N_GALFIT_DEG', -99, $
                             'Q_GALFIT_DEG', -99, $
                             'PA_GALFIT_DEG', -99)
endif

if keyword_set(bd) then begin
    IF file_test(obj) THEN BEGIN
        hd = headfits(obj, exten = 2)
        mag0_b = sxpar(hd, '3_MAG')
        mag_b = float(strmid(mag0_b, 0, strpos(mag0_b, '+/-')))
        magerr_b = float(strmid(mag0_b, strpos(mag0_b, '+/-')+3, strlen(mag0_b)))
        re0_b = sxpar(hd, '3_RE')
        re_b = float(strmid(re0_b, 0, strpos(re0_b, '+/-')))
        reerr_b = float(strmid(re0_b, strpos(re0_b, '+/-')+3, strlen(re0_b)))
        n0_b = sxpar(hd, '3_N')
        n_b = float(strmid(n0_b, 0, strpos(n0_b, '+/-')))
        nerr_b = float(strmid(n0_b, strpos(n0_b, '+/-')+3, strlen(n0_b)))
        q0_b = sxpar(hd, '3_AR')
        q_b = float(strmid(q0_b, 0, strpos(q0_b, '+/-')))
        qerr_b = float(strmid(q0_b, strpos(q0_b, '+/-')+3, strlen(q0_b)))
        pa0_b = sxpar(hd, '3_PA')
        pa_b = float(strmid(pa0_b, 0, strpos(pa0_b, '+/-')))
        paerr_b = float(strmid(pa0_b, strpos(pa0_b, '+/-')+3, strlen(pa0_b)))
        x0_b = sxpar(hd, '3_XC')
        x_b = float(strmid(x0_b, 0, strpos(x0_b, '+/-')))
        xerr_b = float(strmid(x0_b, strpos(x0_b, '+/-')+3, strlen(x0_b)))
        y0_b = sxpar(hd, '3_YC')
        y_b = float(strmid(y0_b, 0, strpos(y0_b, '+/-')))
        yerr_b = float(strmid(y0_b, strpos(y0_b, '+/-')+3, strlen(y0_b)))
; find number of neighbors
        neigh_galfit = comp-4
    ENDIF ELSE BEGIN
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
    
    feedback = create_struct('mag_galfit_d', mag, 'magerr_galfit_d', magerr, $
                             're_galfit_d', re, 'reerr_galfit_d', reerr, $
                             'n_galfit_d', n, 'nerr_galfit_d', nerr, $
                             'q_galfit_d', q, 'qerr_galfit_d', qerr, $
                             'pa_galfit_d', pa, 'paerr_galfit_d', paerr, $
                             'x_galfit_d', x, 'xerr_galfit_d', xerr, $
                             'y_galfit_d', y, 'yerr_galfit_d', yerr, $
                             'mag_galfit_b', mag_b, 'magerr_galfit_b', magerr_b, $
                             're_galfit_b', re_b, 'reerr_galfit_b', reerr_b, $
                             'n_galfit_b', n_b, 'nerr_galfit_b', nerr_b, $
                             'q_galfit_b', q_b, 'qerr_galfit_b', qerr_b, $
                             'pa_galfit_b', pa_b, 'paerr_galfit_b', paerr_b, $
                             'x_galfit_b', x_b, 'xerr_galfit_b', xerr_b, $
                             'y_galfit_b', y_b, 'yerr_galfit_b', yerr_b, $
                             'psf_galfit_bd', psf, 'sky_galfit_bd', sky, $
                             'mag_galfit_band_d', mag, 'magerr_galfit_band_d', magerr, $
                             're_galfit_band_d', re, 'reerr_galfit_band_d', reerr, $
                             'n_galfit_band_d', n, 'nerr_galfit_band_d', nerr, $
                             'q_galfit_band_d', q, 'qerr_galfit_band_d', qerr, $
                             'pa_galfit_band_d', pa, 'paerr_galfit_band_d', paerr, $
                             'x_galfit_band_d', x, 'xerr_galfit_band_d', xerr, $
                             'y_galfit_band_d', y, 'yerr_galfit_band_d', yerr, $
                             'mag_galfit_band_b', mag_b, 'magerr_galfit_band_b', magerr_b, $
                             're_galfit_band_b', re_b, 'reerr_galfit_band_b', reerr_b, $
                             'n_galfit_band_b', n_b, 'nerr_galfit_band_b', nerr_b, $
                             'q_galfit_band_b', q_b, 'qerr_galfit_band_b', qerr_b, $
                             'pa_galfit_band_b', pa_b, 'paerr_galfit_band_b', paerr_b, $
                             'x_galfit_band_b', x_b, 'xerr_galfit_band_b', xerr_b, $
                             'y_galfit_band_b', y_b, 'yerr_galfit_band_b', yerr_b, $
                             'sky_galfit_band_bd', sky, $
                             'mag_galfit_cheb_d', mag, 'magerr_galfit_cheb_d', magerr, $
                             're_galfit_cheb_d', re, 'reerr_galfit_cheb_d', reerr, $
                             'n_galfit_cheb_d', n, 'nerr_galfit_cheb_d', nerr, $
                             'q_galfit_cheb_d', q, 'qerr_galfit_cheb_d', qerr, $
                             'pa_galfit_cheb_d', pa, 'paerr_galfit_cheb_d', paerr, $
                             'x_galfit_cheb_d', x, 'xerr_galfit_cheb_d', xerr, $
                             'y_galfit_cheb_d', y, 'yerr_galfit_cheb_d', yerr, $
                             'mag_galfit_cheb_b', mag_b, 'magerr_galfit_cheb_b', magerr_b, $
                             're_galfit_cheb_b', re_b, 'reerr_galfit_cheb_b', reerr_b, $
                             'n_galfit_cheb_b', n_b, 'nerr_galfit_cheb_b', nerr_b, $
                             'q_galfit_cheb_b', q_b, 'qerr_galfit_cheb_b', qerr_b, $
                             'pa_galfit_cheb_b', pa_b, 'paerr_galfit_cheb_b', paerr_b, $
                             'x_galfit_cheb_b', x_b, 'xerr_galfit_cheb_b', xerr_b, $
                             'y_galfit_cheb_b', y_b, 'yerr_galfit_cheb_b', yerr_b, $
                             'sky_galfit_cheb_bd', sky, $
                             'psf_galfit_band_bd', psf, $
                             'chisq_galfit_bd', chisq_galfit, $
                             'ndof_galfit_bd', ndof_galfit, $
                             'nfree_galfit_bd', nfree_galfit, $
                             'nfix_galfit_bd', nfix_galfit, $
                             'chi2nu_galfit_bd', chi2nu_galfit, $
                             'neigh_galfit_bd', neigh_galfit, 'flag_galfit_bd', flag_galfit, $
                             'X_GALFIT_DEG_B', -99, $
                             'Y_GALFIT_DEG_B', -99, $
                             'MAG_GALFIT_DEG_B', -99, $
                             'RE_GALFIT_DEG_B', -99, $
                             'N_GALFIT_DEG_B', -99, $
                             'Q_GALFIT_DEG_B', -99, $
                             'PA_GALFIT_DEG_B', -99, $
                             'X_GALFIT_DEG_D', -99, $
                             'Y_GALFIT_DEG_D', -99, $
                             'MAG_GALFIT_DEG_D', -99, $
                             'RE_GALFIT_DEG_D', -99, $
                             'N_GALFIT_DEG_D', -99, $
                             'Q_GALFIT_DEG_D', -99, $
                             'PA_GALFIT_DEG_D', -99)
ENDIF       
return, feedback
END