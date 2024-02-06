FUNCTION read_sersic_results_2comp, obj, nband, bd=bd,boxy=boxy
IF file_test(obj[0]) THEN BEGIN
    result = mrdfits(obj[0], 'FINAL_BAND',/silent)
    res_cheb = mrdfits(obj[0], 'FINAL_CHEB',/silent)
    fit_info = mrdfits(obj[0], 'FIT_INFO',/silent)
    band_info = mrdfits(obj[0], 'BAND_INFO',/silent)
;       hd = headfits(obj[0], exten = nband+1)
    comp=1
    repeat comp = comp +1 until tag_exist(result, 'COMP'+strtrim(comp,2)+'_MAG') eq 0
; delete feedback, just in case the format of one is different,
; avoiding crash
    delvarx, feedback
    if not keyword_set(bd) then begin
        feedback = create_struct('mag_galfit_d', result[0].COMP2_MAG, 'magerr_galfit_d',result[0].COMP2_MAG_ERR, $
                                 're_galfit_d', result[0].COMP2_RE, 'reerr_galfit_d', result[0].COMP2_RE_ERR, $
                                 'n_galfit_d', result[0].COMP2_N, 'nerr_galfit_d' ,result[0].COMP2_N_ERR, $
                                 'q_galfit_d', result[0].COMP2_AR, 'qerr_galfit_d', result[0].COMP2_AR_ERR, $
                                 'pa_galfit_d', result[0].COMP2_PA, 'paerr_galfit_d', result[0].COMP2_PA_ERR, $
                                 'x_galfit_d', result[0].COMP2_XC, 'xerr_galfit_d', result[0].COMP2_XC_ERR, $
                                 'y_galfit_d', result[0].COMP2_YC, 'yerr_galfit_d', result[0].COMP2_YC_ERR, $
                                 'psf_galfit_d', strtrim(band_info[0].psf), 'sky_galfit', result[0].COMP1_SKY, $
                                 'mag_galfit_band_d', result.COMP2_MAG, 'magerr_galfit_band_d',result.COMP2_MAG_ERR, $
                                 're_galfit_band_d', result.COMP2_RE, 'reerr_galfit_band_d', result.COMP2_RE_ERR, $
                                 'n_galfit_band_d', result.COMP2_N, 'nerr_galfit_band_d' ,result.COMP2_N_ERR, $
                                 'q_galfit_band_d', result.COMP2_AR, 'qerr_galfit_band_d', result.COMP2_AR_ERR, $
                                 'pa_galfit_band_d', result.COMP2_PA, 'paerr_galfit_band_d', result.COMP2_PA_ERR, $
                                 'x_galfit_band_d', result.COMP2_XC, 'xerr_galfit_band_d', result.COMP2_XC_ERR, $
                                 'y_galfit_band_d', result.COMP2_YC, 'yerr_galfit_band_d', result.COMP2_YC_ERR, $
                                 'sky_galfit_band', result.COMP1_SKY, $
                                 'mag_galfit_cheb_d', res_cheb.COMP2_MAG, 'magerr_galfit_cheb_d',res_cheb.COMP2_MAG_ERR, $
                                 're_galfit_cheb_d', res_cheb.COMP2_RE, 'reerr_galfit_cheb_d', res_cheb.COMP2_RE_ERR, $
                                 'n_galfit_cheb_d', res_cheb.COMP2_N, 'nerr_galfit_cheb_d' ,res_cheb.COMP2_N_ERR, $
                                 'q_galfit_cheb_d', res_cheb.COMP2_AR, 'qerr_galfit_cheb_d', res_cheb.COMP2_AR_ERR, $
                                 'pa_galfit_cheb_d', res_cheb.COMP2_PA, 'paerr_galfit_cheb_d', res_cheb.COMP2_PA_ERR, $
                                 'x_galfit_cheb_d', res_cheb.COMP2_XC, 'xerr_galfit_cheb_d', res_cheb.COMP2_XC_ERR, $
                                 'y_galfit_cheb_d', res_cheb.COMP2_YC, 'yerr_galfit_cheb_d', res_cheb.COMP2_YC_ERR, $
;                                 'x_galfit_band_comp3', result.COMP3_XC, 'xerr_galfit_band_comp3',result.COMP3_XC_ERR, $
;                                 'y_galfit_band_comp3', result.COMP3_YC, 'yerr_galfit_band_comp3',result.COMP3_YC_ERR, $
;                                 'mag_galfit_band_comp3', result.COMP3_MAG, 'magerr_galfit_band_comp3',result.COMP3_MAG_ERR, $
;                                 'x_galfit_cheb_comp3', res_cheb.COMP3_XC, 'xerr_galfit_cheb_comp3',res_cheb.COMP3_XC_ERR, $
;                                 'y_galfit_cheb_comp3', res_cheb.COMP3_YC, 'yerr_galfit_cheb_comp3',res_cheb.COMP3_YC_ERR, $
;                                 'mag_galfit_cheb_comp3', res_cheb.COMP3_MAG, 'magerr_galfit_cheb_comp3',res_cheb.COMP3_MAG_ERR, $
                                 'sky_galfit_cheb', res_cheb.COMP1_SKY, $
                                 'initfile', strtrim(fit_info.initfile,2), $
                                 'constrnt', strtrim(fit_info.constrnt,2), $
                                 'fitsect', strtrim(fit_info.fitsect,2), $
                                 'convbox', strtrim(fit_info.convbox,2), $
                                 'psf_galfit_band', strtrim(band_info.psf, 2), $
                                 'chisq_galfit', fit_info.chisq, $
                                 'ndof_galfit', fit_info.ndof, $
                                 'nfree_galfit', fit_info.nfree, $
                                 'nfix_galfit', fit_info.nfix, $
                                 'cputime_setup_galfit', fit_info.cputime_setup, $
                                 'cputime_fit_galfit', fit_info.cputime_fit, $
                                 'cputime_total_galfit', fit_info.cputime_total, $
                                 'chi2nu_galfit', fit_info.chi2nu, $
                                 'niter_galfit', fit_info.niter, $
                                 'version_galfit', fit_info.version, $
                                 'firstcon_galfit', fit_info.firstcon, $
                                 'lastcon_galfit', fit_info.lastcon, $
                                 'neigh_galfit', comp-3, 'flag_galfit', 2)
; TO BE ADDED:
; fitting time
; NEIGH_GALFIT HAS TO BE ADAPTED! WHY??
    ENDIF
    if keyword_set(bd) and keyword_set(boxy) then begin
        feedback = create_struct('mag_galfit_d', result[0].COMP2_MAG, 'magerr_galfit_d',result[0].COMP2_MAG_ERR, $
                                 're_galfit_d', result[0].COMP2_RE, 'reerr_galfit_d', result[0].COMP2_RE_ERR, $
                                 'n_galfit_d', result[0].COMP2_N, 'nerr_galfit_d' ,result[0].COMP2_N_ERR, $
                                 'q_galfit_d', result[0].COMP2_AR, 'qerr_galfit_d', result[0].COMP2_AR_ERR, $
                                 'pa_galfit_d', result[0].COMP2_PA, 'paerr_galfit_d', result[0].COMP2_PA_ERR, $
                                 'x_galfit_d', result[0].COMP2_XC, 'xerr_galfit_d', result[0].COMP2_XC_ERR, $
                                 'y_galfit_d', result[0].COMP2_YC, 'yerr_galfit_d', result[0].COMP2_YC_ERR, $
                                 'mag_galfit_b', result[0].COMP3_MAG, 'magerr_galfit_b',result[0].COMP3_MAG_ERR, $
                                 're_galfit_b', result[0].COMP3_RE, 'reerr_galfit_b', result[0].COMP3_RE_ERR, $
                                 'n_galfit_b', result[0].COMP3_N, 'nerr_galfit_b' ,result[0].COMP3_N_ERR, $
                                 'q_galfit_b', result[0].COMP3_AR, 'qerr_galfit_b', result[0].COMP3_AR_ERR, $
                                 'pa_galfit_b', result[0].COMP3_PA, 'paerr_galfit_b', result[0].COMP3_PA_ERR, $
                                 'boxy_galfit_b', result[0].COMP3_C0, 'boxyerr_galfit_b', result[0].COMP3_C0_ERR, $
                                 'x_galfit_b', result[0].COMP3_XC, 'xerr_galfit_b', result[0].COMP3_XC_ERR, $
                                 'y_galfit_b', result[0].COMP3_YC, 'yerr_galfit_b', result[0].COMP3_YC_ERR, $
                                 'psf_galfit', strtrim(band_info[0].psf,2), 'sky_galfit', result[0].COMP1_SKY, $
                                 'mag_galfit_band_d', result.COMP2_MAG, 'magerr_galfit_band_d',result.COMP2_MAG_ERR, $
                                 're_galfit_band_d', result.COMP2_RE, 'reerr_galfit_band_d', result.COMP2_RE_ERR, $
                                 'n_galfit_band_d', result.COMP2_N, 'nerr_galfit_band_d' ,result.COMP2_N_ERR, $
                                 'q_galfit_band_d', result.COMP2_AR, 'qerr_galfit_band_d', result.COMP2_AR_ERR, $
                                 'pa_galfit_band_d', result.COMP2_PA, 'paerr_galfit_band_d', result.COMP2_PA_ERR, $
                                 'x_galfit_band_d', result.COMP2_XC, 'xerr_galfit_band_d', result.COMP2_XC_ERR, $
                                 'y_galfit_band_d', result.COMP2_YC, 'yerr_galfit_band_d', result.COMP2_YC_ERR, $
                                 'mag_galfit_band_b', result.COMP3_MAG, 'magerr_galfit_band_b',result.COMP3_MAG_ERR, $
                                 're_galfit_band_b', result.COMP3_RE, 'reerr_galfit_band_b', result.COMP3_RE_ERR, $
                                 'n_galfit_band_b', result.COMP3_N, 'nerr_galfit_band_b' ,result.COMP3_N_ERR, $
                                 'q_galfit_band_b', result.COMP3_AR, 'qerr_galfit_band_b', result.COMP3_AR_ERR, $
                                 'pa_galfit_band_b', result.COMP3_PA, 'paerr_galfit_band_b', result.COMP3_PA_ERR, $
                                 'boxy_galfit_band_b', result.COMP3_C0, 'boxyerr_galfit_band_b', result.COMP3_C0_ERR, $
                                 'x_galfit_band_b', result.COMP3_XC, 'xerr_galfit_band_b', result.COMP3_XC_ERR, $
                                 'y_galfit_band_b', result.COMP3_YC, 'yerr_galfit_band_b', result.COMP3_YC_ERR, $
                                 'sky_galfit_band', result.COMP1_SKY, $
;                                 'mag_galfit_band_psf', result.COMP4_MAG, 'magerr_galfit_band_psf',result.COMP3_MAG_ERR, $
;                                 'x_galfit_band_comp3', result.COMP4_XC, 'xerr_galfit_band_comp3',result.COMP4_XC_ERR, $
;                                 'y_galfit_band_comp3', result.COMP4_YC, 'yerr_galfit_band_comp3',result.COMP4_YC_ERR, $
;                                 'mag_galfit_band_comp3', result.COMP4_MAG, 'magerr_galfit_band_comp3',result.COMP4_MAG_ERR, $
                                 'mag_galfit_cheb_d', res_cheb.COMP2_MAG, 'magerr_galfit_cheb_d',res_cheb.COMP2_MAG_ERR, $
                                 're_galfit_cheb_d', res_cheb.COMP2_RE, 'reerr_galfit_cheb_d', res_cheb.COMP2_RE_ERR, $
                                 'n_galfit_cheb_d', res_cheb.COMP2_N, 'nerr_galfit_cheb_d' ,res_cheb.COMP2_N_ERR, $
                                 'q_galfit_cheb_d', res_cheb.COMP2_AR, 'qerr_galfit_cheb_d', res_cheb.COMP2_AR_ERR, $
                                 'pa_galfit_cheb_d', res_cheb.COMP2_PA, 'paerr_galfit_cheb_d', res_cheb.COMP2_PA_ERR, $
                                 'x_galfit_cheb_d', res_cheb.COMP2_XC, 'xerr_galfit_cheb_d', res_cheb.COMP2_XC_ERR, $
                                 'y_galfit_cheb_d', res_cheb.COMP2_YC, 'yerr_galfit_cheb_d', res_cheb.COMP2_YC_ERR, $
                                 'mag_galfit_cheb_b', res_cheb.COMP3_MAG, 'magerr_galfit_cheb_b',res_cheb.COMP3_MAG_ERR, $
                                 're_galfit_cheb_b', res_cheb.COMP3_RE, 'reerr_galfit_cheb_b', res_cheb.COMP3_RE_ERR, $
                                 'n_galfit_cheb_b', res_cheb.COMP3_N, 'nerr_galfit_cheb_b' ,res_cheb.COMP3_N_ERR, $
                                 'q_galfit_cheb_b', res_cheb.COMP3_AR, 'qerr_galfit_cheb_b', res_cheb.COMP3_AR_ERR, $
                                 'pa_galfit_cheb_b', res_cheb.COMP3_PA, 'paerr_galfit_cheb_b', res_cheb.COMP3_PA_ERR, $
                                 'x_galfit_cheb_b', res_cheb.COMP3_XC, 'xerr_galfit_cheb_b', res_cheb.COMP3_XC_ERR, $
                                 'y_galfit_cheb_b', res_cheb.COMP3_YC, 'yerr_galfit_cheb_b', res_cheb.COMP3_YC_ERR, $
;                                 'x_galfit_cheb_comp3', res_cheb.COMP4_XC, 'xerr_galfit_cheb_comp3',res_cheb.COMP4_XC_ERR, $
;                                 'y_galfit_cheb_comp3', res_cheb.COMP4_YC, 'yerr_galfit_cheb_comp3',res_cheb.COMP4_YC_ERR, $
;                                 'mag_galfit_cheb_comp3', res_cheb.COMP4_MAG, 'magerr_galfit_cheb_comp3',res_cheb.COMP4_MAG_ERR, $
                                 'sky_galfit_cheb', res_cheb.COMP1_SKY, $
                                 'initfile_bd', strtrim(fit_info.initfile,2), $
                                 'constrnt_bd', strtrim(fit_info.constrnt,2), $
                                 'psf_galfit_band', strtrim(band_info.psf, 2), $
                                 'chisq_galfit_bd', fit_info.chisq, $
                                 'ndof_galfit_bd', fit_info.ndof, $
                                 'nfree_galfit_bd', fit_info.nfree, $
                                 'nfix_galfit_bd', fit_info.nfix, $
                                 'cputime_setup_galfit_bd', fit_info.cputime_setup, $
                                 'cputime_fit_galfit_bd', fit_info.cputime_fit, $
                                 'cputime_total_galfit_bd', fit_info.cputime_total, $
                                 'chi2nu_galfit_bd', fit_info.chi2nu, $
                                 'niter_galfit_bd', fit_info.niter, $
                                 'version_galfit_bd', fit_info.version, $
                                 'firstcon_galfit_bd', fit_info.firstcon, $
                                 'lastcon_galfit_bd', fit_info.lastcon, $
                                 'neigh_galfit_bd', comp-4, 'flag_galfit_bd', 2)
                                 
    ENDIF
    if keyword_set(bd) and ~keyword_set(boxy) then begin
      feedback = create_struct('mag_galfit_d', result[0].COMP2_MAG, 'magerr_galfit_d',result[0].COMP2_MAG_ERR, $
        're_galfit_d', result[0].COMP2_RE, 'reerr_galfit_d', result[0].COMP2_RE_ERR, $
        'n_galfit_d', result[0].COMP2_N, 'nerr_galfit_d' ,result[0].COMP2_N_ERR, $
        'q_galfit_d', result[0].COMP2_AR, 'qerr_galfit_d', result[0].COMP2_AR_ERR, $
        'pa_galfit_d', result[0].COMP2_PA, 'paerr_galfit_d', result[0].COMP2_PA_ERR, $
        'x_galfit_d', result[0].COMP2_XC, 'xerr_galfit_d', result[0].COMP2_XC_ERR, $
        'y_galfit_d', result[0].COMP2_YC, 'yerr_galfit_d', result[0].COMP2_YC_ERR, $
        'mag_galfit_b', result[0].COMP3_MAG, 'magerr_galfit_b',result[0].COMP3_MAG_ERR, $
        're_galfit_b', result[0].COMP3_RE, 'reerr_galfit_b', result[0].COMP3_RE_ERR, $
        'n_galfit_b', result[0].COMP3_N, 'nerr_galfit_b' ,result[0].COMP3_N_ERR, $
        'q_galfit_b', result[0].COMP3_AR, 'qerr_galfit_b', result[0].COMP3_AR_ERR, $
        'pa_galfit_b', result[0].COMP3_PA, 'paerr_galfit_b', result[0].COMP3_PA_ERR, $
        'x_galfit_b', result[0].COMP3_XC, 'xerr_galfit_b', result[0].COMP3_XC_ERR, $
        'y_galfit_b', result[0].COMP3_YC, 'yerr_galfit_b', result[0].COMP3_YC_ERR, $
        'psf_galfit', strtrim(band_info[0].psf,2), 'sky_galfit', result[0].COMP1_SKY, $
        'mag_galfit_band_d', result.COMP2_MAG, 'magerr_galfit_band_d',result.COMP2_MAG_ERR, $
        're_galfit_band_d', result.COMP2_RE, 'reerr_galfit_band_d', result.COMP2_RE_ERR, $
        'n_galfit_band_d', result.COMP2_N, 'nerr_galfit_band_d' ,result.COMP2_N_ERR, $
        'q_galfit_band_d', result.COMP2_AR, 'qerr_galfit_band_d', result.COMP2_AR_ERR, $
        'pa_galfit_band_d', result.COMP2_PA, 'paerr_galfit_band_d', result.COMP2_PA_ERR, $
        'x_galfit_band_d', result.COMP2_XC, 'xerr_galfit_band_d', result.COMP2_XC_ERR, $
        'y_galfit_band_d', result.COMP2_YC, 'yerr_galfit_band_d', result.COMP2_YC_ERR, $
        'mag_galfit_band_b', result.COMP3_MAG, 'magerr_galfit_band_b',result.COMP3_MAG_ERR, $
        're_galfit_band_b', result.COMP3_RE, 'reerr_galfit_band_b', result.COMP3_RE_ERR, $
        'n_galfit_band_b', result.COMP3_N, 'nerr_galfit_band_b' ,result.COMP3_N_ERR, $
        'q_galfit_band_b', result.COMP3_AR, 'qerr_galfit_band_b', result.COMP3_AR_ERR, $
        'pa_galfit_band_b', result.COMP3_PA, 'paerr_galfit_band_b', result.COMP3_PA_ERR, $
        'x_galfit_band_b', result.COMP3_XC, 'xerr_galfit_band_b', result.COMP3_XC_ERR, $
        'y_galfit_band_b', result.COMP3_YC, 'yerr_galfit_band_b', result.COMP3_YC_ERR, $
        'sky_galfit_band', result.COMP1_SKY, $
        ;                                 'mag_galfit_band_psf', result.COMP4_MAG, 'magerr_galfit_band_psf',result.COMP3_MAG_ERR, $
        ;                                 'x_galfit_band_comp3', result.COMP4_XC, 'xerr_galfit_band_comp3',result.COMP4_XC_ERR, $
        ;                                 'y_galfit_band_comp3', result.COMP4_YC, 'yerr_galfit_band_comp3',result.COMP4_YC_ERR, $
        ;                                 'mag_galfit_band_comp3', result.COMP4_MAG, 'magerr_galfit_band_comp3',result.COMP4_MAG_ERR, $
        'mag_galfit_cheb_d', res_cheb.COMP2_MAG, 'magerr_galfit_cheb_d',res_cheb.COMP2_MAG_ERR, $
        're_galfit_cheb_d', res_cheb.COMP2_RE, 'reerr_galfit_cheb_d', res_cheb.COMP2_RE_ERR, $
        'n_galfit_cheb_d', res_cheb.COMP2_N, 'nerr_galfit_cheb_d' ,res_cheb.COMP2_N_ERR, $
        'q_galfit_cheb_d', res_cheb.COMP2_AR, 'qerr_galfit_cheb_d', res_cheb.COMP2_AR_ERR, $
        'pa_galfit_cheb_d', res_cheb.COMP2_PA, 'paerr_galfit_cheb_d', res_cheb.COMP2_PA_ERR, $
        'x_galfit_cheb_d', res_cheb.COMP2_XC, 'xerr_galfit_cheb_d', res_cheb.COMP2_XC_ERR, $
        'y_galfit_cheb_d', res_cheb.COMP2_YC, 'yerr_galfit_cheb_d', res_cheb.COMP2_YC_ERR, $
        'mag_galfit_cheb_b', res_cheb.COMP3_MAG, 'magerr_galfit_cheb_b',res_cheb.COMP3_MAG_ERR, $
        're_galfit_cheb_b', res_cheb.COMP3_RE, 'reerr_galfit_cheb_b', res_cheb.COMP3_RE_ERR, $
        'n_galfit_cheb_b', res_cheb.COMP3_N, 'nerr_galfit_cheb_b' ,res_cheb.COMP3_N_ERR, $
        'q_galfit_cheb_b', res_cheb.COMP3_AR, 'qerr_galfit_cheb_b', res_cheb.COMP3_AR_ERR, $
        'pa_galfit_cheb_b', res_cheb.COMP3_PA, 'paerr_galfit_cheb_b', res_cheb.COMP3_PA_ERR, $
        'x_galfit_cheb_b', res_cheb.COMP3_XC, 'xerr_galfit_cheb_b', res_cheb.COMP3_XC_ERR, $
        'y_galfit_cheb_b', res_cheb.COMP3_YC, 'yerr_galfit_cheb_b', res_cheb.COMP3_YC_ERR, $
        ;                                 'x_galfit_cheb_comp3', res_cheb.COMP4_XC, 'xerr_galfit_cheb_comp3',res_cheb.COMP4_XC_ERR, $
        ;                                 'y_galfit_cheb_comp3', res_cheb.COMP4_YC, 'yerr_galfit_cheb_comp3',res_cheb.COMP4_YC_ERR, $
        ;                                 'mag_galfit_cheb_comp3', res_cheb.COMP4_MAG, 'magerr_galfit_cheb_comp3',res_cheb.COMP4_MAG_ERR, $
        'sky_galfit_cheb', res_cheb.COMP1_SKY, $
        'initfile_bd', strtrim(fit_info.initfile,2), $
        'constrnt_bd', strtrim(fit_info.constrnt,2), $
        'psf_galfit_band', strtrim(band_info.psf, 2), $
        'chisq_galfit_bd', fit_info.chisq, $
        'ndof_galfit_bd', fit_info.ndof, $
        'nfree_galfit_bd', fit_info.nfree, $
        'nfix_galfit_bd', fit_info.nfix, $
        'cputime_setup_galfit_bd', fit_info.cputime_setup, $
        'cputime_fit_galfit_bd', fit_info.cputime_fit, $
        'cputime_total_galfit_bd', fit_info.cputime_total, $
        'chi2nu_galfit_bd', fit_info.chi2nu, $
        'niter_galfit_bd', fit_info.niter, $
        'version_galfit_bd', fit_info.version, $
        'firstcon_galfit_bd', fit_info.firstcon, $
        'lastcon_galfit_bd', fit_info.lastcon, $
        'neigh_galfit_bd', comp-4, 'flag_galfit_bd', 2)

    ENDIF
; to include:
; there is more band_info which is not used yet (band, wl, datain,
; sigma, MASL, magzpt) Not sure we'll need them!

ENDIF ELSE BEGIN
    psf=strarr(nband)
    for n=0,nband-1 do psf[n]='none'
    
    if not keyword_set(bd) then begin
        feedback = create_struct('mag_galfit_d', -999., 'magerr_galfit_d',99999., $
                                 're_galfit_d', -99., 'reerr_galfit_d', 99999., $
                                 'n_galfit_d', -99., 'nerr_galfit_d' ,99999., $
                                 'q_galfit_d', -99., 'qerr_galfit_d', 99999., $
                                 'pa_galfit_d', 0., 'paerr_galfit_d', 99999., $
                                 'x_galfit_d', 0., 'xerr_galfit_d', 99999., $
                                 'y_galfit_d', 0., 'yerr_galfit_d', 99999., $
                                 'psf_galfit', 'none', 'sky_galfit', -999., $
                                 'mag_galfit_band_d', fltarr(nband)-999., 'magerr_galfit_band_d',fltarr(nband)+99999., $
                                 're_galfit_band_d', fltarr(nband)-99., 'reerr_galfit_band_d', fltarr(nband)+99999., $
                                 'n_galfit_band_d', fltarr(nband)-99., 'nerr_galfit_band_d' ,fltarr(nband)+99999., $
                                 'q_galfit_band_d', fltarr(nband)-99., 'qerr_galfit_band_d', fltarr(nband)+99999., $
                                 'pa_galfit_band_d', fltarr(nband), 'paerr_galfit_band_d', fltarr(nband)+99999., $
                                 'x_galfit_band_d', fltarr(nband), 'xerr_galfit_band_d', fltarr(nband)+99999., $
                                 'y_galfit_band_d', fltarr(nband), 'yerr_galfit_band_d', fltarr(nband)+99999., $ 
                                 'sky_galfit_band', fltarr(nband)-999., $
                                 'mag_galfit_cheb_d', fltarr(nband)-999., 'magerr_galfit_cheb_d',fltarr(nband)+99999., $
                                 're_galfit_cheb_d', fltarr(nband)-99., 'reerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'n_galfit_cheb_d', fltarr(nband)-99., 'nerr_galfit_cheb_d' ,fltarr(nband)+99999., $
                                 'q_galfit_cheb_d', fltarr(nband)-99., 'qerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'pa_galfit_cheb_d', fltarr(nband), 'paerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'x_galfit_cheb_d', fltarr(nband), 'xerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'y_galfit_cheb_d', fltarr(nband), 'yerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'sky_galfit_cheb', fltarr(nband)-999., $
                                 'initfile_bd', ' ', $
                                 'constrnt_bd', ' ', $
                                 'psf_galfit_band', psf, $
                                 'chisq_galfit', -99., $
                                 'ndof_galfit', -99l, $
                                 'nfree_galfit', -99l, $
                                 'nfix_galfit', -99l, $
                                 'cputime_setup_galfit', -99., $
                                 'cputime_fit_galfit', -99., $
                                 'cputime_total_galfit', -99., $
                                 'chi2nu_galfit', -99., $
                                 'niter_galfit', -99, $
                                 'version_galfit', 'crash', $
                                 'firstcon_galfit', -99, $
                                 'lastcon_galfit', -99, $
                                 'neigh_galfit', -99, 'flag_galfit', 1)
                                ; TO BE ADDED:
; fitting time
; NEIGH_GALFIT HAS TO BE ADAPTED!
    ENDIF
    if keyword_set(bd) then begin
        feedback = create_struct('mag_galfit_d', -999., 'magerr_galfit_d',99999., $
                                 're_galfit_d', -99., 'reerr_galfit_d', 99999., $
                                 'n_galfit_d', -99., 'nerr_galfit_d' ,99999., $
                                 'q_galfit_d', -99., 'qerr_galfit_d', 99999., $
                                 'pa_galfit_d', 0., 'paerr_galfit_d', 99999., $
                                 'x_galfit_d', 0., 'xerr_galfit_d', 99999., $
                                 'y_galfit_d', 0., 'yerr_galfit_d', 99999., $
                                 'mag_galfit_b', -999., 'magerr_galfit_b',99999., $
                                 're_galfit_b', -99., 'reerr_galfit_b', 99999., $
                                 'n_galfit_b', -99., 'nerr_galfit_b' ,99999., $
                                 'q_galfit_b', -99., 'qerr_galfit_b', 99999., $
                                 'pa_galfit_b', 0., 'paerr_galfit_b', 99999., $
                                 'x_galfit_b', 0., 'xerr_galfit_b', 99999., $
                                 'y_galfit_b', 0., 'yerr_galfit_b', 99999., $
                                 'psf_galfit', 'none', 'sky_galfit', -999., $
                                 'mag_galfit_band_d', fltarr(nband)-999., 'magerr_galfit_band_d',fltarr(nband)+99999., $
                                 're_galfit_band_d', fltarr(nband)-99., 'reerr_galfit_band_d', fltarr(nband)+99999., $
                                 'n_galfit_band_d', fltarr(nband)-99., 'nerr_galfit_band_d' ,fltarr(nband)+99999., $
                                 'q_galfit_band_d', fltarr(nband)-99., 'qerr_galfit_band_d', fltarr(nband)+99999., $
                                 'pa_galfit_band_d', fltarr(nband), 'paerr_galfit_band_d', fltarr(nband)+99999., $
                                 'x_galfit_band_d', fltarr(nband), 'xerr_galfit_band_d', fltarr(nband)+99999., $
                                 'y_galfit_band_d', fltarr(nband), 'yerr_galfit_band_d', fltarr(nband)+99999., $ 
                                 'mag_galfit_band_b', fltarr(nband)-999., 'magerr_galfit_band_b',fltarr(nband)+99999., $
                                 're_galfit_band_b', fltarr(nband)-99., 'reerr_galfit_band_b', fltarr(nband)+99999., $
                                 'n_galfit_band_b', fltarr(nband)-99., 'nerr_galfit_band_b' ,fltarr(nband)+99999., $
                                 'q_galfit_band_b', fltarr(nband)-99., 'qerr_galfit_band_b', fltarr(nband)+99999., $
                                 'pa_galfit_band_b', fltarr(nband), 'paerr_galfit_band_b', fltarr(nband)+99999., $
                                 'x_galfit_band_b', fltarr(nband), 'xerr_galfit_band_b', fltarr(nband)+99999., $
                                 'y_galfit_band_b', fltarr(nband), 'yerr_galfit_band_b', fltarr(nband)+99999., $ 
                                 'sky_galfit_band', fltarr(nband)-999., $
;                                 'x_galfit_band_comp3', fltarr(nband), 'xerr_galfit_band_comp3', fltarr(nband)+99999., $
;                                 'y_galfit_band_comp3', fltarr(nband), 'yerr_galfit_band_comp3', fltarr(nband)+99999., $ 
;                                 'mag_galfit_band_comp3', fltarr(nband), 'magerr_galfit_band_comp3', fltarr(nband)+99999., $ 
                                 'mag_galfit_cheb_d', fltarr(nband)-999., 'magerr_galfit_cheb_d',fltarr(nband)+99999., $
                                 're_galfit_cheb_d', fltarr(nband)-99., 'reerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'n_galfit_cheb_d', fltarr(nband)-99., 'nerr_galfit_cheb_d' ,fltarr(nband)+99999., $
                                 'q_galfit_cheb_d', fltarr(nband)-99., 'qerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'pa_galfit_cheb_d', fltarr(nband), 'paerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'x_galfit_cheb_d', fltarr(nband), 'xerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'y_galfit_cheb_d', fltarr(nband), 'yerr_galfit_cheb_d', fltarr(nband)+99999., $
                                 'mag_galfit_cheb_b', fltarr(nband)-999., 'magerr_galfit_cheb_b',fltarr(nband)+99999., $
                                 're_galfit_cheb_b', fltarr(nband)-99., 'reerr_galfit_cheb_b', fltarr(nband)+99999., $
                                 'n_galfit_cheb_b', fltarr(nband)-99., 'nerr_galfit_cheb_b' ,fltarr(nband)+99999., $
                                 'q_galfit_cheb_b', fltarr(nband)-99., 'qerr_galfit_cheb_b', fltarr(nband)+99999., $
                                 'pa_galfit_cheb_b', fltarr(nband), 'paerr_galfit_cheb_b', fltarr(nband)+99999., $
                                 'x_galfit_cheb_b', fltarr(nband), 'xerr_galfit_cheb_b', fltarr(nband)+99999., $
                                 'y_galfit_cheb_b', fltarr(nband), 'yerr_galfit_cheb_b', fltarr(nband)+99999., $
                                 'sky_galfit_cheb', fltarr(nband)-999., $
                                 'initfile_bd', ' ', $
                                 'constrnt_bd', ' ', $
                                 'psf_galfit_band', psf, $
                                 'chisq_galfit_bd', -99., $
                                 'ndof_galfit_bd', -99l, $
                                 'nfree_galfit_bd', -99l, $
                                 'nfix_galfit_bd', -99l, $
                                 'cputime_setup_galfit_bd', -99., $
                                 'cputime_fit_galfit_bd', -99., $
                                 'cputime_total_galfit_bd', -99., $
                                 'chi2nu_galfit_bd', -99., $
                                 'niter_galfit_bd', -99, $
                                 'version_galfit_bd', 'crash', $
                                 'firstcon_galfit_bd', -99, $
                                 'lastcon_galfit_bd', -99, $
                                 'neigh_galfit_bd', -99, 'flag_galfit_bd', 1)
                                ; TO BE ADDED:
; fitting time
; NEIGH_GALFIT HAS TO BE ADAPTED!
    ENDIF
    
ENDELSE
return, feedback
;[mag, magerr, re, reerr, n, nerr, q, qerr, pa, paerr, $
;            x, xerr, y, yerr, sky, neigh_galfit, chisq_galfit, ndof_galfit, $
;            nfree_galfit, nfix_galfit, chi2nu_galfit]
END
