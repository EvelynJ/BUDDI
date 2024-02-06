; n1 and n2 define where the first and last GC/PSF profiles start and end

FUNCTION read_sersic_results_n, obj, nband, n1, n2;, bd=bd
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

    ;make a separate structure for each component, and then
    ;concatenate them at the end

    if n1 le 3 then feedback_3 = create_struct($
      'x_galfit_band_comp3', result.COMP3_XC, 'xerr_galfit_band_comp3',result.COMP3_XC_ERR, $
      'y_galfit_band_comp3', result.COMP3_YC, 'yerr_galfit_band_comp3',result.COMP3_YC_ERR, $
      'mag_galfit_band_comp3', result.COMP3_MAG, 'magerr_galfit_band_comp3',result.COMP3_MAG_ERR, $
      'x_galfit_cheb_comp3', res_cheb.COMP3_XC, 'xerr_galfit_cheb_comp3',res_cheb.COMP3_XC_ERR, $
      'y_galfit_cheb_comp3', res_cheb.COMP3_YC, 'yerr_galfit_cheb_comp3',res_cheb.COMP3_YC_ERR, $
      'mag_galfit_cheb_comp3', res_cheb.COMP3_MAG, 'magerr_galfit_cheb_comp3',res_cheb.COMP3_MAG_ERR)
    if n1 eq 3 then feedback=feedback_3


    if n1 le 4 and n2 ge 4 then feedback_4 = create_struct($
      'x_galfit_band_comp4', result.COMP4_XC, 'xerr_galfit_band_comp4',result.COMP4_XC_ERR, $
      'y_galfit_band_comp4', result.COMP4_YC, 'yerr_galfit_band_comp4',result.COMP4_YC_ERR, $
      'mag_galfit_band_comp4', result.COMP4_MAG, 'magerr_galfit_band_comp4',result.COMP4_MAG_ERR, $
      'x_galfit_cheb_comp4', res_cheb.COMP4_XC, 'xerr_galfit_cheb_comp4',res_cheb.COMP4_XC_ERR, $
      'y_galfit_cheb_comp4', res_cheb.COMP4_YC, 'yerr_galfit_cheb_comp4',res_cheb.COMP4_YC_ERR, $
      'mag_galfit_cheb_comp4', res_cheb.COMP4_MAG, 'magerr_galfit_cheb_comp4',res_cheb.COMP4_MAG_ERR)
    if n1 eq 4 then feedback=feedback_4

    if n1 le 5 and n2 ge 5 then feedback_5 = create_struct($
      'x_galfit_band_comp5', result.COMP5_XC, 'xerr_galfit_band_comp5',result.COMP5_XC_ERR, $
      'y_galfit_band_comp5', result.COMP5_YC, 'yerr_galfit_band_comp5',result.COMP5_YC_ERR, $
      'mag_galfit_band_comp5', result.COMP5_MAG, 'magerr_galfit_band_comp5',result.COMP5_MAG_ERR, $
      'x_galfit_cheb_comp5', res_cheb.COMP5_XC, 'xerr_galfit_cheb_comp5',res_cheb.COMP5_XC_ERR, $
      'y_galfit_cheb_comp5', res_cheb.COMP5_YC, 'yerr_galfit_cheb_comp5',res_cheb.COMP5_YC_ERR, $
      'mag_galfit_cheb_comp5', res_cheb.COMP5_MAG, 'magerr_galfit_cheb_comp5',res_cheb.COMP5_MAG_ERR)
    if n1 eq 5 then feedback=feedback_5

    if n1 le 6 and n2 ge 6 then feedback_6 = create_struct($
      'x_galfit_band_comp6', result.COMP6_XC, 'xerr_galfit_band_comp6',result.COMP6_XC_ERR, $
      'y_galfit_band_comp6', result.COMP6_YC, 'yerr_galfit_band_comp6',result.COMP6_YC_ERR, $
      'mag_galfit_band_comp6', result.COMP6_MAG, 'magerr_galfit_band_comp6',result.COMP6_MAG_ERR, $
      'x_galfit_cheb_comp6', res_cheb.COMP6_XC, 'xerr_galfit_cheb_comp6',res_cheb.COMP6_XC_ERR, $
      'y_galfit_cheb_comp6', res_cheb.COMP6_YC, 'yerr_galfit_cheb_comp6',res_cheb.COMP6_YC_ERR, $
      'mag_galfit_cheb_comp6', res_cheb.COMP6_MAG, 'magerr_galfit_cheb_comp6',res_cheb.COMP6_MAG_ERR)
    if n1 eq 6 then feedback=feedback_6

    if n1 le 7 and n2 ge 7 then feedback_7 = create_struct($
      'x_galfit_band_comp7', result.COMP7_XC, 'xerr_galfit_band_comp7',result.COMP7_XC_ERR, $
      'y_galfit_band_comp7', result.COMP7_YC, 'yerr_galfit_band_comp7',result.COMP7_YC_ERR, $
      'mag_galfit_band_comp7', result.COMP7_MAG, 'magerr_galfit_band_comp7',result.COMP7_MAG_ERR, $
      'x_galfit_cheb_comp7', res_cheb.COMP7_XC, 'xerr_galfit_cheb_comp7',res_cheb.COMP7_XC_ERR, $
      'y_galfit_cheb_comp7', res_cheb.COMP7_YC, 'yerr_galfit_cheb_comp7',res_cheb.COMP7_YC_ERR, $
      'mag_galfit_cheb_comp7', res_cheb.COMP7_MAG, 'magerr_galfit_cheb_comp7',res_cheb.COMP7_MAG_ERR)
    if n1 eq 7 then feedback=feedback_7

    if n1 le 8 and n2 ge 8 then feedback_8 = create_struct($
      'x_galfit_band_comp8', result.COMP8_XC, 'xerr_galfit_band_comp8',result.COMP8_XC_ERR, $
      'y_galfit_band_comp8', result.COMP8_YC, 'yerr_galfit_band_comp8',result.COMP8_YC_ERR, $
      'mag_galfit_band_comp8', result.COMP8_MAG, 'magerr_galfit_band_comp8',result.COMP8_MAG_ERR, $
      'x_galfit_cheb_comp8', res_cheb.COMP8_XC, 'xerr_galfit_cheb_comp8',res_cheb.COMP8_XC_ERR, $
      'y_galfit_cheb_comp8', res_cheb.COMP8_YC, 'yerr_galfit_cheb_comp8',res_cheb.COMP8_YC_ERR, $
      'mag_galfit_cheb_comp8', res_cheb.COMP8_MAG, 'magerr_galfit_cheb_comp8',res_cheb.COMP8_MAG_ERR)
    if n1 eq 8 then feedback=feedback_8

    if n1 le 9 and n2 ge 9 then feedback_9 = create_struct($
      'x_galfit_band_comp9', result.COMP9_XC, 'xerr_galfit_band_comp9',result.COMP9_XC_ERR, $
      'y_galfit_band_comp9', result.COMP9_YC, 'yerr_galfit_band_comp9',result.COMP9_YC_ERR, $
      'mag_galfit_band_comp9', result.COMP9_MAG, 'magerr_galfit_band_comp9',result.COMP9_MAG_ERR, $
      'x_galfit_cheb_comp9', res_cheb.COMP9_XC, 'xerr_galfit_cheb_comp9',res_cheb.COMP9_XC_ERR, $
      'y_galfit_cheb_comp9', res_cheb.COMP9_YC, 'yerr_galfit_cheb_comp9',res_cheb.COMP9_YC_ERR, $
      'mag_galfit_cheb_comp9', res_cheb.COMP9_MAG, 'magerr_galfit_cheb_comp9',res_cheb.COMP9_MAG_ERR)
    if n1 eq 9 then feedback=feedback_9

    if n1 le 10 and n2 ge 10 then feedback_10 = create_struct($
      'x_galfit_band_comp10', result.COMP10_XC, 'xerr_galfit_band_comp10',result.COMP10_XC_ERR, $
      'y_galfit_band_comp10', result.COMP10_YC, 'yerr_galfit_band_comp10',result.COMP10_YC_ERR, $
      'mag_galfit_band_comp10', result.COMP10_MAG, 'magerr_galfit_band_comp10',result.COMP10_MAG_ERR, $
      'x_galfit_cheb_comp10', res_cheb.COMP10_XC, 'xerr_galfit_cheb_comp10',res_cheb.COMP10_XC_ERR, $
      'y_galfit_cheb_comp10', res_cheb.COMP10_YC, 'yerr_galfit_cheb_comp10',res_cheb.COMP10_YC_ERR, $
      'mag_galfit_cheb_comp10', res_cheb.COMP10_MAG, 'magerr_galfit_cheb_comp10',res_cheb.COMP10_MAG_ERR)
    if n1 eq 9 then feedback=feedback_10

    if n1 le 11 and n2 ge 11 then feedback_11 = create_struct($
      'x_galfit_band_comp11', result.COMP11_XC, 'xerr_galfit_band_comp11',result.COMP11_XC_ERR, $
      'y_galfit_band_comp11', result.COMP11_YC, 'yerr_galfit_band_comp11',result.COMP11_YC_ERR, $
      'mag_galfit_band_comp11', result.COMP11_MAG, 'magerr_galfit_band_comp11',result.COMP11_MAG_ERR, $
      'x_galfit_cheb_comp11', res_cheb.COMP11_XC, 'xerr_galfit_cheb_comp11',res_cheb.COMP11_XC_ERR, $
      'y_galfit_cheb_comp11', res_cheb.COMP11_YC, 'yerr_galfit_cheb_comp11',res_cheb.COMP11_YC_ERR, $
      'mag_galfit_cheb_comp11', res_cheb.COMP11_MAG, 'magerr_galfit_cheb_comp11',res_cheb.COMP11_MAG_ERR)
    if n1 eq 11 then feedback=feedback_11

    if n1 le 12 and n2 ge 12 then feedback_12 = create_struct($
      'x_galfit_band_comp12', result.COMP12_XC, 'xerr_galfit_band_comp12',result.COMP12_XC_ERR, $
      'y_galfit_band_comp12', result.COMP12_YC, 'yerr_galfit_band_comp12',result.COMP12_YC_ERR, $
      'mag_galfit_band_comp12', result.COMP12_MAG, 'magerr_galfit_band_comp12',result.COMP12_MAG_ERR, $
      'x_galfit_cheb_comp12', res_cheb.COMP12_XC, 'xerr_galfit_cheb_comp12',res_cheb.COMP12_XC_ERR, $
      'y_galfit_cheb_comp12', res_cheb.COMP12_YC, 'yerr_galfit_cheb_comp12',res_cheb.COMP12_YC_ERR, $
      'mag_galfit_cheb_comp12', res_cheb.COMP12_MAG, 'magerr_galfit_cheb_comp12',res_cheb.COMP12_MAG_ERR)
    if n1 eq 12 then feedback=feedback_12

    if n1 le 13 and n2 ge 13 then feedback_13 = create_struct($
      'x_galfit_band_comp13', result.COMP13_XC, 'xerr_galfit_band_comp13',result.COMP13_XC_ERR, $
      'y_galfit_band_comp13', result.COMP13_YC, 'yerr_galfit_band_comp13',result.COMP13_YC_ERR, $
      'mag_galfit_band_comp13', result.COMP13_MAG, 'magerr_galfit_band_comp13',result.COMP13_MAG_ERR, $
      'x_galfit_cheb_comp13', res_cheb.COMP13_XC, 'xerr_galfit_cheb_comp13',res_cheb.COMP13_XC_ERR, $
      'y_galfit_cheb_comp13', res_cheb.COMP13_YC, 'yerr_galfit_cheb_comp13',res_cheb.COMP13_YC_ERR, $
      'mag_galfit_cheb_comp13', res_cheb.COMP13_MAG, 'magerr_galfit_cheb_comp13',res_cheb.COMP13_MAG_ERR)
    if n1 eq 13 then feedback=feedback_13

    if n1 le 14 and n2 ge 14 then feedback_14 = create_struct($
      'x_galfit_band_comp14', result.COMP14_XC, 'xerr_galfit_band_comp14',result.COMP14_XC_ERR, $
      'y_galfit_band_comp14', result.COMP14_YC, 'yerr_galfit_band_comp14',result.COMP14_YC_ERR, $
      'mag_galfit_band_comp14', result.COMP14_MAG, 'magerr_galfit_band_comp14',result.COMP14_MAG_ERR, $
      'x_galfit_cheb_comp14', res_cheb.COMP14_XC, 'xerr_galfit_cheb_comp14',res_cheb.COMP14_XC_ERR, $
      'y_galfit_cheb_comp14', res_cheb.COMP14_YC, 'yerr_galfit_cheb_comp14',res_cheb.COMP14_YC_ERR, $
      'mag_galfit_cheb_comp14', res_cheb.COMP14_MAG, 'magerr_galfit_cheb_comp14',res_cheb.COMP14_MAG_ERR)
    if n1 eq 14 then feedback=feedback_14

    if n1 le 15 and n2 ge 15 then feedback_15 = create_struct($
      'x_galfit_band_comp15', result.COMP15_XC, 'xerr_galfit_band_comp15',result.COMP15_XC_ERR, $
      'y_galfit_band_comp15', result.COMP15_YC, 'yerr_galfit_band_comp15',result.COMP15_YC_ERR, $
      'mag_galfit_band_comp15', result.COMP15_MAG, 'magerr_galfit_band_comp15',result.COMP15_MAG_ERR, $
      'x_galfit_cheb_comp15', res_cheb.COMP15_XC, 'xerr_galfit_cheb_comp15',res_cheb.COMP15_XC_ERR, $
      'y_galfit_cheb_comp15', res_cheb.COMP15_YC, 'yerr_galfit_cheb_comp15',res_cheb.COMP15_YC_ERR, $
      'mag_galfit_cheb_comp15', res_cheb.COMP15_MAG, 'magerr_galfit_cheb_comp15',res_cheb.COMP15_MAG_ERR)
    if n1 eq 15 then feedback=feedback_15

    if n1 le 16 and n2 ge 16 then feedback_16 = create_struct($
      'x_galfit_band_comp16', result.COMP16_XC, 'xerr_galfit_band_comp16',result.COMP16_XC_ERR, $
      'y_galfit_band_comp16', result.COMP16_YC, 'yerr_galfit_band_comp16',result.COMP16_YC_ERR, $
      'mag_galfit_band_comp16', result.COMP16_MAG, 'magerr_galfit_band_comp16',result.COMP16_MAG_ERR, $
      'x_galfit_cheb_comp16', res_cheb.COMP16_XC, 'xerr_galfit_cheb_comp16',res_cheb.COMP16_XC_ERR, $
      'y_galfit_cheb_comp16', res_cheb.COMP16_YC, 'yerr_galfit_cheb_comp16',res_cheb.COMP16_YC_ERR, $
      'mag_galfit_cheb_comp16', res_cheb.COMP16_MAG, 'magerr_galfit_cheb_comp16',res_cheb.COMP16_MAG_ERR)
    if n1 eq 16 then feedback=feedback_16

    if n1 le 17 and n2 ge 17 then feedback_17 = create_struct($
      'x_galfit_band_comp17', result.COMP17_XC, 'xerr_galfit_band_comp17',result.COMP17_XC_ERR, $
      'y_galfit_band_comp17', result.COMP17_YC, 'yerr_galfit_band_comp17',result.COMP17_YC_ERR, $
      'mag_galfit_band_comp17', result.COMP17_MAG, 'magerr_galfit_band_comp17',result.COMP17_MAG_ERR, $
      'x_galfit_cheb_comp17', res_cheb.COMP17_XC, 'xerr_galfit_cheb_comp17',res_cheb.COMP17_XC_ERR, $
      'y_galfit_cheb_comp17', res_cheb.COMP17_YC, 'yerr_galfit_cheb_comp17',res_cheb.COMP17_YC_ERR, $
      'mag_galfit_cheb_comp17', res_cheb.COMP17_MAG, 'magerr_galfit_cheb_comp17',res_cheb.COMP17_MAG_ERR)
    if n1 eq 17 then feedback=feedback_17

    if n1 le 18 and n2 ge 18 then feedback_18 = create_struct($
      'x_galfit_band_comp18', result.COMP18_XC, 'xerr_galfit_band_comp18',result.COMP18_XC_ERR, $
      'y_galfit_band_comp18', result.COMP18_YC, 'yerr_galfit_band_comp18',result.COMP18_YC_ERR, $
      'mag_galfit_band_comp18', result.COMP18_MAG, 'magerr_galfit_band_comp18',result.COMP18_MAG_ERR, $
      'x_galfit_cheb_comp18', res_cheb.COMP18_XC, 'xerr_galfit_cheb_comp18',res_cheb.COMP18_XC_ERR, $
      'y_galfit_cheb_comp18', res_cheb.COMP18_YC, 'yerr_galfit_cheb_comp18',res_cheb.COMP18_YC_ERR, $
      'mag_galfit_cheb_comp18', res_cheb.COMP18_MAG, 'magerr_galfit_cheb_comp18',res_cheb.COMP18_MAG_ERR)
    if n1 eq 18 then feedback=feedback_18

    if n1 le 19 and n2 ge 19 then feedback_19 = create_struct($
      'x_galfit_band_comp19', result.COMP19_XC, 'xerr_galfit_band_comp19',result.COMP19_XC_ERR, $
      'y_galfit_band_comp19', result.COMP19_YC, 'yerr_galfit_band_comp19',result.COMP19_YC_ERR, $
      'mag_galfit_band_comp19', result.COMP19_MAG, 'magerr_galfit_band_comp19',result.COMP19_MAG_ERR, $
      'x_galfit_cheb_comp19', res_cheb.COMP19_XC, 'xerr_galfit_cheb_comp19',res_cheb.COMP19_XC_ERR, $
      'y_galfit_cheb_comp19', res_cheb.COMP19_YC, 'yerr_galfit_cheb_comp19',res_cheb.COMP19_YC_ERR, $
      'mag_galfit_cheb_comp19', res_cheb.COMP19_MAG, 'magerr_galfit_cheb_comp19',res_cheb.COMP19_MAG_ERR)
    if n1 eq 19 then feedback=feedback_19

    if n1 le 20 and n2 ge 20 then feedback_20 = create_struct($
      'x_galfit_band_comp20', result.COMP20_XC, 'xerr_galfit_band_comp20',result.COMP20_XC_ERR, $
      'y_galfit_band_comp20', result.COMP20_YC, 'yerr_galfit_band_comp20',result.COMP20_YC_ERR, $
      'mag_galfit_band_comp20', result.COMP20_MAG, 'magerr_galfit_band_comp20',result.COMP20_MAG_ERR, $
      'x_galfit_cheb_comp20', res_cheb.COMP20_XC, 'xerr_galfit_cheb_comp20',res_cheb.COMP20_XC_ERR, $
      'y_galfit_cheb_comp20', res_cheb.COMP20_YC, 'yerr_galfit_cheb_comp20',res_cheb.COMP20_YC_ERR, $
      'mag_galfit_cheb_comp20', res_cheb.COMP20_MAG, 'magerr_galfit_cheb_comp20',res_cheb.COMP20_MAG_ERR)
    if n1 eq 20 then feedback=feedback_20

    if n1 le 21 and n2 ge 21 then feedback_21 = create_struct($
      'x_galfit_band_comp21', result.COMP21_XC, 'xerr_galfit_band_comp21',result.COMP21_XC_ERR, $
      'y_galfit_band_comp21', result.COMP21_YC, 'yerr_galfit_band_comp21',result.COMP21_YC_ERR, $
      'mag_galfit_band_comp21', result.COMP21_MAG, 'magerr_galfit_band_comp21',result.COMP21_MAG_ERR, $
      'x_galfit_cheb_comp21', res_cheb.COMP21_XC, 'xerr_galfit_cheb_comp21',res_cheb.COMP21_XC_ERR, $
      'y_galfit_cheb_comp21', res_cheb.COMP21_YC, 'yerr_galfit_cheb_comp21',res_cheb.COMP21_YC_ERR, $
      'mag_galfit_cheb_comp21', res_cheb.COMP21_MAG, 'magerr_galfit_cheb_comp21',res_cheb.COMP21_MAG_ERR)
    if n1 eq 21 then feedback=feedback_21

    if n1 le 22 and n2 ge 22 then feedback_22 = create_struct($
      'x_galfit_band_comp22', result.COMP22_XC, 'xerr_galfit_band_comp22',result.COMP22_XC_ERR, $
      'y_galfit_band_comp22', result.COMP22_YC, 'yerr_galfit_band_comp22',result.COMP22_YC_ERR, $
      'mag_galfit_band_comp22', result.COMP22_MAG, 'magerr_galfit_band_comp22',result.COMP22_MAG_ERR, $
      'x_galfit_cheb_comp22', res_cheb.COMP22_XC, 'xerr_galfit_cheb_comp22',res_cheb.COMP22_XC_ERR, $
      'y_galfit_cheb_comp22', res_cheb.COMP22_YC, 'yerr_galfit_cheb_comp22',res_cheb.COMP22_YC_ERR, $
      'mag_galfit_cheb_comp22', res_cheb.COMP22_MAG, 'magerr_galfit_cheb_comp22',res_cheb.COMP22_MAG_ERR)
    if n1 eq 22 then feedback=feedback_22

    if n1 le 23 and n2 ge 23 then feedback_23 = create_struct($
      'x_galfit_band_comp23', result.COMP23_XC, 'xerr_galfit_band_comp23',result.COMP23_XC_ERR, $
      'y_galfit_band_comp23', result.COMP23_YC, 'yerr_galfit_band_comp23',result.COMP23_YC_ERR, $
      'mag_galfit_band_comp23', result.COMP23_MAG, 'magerr_galfit_band_comp23',result.COMP23_MAG_ERR, $
      'x_galfit_cheb_comp23', res_cheb.COMP23_XC, 'xerr_galfit_cheb_comp23',res_cheb.COMP23_XC_ERR, $
      'y_galfit_cheb_comp23', res_cheb.COMP23_YC, 'yerr_galfit_cheb_comp23',res_cheb.COMP23_YC_ERR, $
      'mag_galfit_cheb_comp23', res_cheb.COMP23_MAG, 'magerr_galfit_cheb_comp23',res_cheb.COMP23_MAG_ERR)
    if n1 eq 23 then feedback=feedback_23

    if n1 le 24 and n2 ge 24 then feedback_24 = create_struct($
      'x_galfit_band_comp24', result.COMP24_XC, 'xerr_galfit_band_comp24',result.COMP24_XC_ERR, $
      'y_galfit_band_comp24', result.COMP24_YC, 'yerr_galfit_band_comp24',result.COMP24_YC_ERR, $
      'mag_galfit_band_comp24', result.COMP24_MAG, 'magerr_galfit_band_comp24',result.COMP24_MAG_ERR, $
      'x_galfit_cheb_comp24', res_cheb.COMP24_XC, 'xerr_galfit_cheb_comp24',res_cheb.COMP24_XC_ERR, $
      'y_galfit_cheb_comp24', res_cheb.COMP24_YC, 'yerr_galfit_cheb_comp24',res_cheb.COMP24_YC_ERR, $
      'mag_galfit_cheb_comp24', res_cheb.COMP24_MAG, 'magerr_galfit_cheb_comp24',res_cheb.COMP24_MAG_ERR)
    if n1 eq 24 then feedback=feedback_24


    ;feedback=create_struct('temp','blank')
    if n1 lt 3 then feedback=[feedback,feedback_3]
    if n1 lt 4 and n2 ge 4 then feedback=[feedback,feedback_4]
    if n1 lt 5 and n2 ge 5 then feedback=[feedback,feedback_5]
    if n1 lt 6 and n2 ge 6 then feedback=[feedback,feedback_6]
    if n1 lt 7 and n2 ge 7 then feedback=[feedback,feedback_7]
    if n1 lt 8 and n2 ge 8 then feedback=[feedback,feedback_8]
    if n1 lt 9 and n2 ge 9 then feedback=[feedback,feedback_9]
    if n1 lt 10 and n2 ge 10 then feedback=[feedback,feedback_10]
    if n1 lt 11 and n2 ge 11 then feedback=[feedback,feedback_11]
    if n1 lt 12 and n2 ge 12 then feedback=[feedback,feedback_12]
    if n1 lt 13 and n2 ge 13 then feedback=[feedback,feedback_13]
    if n1 lt 14 and n2 ge 14 then feedback=[feedback,feedback_14]
    if n1 lt 15 and n2 ge 15 then feedback=[feedback,feedback_15]
    if n1 lt 16 and n2 ge 16 then feedback=[feedback,feedback_16]
    if n1 lt 17 and n2 ge 17 then feedback=[feedback,feedback_17]
    if n1 lt 18 and n2 ge 18 then feedback=[feedback,feedback_18]
    if n1 lt 19 and n2 ge 19 then feedback=[feedback,feedback_19]
    if n1 lt 20 and n2 ge 20 then feedback=[feedback,feedback_20]
    if n1 lt 21 and n2 ge 21 then feedback=[feedback,feedback_21]
    if n1 lt 22 and n2 ge 22 then feedback=[feedback,feedback_22]
    if n1 lt 23 and n2 ge 23 then feedback=[feedback,feedback_23]
    if n1 lt 24 and n2 ge 24 then feedback=[feedback,feedback_24]

    ; TO BE ADDED:
    ; fitting time
    ; NEIGH_GALFIT HAS TO BE ADAPTED! WHY??


    ; to include:
    ; there is more band_info which is not used yet (band, wl, datain,
    ; sigma, MASL, magzpt) Not sure we'll need them!

  ENDIF ELSE BEGIN
    psf=strarr(nband)
    for n=0,nband-1 do psf[n]='none'


    ; TO BE ADDED:
    ; fitting time
    ; NEIGH_GALFIT HAS TO BE ADAPTED!


  ENDELSE
  return, feedback

  ;[mag, magerr, re, reerr, n, nerr, q, qerr, pa, paerr, $
  ;            x, xerr, y, yerr, sky, neigh_galfit, chisq_galfit, ndof_galfit, $
  ;            nfree_galfit, nfix_galfit, chi2nu_galfit]
END
