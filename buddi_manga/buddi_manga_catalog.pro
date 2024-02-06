pro BUDDI_MaNGA_catalog
;set up directories for fits and results
root='/home/kjegathe/BUDDI_out/DR17_fits/'
root_files='/home/kjegathe/BUDDI/static_files/'
root_fits=root+'fits/'
root_output=root+'output/'

readcol,root_fits+'galaxies_in_progress.txt',format='a',in_prog,comment='#',/silent
readcol,root_fits+'galaxies_completed_SS.txt',format='a',completed_SS,comment='#',/silent
readcol,root_fits+'galaxies_completed_SE.txt',format='a',completed_SE,comment='#',/silent

Manga_list=mrdfits(root_files+'dapall-v3_1_1-3.1.0.fits',1)
Pymorph_list=mrdfits(root_files+'manga-pymorph-DR17.fits',2) ;extensions 1,2,3 are g,r,i band results

Manga_RA=Manga_list.OBJRA
Manga_Dec=Manga_list.OBJDEC

;avoid mismatch by rounding to 3dp
Manga_RA = Float(Round(Manga_RA*1000)/1000.0d)
Manga_Dec= Float(Round(Manga_Dec*1000)/1000.0d)

Manga_plate=Manga_list.plate
Manga_ifu=Manga_list.IFUDESIGN
Manga_plate_ifu=Manga_list.PLATEIFU

Pymorph_plate_ifu=Pymorph_list.PLATEIFU
Manga_ID_PYM=Pymorph_list.MANGA_ID
Redshift=Pymorph_list.Z
OBJID_PYM=Pymorph_list.OBJID
RA_PYM=Pymorph_list.RA
DEC_PYM=Pymorph_list.DEC


ss=0
se=0

for n=0,n_elements(in_prog)-1,1 do begin

  val_ss=where(completed_SS eq in_prog[n])
  val_se=where(completed_SE eq in_prog[n])
  if file_test(root_fits+in_prog[n]+'/IFU_decomp_1comp/decomposed_data/',/DIRECTORY) eq 1 then val_ss = 1 else val_ss = -1
  if file_test(root_fits+in_prog[n]+'/IFU_decomp_2comp/decomposed_data/',/DIRECTORY) eq 1 then val_se = 1 else val_se = -1
  print,n,' ',in_prog[n];,val_ss,val_se
  
  

  if val_ss ne -1 then begin
    res=read_sersic_results_2comp(root_fits+in_prog[n]+'/IFU_decomp_1comp/binned_images/imgblock.fits',12,bd=0)
    MAG_galfitm_SS=median(res.mag_galfit_band_d)
    MAG_galfitm_SS_g=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[2])+10^(-0.4*res.mag_galfit_band_d[3])+10^(-0.4*res.mag_galfit_band_d[4]))
    MAG_galfitm_SS_r=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[5])+10^(-0.4*res.mag_galfit_band_d[6]))
    MAG_galfitm_SS_i=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[7])+10^(-0.4*res.mag_galfit_band_d[8]))
    MAG_galfitm_SS_z=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[9])+10^(-0.4*res.mag_galfit_band_d[10]))
    RE_galfitm_cheb_A_SS=res.re_galfit_cheb_d[0]
    RE_galfitm_cheb_B_SS=res.re_galfit_cheb_d[1]
    N_galfitm_cheb_A_SS=res.n_galfit_cheb_d[0]
    N_galfitm_cheb_B_SS=res.n_galfit_cheb_d[1]
    PA_galfitm_SS=res.pa_galfit_cheb_d[0]
    Q_galfitm_SS=res.q_galfit_cheb_d[0]

    MAG_galfitm_ERR_SS=median(res.magerr_galfit_band_d)
    RE_galfitm_cheb_A_ERR_SS=res.reerr_galfit_cheb_d[0]
    RE_galfitm_cheb_B_ERR_SS=res.reerr_galfit_cheb_d[1]
    N_galfitm_cheb_A_ERR_SS=res.nerr_galfit_cheb_d[0]
    N_galfitm_cheb_B_ERR_SS=res.nerr_galfit_cheb_d[1]
    PA_galfitm_ERR_SS=res.paerr_galfit_cheb_d[0]
    Q_galfitm_ERR_SS=res.qerr_galfit_cheb_d[0]
    
  endif else begin
    MAG_galfitm_SS=999.
    MAG_galfitm_SS_g=999.
    MAG_galfitm_SS_r=999.
    MAG_galfitm_SS_i=999.
    MAG_galfitm_SS_z=999.
    RE_galfitm_cheb_A_SS=999.
    RE_galfitm_cheb_B_SS=999.
    N_galfitm_cheb_A_SS=999.
    N_galfitm_cheb_B_SS=999.
    PA_galfitm_SS=999.
    Q_galfitm_SS=999.

    MAG_galfitm_ERR_SS=999.
    RE_galfitm_cheb_A_ERR_SS=999.
    RE_galfitm_cheb_B_ERR_SS=999.
    N_galfitm_cheb_A_ERR_SS=999.
    N_galfitm_cheb_B_ERR_SS=999.
    PA_galfitm_ERR_SS=999.
    Q_galfitm_ERR_SS=999.
   
  endelse



  if val_se ne -1 then begin
    res=read_sersic_results_2comp(root_fits+in_prog[n]+'/IFU_decomp_2comp/binned_images/imgblock.fits',12,bd=1)
    MAG_galfitm_COMP1_SE=median(res.mag_galfit_band_d)
    MAG_galfitm_COMP1_SE_g=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[2])+10^(-0.4*res.mag_galfit_band_d[3])+10^(-0.4*res.mag_galfit_band_d[4]))
    MAG_galfitm_COMP1_SE_r=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[5])+10^(-0.4*res.mag_galfit_band_d[6]))
    MAG_galfitm_COMP1_SE_i=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[7])+10^(-0.4*res.mag_galfit_band_d[8]))
    MAG_galfitm_COMP1_SE_z=-2.5*alog10(10^(-0.4*res.mag_galfit_band_d[9])+10^(-0.4*res.mag_galfit_band_d[10]))
    RE_galfitm_COMP1_cheb_A_SE=res.re_galfit_cheb_d[0]
    RE_galfitm_COMP1_cheb_B_SE=res.re_galfit_cheb_d[1]
    N_galfitm_COMP1_SE=res.n_galfit_cheb_d[0]
    PA_galfitm_COMP1_SE=res.pa_galfit_cheb_d[0]
    Q_galfitm_COMP1_SE=res.q_galfit_cheb_d[0]

    MAG_galfitm_COMP2_SE=median(res.mag_galfit_band_b)
    MAG_galfitm_COMP2_SE_g=-2.5*alog10(10^(-0.4*res.mag_galfit_band_b[2])+10^(-0.4*res.mag_galfit_band_b[3])+10^(-0.4*res.mag_galfit_band_b[4]))
    MAG_galfitm_COMP2_SE_r=-2.5*alog10(10^(-0.4*res.mag_galfit_band_b[5])+10^(-0.4*res.mag_galfit_band_b[6]))
    MAG_galfitm_COMP2_SE_i=-2.5*alog10(10^(-0.4*res.mag_galfit_band_b[7])+10^(-0.4*res.mag_galfit_band_b[8]))
    MAG_galfitm_COMP2_SE_z=-2.5*alog10(10^(-0.4*res.mag_galfit_band_b[9])+10^(-0.4*res.mag_galfit_band_b[10]))
    RE_galfitm_COMP2_cheb_A_SE=res.re_galfit_cheb_b[0]
    RE_galfitm_COMP2_cheb_B_SE=res.re_galfit_cheb_b[1]
    N_galfitm_COMP2_cheb_A_SE=res.n_galfit_cheb_b[0]
    N_galfitm_COMP2_cheb_B_SE=res.n_galfit_cheb_b[1]
    PA_galfitm_COMP2_SE=res.pa_galfit_cheb_b[0]
    Q_galfitm_COMP2_SE=res.q_galfit_cheb_b[0]

    MAG_galfitm_COMP1_ERR_SE=median(res.magerr_galfit_band_d)
    RE_galfitm_COMP1_cheb_A_ERR_SE=res.reerr_galfit_cheb_d[0]
    RE_galfitm_COMP1_cheb_B_ERR_SE=res.reerr_galfit_cheb_d[1]
    N_galfitm_COMP1_ERR_SE=res.nerr_galfit_cheb_d[0]
    PA_galfitm_COMP1_ERR_SE=res.qerr_galfit_cheb_d[0]
    Q_galfitm_COMP1_ERR_SE=res.paerr_galfit_cheb_d[0]
    MAG_galfitm_COMP2_ERR_SE=median(res.magerr_galfit_band_b)
    RE_galfitm_COMP2_cheb_A_ERR_SE=res.reerr_galfit_cheb_b[0]
    RE_galfitm_COMP2_cheb_B_ERR_SE=res.reerr_galfit_cheb_b[1]
    N_galfitm_COMP2_cheb_A_ERR_SE=res.nerr_galfit_cheb_b[0]
    N_galfitm_COMP2_cheb_B_ERR_SE=res.nerr_galfit_cheb_b[1]
    PA_galfitm_COMP2_ERR_SE=res.paerr_galfit_cheb_b[0]
    Q_galfitm_COMP2_ERR_SE=res.qerr_galfit_cheb_b[0]

  endif else begin
    MAG_galfitm_COMP1_SE=999.
    MAG_galfitm_COMP1_SE_g=999.
    MAG_galfitm_COMP1_SE_r=999.
    MAG_galfitm_COMP1_SE_i=999.
    MAG_galfitm_COMP1_SE_z=999.
    RE_galfitm_COMP1_cheb_A_SE=999.
    RE_galfitm_COMP1_cheb_B_SE=999.
    N_galfitm_COMP1_SE=999.
    PA_galfitm_COMP1_SE=999.
    Q_galfitm_COMP1_SE=999.

    MAG_galfitm_COMP2_SE=999.
    MAG_galfitm_COMP2_SE_g=999.
    MAG_galfitm_COMP2_SE_r=999.
    MAG_galfitm_COMP2_SE_i=999.
    MAG_galfitm_COMP2_SE_z=999.
    RE_galfitm_COMP2_cheb_A_SE=999.
    RE_galfitm_COMP2_cheb_B_SE=999.
    N_galfitm_COMP2_cheb_A_SE=999.
    N_galfitm_COMP2_cheb_B_SE=999.
    PA_galfitm_COMP2_SE=999.
    Q_galfitm_COMP2_SE=999.

    MAG_galfitm_COMP1_ERR_SE=999.
    RE_galfitm_COMP1_cheb_A_ERR_SE=999.
    RE_galfitm_COMP1_cheb_B_ERR_SE=999.
    N_galfitm_COMP1_ERR_SE=999.
    PA_galfitm_COMP1_ERR_SE=999.
    Q_galfitm_COMP1_ERR_SE=999.

    MAG_galfitm_COMP2_ERR_SE=999.
    RE_galfitm_COMP2_cheb_A_ERR_SE=999.
    RE_galfitm_COMP2_cheb_B_ERR_SE=999.
    N_galfitm_COMP2_cheb_A_ERR_SE=999.
    N_galfitm_COMP2_cheb_B_ERR_SE=999.
    PA_galfitm_COMP2_ERR_SE=999.
    Q_galfitm_COMP2_ERR_SE=999.
    
  endelse






  if val_ss ne -1 or val_se ne -1 then begin
    element=where(strtrim(Pymorph_plate_ifu,2) eq in_prog[n])
    
    MANGA_ID=Manga_ID_PYM[element]
    PLATEIFU=in_prog[n]
    OBJID=OBJID_PYM[element]
    RA=RA_PYM[element]
    DEC=DEC_PYM[element]
    Z=Redshift[element]

    results = create_struct('MANGA_ID', MANGA_ID, $
      'PLATEIFU', PLATEIFU, $
      'OBJID', OBJID, $
      'RA', RA, $
      'DEC', DEC, $
      'Z', Z, $
      'POLY_WAVE_START', 3.621600e3, $
      'POLY_WAVE_END', 1.035380e4, $
      'MAG_galfitm_SS', MAG_galfitm_SS, $
      'MAG_galfitm_SS_g', MAG_galfitm_SS_g, $
      'MAG_galfitm_SS_r', MAG_galfitm_SS_r, $
      'MAG_galfitm_SS_i', MAG_galfitm_SS_i, $
      'MAG_galfitm_SS_z', MAG_galfitm_SS_z, $
      'RE_galfitm_cheb_A_SS', RE_galfitm_cheb_A_SS, $
      'RE_galfitm_cheb_B_SS', RE_galfitm_cheb_B_SS, $
      'N_galfitm_cheb_A_SS', N_galfitm_cheb_A_SS, $
      'N_galfitm_cheb_B_SS', N_galfitm_cheb_B_SS, $
      'PA_galfitm_SS', PA_galfitm_SS, $
      'Q_galfitm_SS', Q_galfitm_SS, $
      'MAG_galfitm_ERR_SS', MAG_galfitm_ERR_SS, $
      'RE_galfitm_cheb_A_ERR_SS', RE_galfitm_cheb_A_ERR_SS, $
      'RE_galfitm_cheb_B_ERR_SS', RE_galfitm_cheb_B_ERR_SS, $
      'N_galfitm_cheb_A_ERR_SS', N_galfitm_cheb_A_ERR_SS, $
      'N_galfitm_cheb_B_ERR_SS', N_galfitm_cheb_B_ERR_SS, $
      'PA_galfitm_ERR_SS', PA_galfitm_ERR_SS, $
      'Q_galfitm_ERR_SS', Q_galfitm_ERR_SS, $

      'MAG_galfitm_COMP1_SE', MAG_galfitm_COMP1_SE, $
      'MAG_galfitm_COMP1_SE_g', MAG_galfitm_COMP1_SE_g, $
      'MAG_galfitm_COMP1_SE_r', MAG_galfitm_COMP1_SE_r, $
      'MAG_galfitm_COMP1_SE_i', MAG_galfitm_COMP1_SE_i, $
      'MAG_galfitm_COMP1_SE_z', MAG_galfitm_COMP1_SE_z, $
      'RE_galfitm_COMP1_cheb_A_SE', RE_galfitm_COMP1_cheb_A_SE, $
      'RE_galfitm_COMP1_cheb_B_SE', RE_galfitm_COMP1_cheb_B_SE, $
      'N_galfitm_COMP1_SE', N_galfitm_COMP1_SE, $
      'PA_galfitm_COMP1_SE', PA_galfitm_COMP1_SE, $
      'Q_galfitm_COMP1_SE', Q_galfitm_COMP1_SE, $
      'MAG_galfitm_COMP2_SE', MAG_galfitm_COMP2_SE, $
      'MAG_galfitm_COMP2_SE_g', MAG_galfitm_COMP2_SE_g, $
      'MAG_galfitm_COMP2_SE_r', MAG_galfitm_COMP2_SE_r, $
      'MAG_galfitm_COMP2_SE_i', MAG_galfitm_COMP2_SE_i, $
      'MAG_galfitm_COMP2_SE_z', MAG_galfitm_COMP2_SE_z, $
      'RE_galfitm_COMP2_cheb_A_SE', RE_galfitm_COMP2_cheb_A_SE, $
      'RE_galfitm_COMP2_cheb_B_SE', RE_galfitm_COMP2_cheb_B_SE, $
      'N_galfitm_COMP2_cheb_A_SE', N_galfitm_COMP2_cheb_A_SE, $
      'N_galfitm_COMP2_cheb_B_SE', N_galfitm_COMP2_cheb_B_SE, $
      'PA_galfitm_COMP2_SE', PA_galfitm_COMP2_SE, $
      'Q_galfitm_COMP2_SE', Q_galfitm_COMP2_SE, $
      'MAG_galfitm_COMP1_ERR_SE', MAG_galfitm_COMP1_ERR_SE, $
      'RE_galfitm_COMP1_cheb_A_ERR_SE', RE_galfitm_COMP1_cheb_A_ERR_SE, $
      'RE_galfitm_COMP1_cheb_B_ERR_SE', RE_galfitm_COMP1_cheb_B_ERR_SE, $
      'N_galfitm_COMP1_ERR_SE', N_galfitm_COMP1_ERR_SE, $
      'PA_galfitm_COMP1_ERR_SE', PA_galfitm_COMP1_ERR_SE, $
      'Q_galfitm_COMP1_ERR_SE', Q_galfitm_COMP1_ERR_SE, $
      'MAG_galfitm_COMP2_ERR_SE', MAG_galfitm_COMP2_ERR_SE, $
      'RE_galfitm_COMP2_cheb_A_ERR_SE', RE_galfitm_COMP2_cheb_A_ERR_SE, $
      'RE_galfitm_COMP2_cheb_B_ERR_SE', RE_galfitm_COMP2_cheb_B_ERR_SE, $
      'N_galfitm_COMP2_cheb_A_ERR_SE', N_galfitm_COMP2_cheb_A_ERR_SE, $
      'N_galfitm_COMP2_cheb_B_ERR_SE', N_galfitm_COMP2_cheb_B_ERR_SE, $
      'PA_galfitm_COMP2_ERR_SE', PA_galfitm_COMP2_ERR_SE, $
      'Q_galfitm_COMP2_ERR_SE', Q_galfitm_COMP2_ERR_SE)

;stop
    if isa(catalog) then begin
;      input_cat=mrdfits(root_output+'BUDDI_MaNGA_cat.fits',1)
;      output_cat=[input_cat,results]
;      spawn,'mv '+root_output+'BUDDI_MaNGA_cat.fits '+root_output+'BUDDI_MaNGA_cat_old.fits'
      catalog=[catalog,results]
    endif else catalog=results


    
  endif
    
endfor
MWRFITS, catalog, root_output+'BUDDI_MaNGA_cat_final_DR17.fits'
  
end
