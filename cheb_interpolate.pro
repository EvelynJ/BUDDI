function cheb_interpolate, cheb, wl, wlout, ordmax=ordmax
  ; .run /home/boris/megamorph/astro-megamorph/scripts_boris/megamorph/cheb_interpolate.pro
  ; cheb_interpolate, cat[0].n_galfit_cheb, wl, z
  ; calling sequence e.g.
  ; cat = mrdfits('~/GAMA/galapagos/GAMA_benedetta_ffvqqff.fits', 1)
  ; wl=[6231, 3543, 4770, 7625, 9134, 10305, 12483, 16313, 22010]
  ; cat = mrdfits('~/GAMA/galapagos/GAMA_benedetta_ffvqqff.fits', 1, /silent)
  ; int=cheb_interpolate(cat.n_galfit_cheb, wl, (1+cat.redshift)##wl, ordmax=3)
  
  ; int=cheb_interpolate(cat.n_galfit_cheb, wl, (1+cat.redshift)##wl, ordmax=3)
  
  ; this routine does NOT work for 0-th order!! It then returns the full order again!
  
  wl = float(wl)
  wlout = float(wlout)
  if not keyword_set(ordmax) then ordmax=n_elements(cheb[0,*])-1
  
  out=fltarr(n_elements(wlout[*,0]),n_elements(cheb[0,*]))
  for i=0, n_elements(cheb[0,*])-1 do begin
    out[*,i] = chebeval(wlout[*,i], cheb[0:ordmax,i], interval=[min(wl),max(wl)])
  endfor
  
  return, out
end
