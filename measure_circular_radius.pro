; 
; this code will read in the centres of a series of binned 
; spectra along with the x and y position of the centre of 
; the galaxy and the position angle. Using this information, 
; it assumes that the true radius of the bin can be calculated 
; by solving the equation for the ellipse
; 
; 

pro measure_circular_radius,bin,x_bin,y_bin,x0,y0,PA, radii

radii=fltarr(n_elements(bin))
r_temp=fltarr(n_elements(bin))

; start by identifying orientation: NE-SW or NE-SW
if PA lt 0 then PA+=180
if PA gt 0   and PA le 90  then orientation='NESW'
if PA gt 180 and PA le 270 then orientation='NESW'
if PA gt 90  and PA le 180 then orientation='NWSE'
if PA gt 270 and PA le 360 then orientation='NWSE'

if PA gt 180 then PA_prime=PA-180 else PA_prime=PA
for n=0,n_elements(bin)-1,1 do begin
  x=x_bin[n]-x0
  y=y_bin[n]-y0
  
  ;find r, distance of bin from the centre as seen on the image
  ; and identify which side of the semi-major axis it lies
  r=sqrt(x*x+y*y)
  r_temp[n]=r
  
  if orientation eq 'NESW' then begin
    
    PA_r=PA/360*2*!pi                 ;convert PA to radians
    PA_prime_r=PA_prime/360*2*!pi
    y_minor=x*tan(PA_prime_r)               ;y value of minor axis at the x-posiiton of the bin
;    if y lt y_minor then radii[n]=-r else radii[n]=r
    
    ;find theta, the angle between the major axis and the bin
    y_major=x*tan(PA_prime_r)               ;y value of minor axis at the x-posiiton of the bin
    x_major=y/tan(PA_prime_r)               ;y value of minor axis at the x-posiiton of the bin
    if (x lt 0 and y ge y_major) or (x ge 0 and y lt y_major) then begin
      phi=acos(x/r)
      theta=phi-((!pi/2)-PA_prime_r)
    endif else if (y ge 0 and y lt y_major) or (y lt 0 and y ge y_major) then begin
      phi=acos(x/r)
      theta=((!pi/2)-PA_prime_r)-phi
    endif else if (x lt 0 and y ge y_minor) or (x ge 0 and y lt y_minor) then begin
      theta_prime=atan(y/x)       ;angle between bin and x-axis
      theta=((!pi/2)-PA_prime_r)+theta_prime
    endif else if (x ge 0 and y ge y_minor) or (x lt 0 and y lt y_minor) then begin
      theta_prime=atan(x/y)
      theta=PA_prime_r+theta_prime
    endif   
    
    ;find major axis
    x_prime=r*cos(theta)    ;x-distance of bin along major axis
    a=(x_prime)/cos(theta)
 
    if y lt y_minor then radii[n]=a else radii[n]=-a
       
  endif else begin
    PA_r=PA/360*2*!pi                 ;convert PA to radians
    PA_prime_r=PA_prime/360*2*!pi
    y_minor=-x/tan(PA_prime_r-(!pi/2))               ;y value of minor axis at the x-posiiton of the bin
;    if y lt y_minor then radii[n]=-r else radii[n]=r
    
    ;find theta, the angle between the major axis and the bin
    y_major=x*tan(PA_prime_r-(!pi/2))               ;y value of minor axis at the x-posiiton of the bin
    x_major=y/tan(PA_prime_r-(!pi/2))               ;y value of minor axis at the x-posiiton of the bin
    if (x ge 0 and y ge y_major) or (x lt 0 and y lt y_major) then begin
      phi=acos(x/r)
      theta=phi-(PA_prime_r-(!pi/2))
    endif else if (y ge 0 and y lt y_major) or (y lt 0 and y ge y_major) then begin
      phi=acos(x/r)
      theta=(PA_prime_r-(!pi/2))-phi
    endif else if (x lt 0 and y ge y_minor) or (x ge 0 and y lt y_minor) then begin
      theta_prime=atan(y/x)       ;angle between bin and x-axis
      theta=(PA_prime_r-(!pi/2))+theta_prime
    endif else if (x ge 0 and y ge y_minor) or (x lt 0 and y lt y_minor) then begin
      theta_prime=atan(x/y)
      theta=(!pi-PA_prime_r)+theta_prime
    endif  
    
    ;find major axis
    x_prime=r*cos(theta)    ;x-distance of bin along major axis
    a=(x_prime)/cos(theta)  ;length of major axis
 
    if y lt y_minor then radii[n]=-a else radii[n]=a
       
    
  endelse
  ;print,x,y, r,radii[n]
endfor





end