function linear_interpolate,x_out,x_in,y_in

x_out_total=n_elements(x_out)-1
x_in_total=n_elements(x_in)-1
y_out=fltarr(x_out_total+1)

for n=0,x_out_total,1 do begin
  temp_val=where(x_in gt x_out[n],count)
  if count eq n_elements(x_in) then begin
    x0=x_in[0]
    x1=x_in[1]
    y0=y_in[0]
    y1=y_in[1]
  endif else if count eq 0 then begin
    x0=x_in[x_in_total-1]
    x1=x_in[x_in_total]
    y0=y_in[x_in_total-1]
    y1=y_in[x_in_total]
  endif else begin
    x0=x_in[temp_val[0]-1]
    x1=x_in[temp_val[0]]
    y0=y_in[temp_val[0]-1]
    y1=y_in[temp_val[0]]
  endelse
  
  y_out[n]=y0+((y1-y0)*(x_out[n]-x0)/(x1-x0)) 

endfor

return,y_out

end


