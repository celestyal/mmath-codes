;+
; :Description:
;   Calculates the amount of flux across a synoptic map in bins of equal
;   longitude and latitude
;
; :Returns: bin_flux
;
; :Arguments:
;   data: in, Array
;     Synoptic map data.
;   bin_size: in, Number
;     Size of the bin in each direction, given in degrees.
;
; :Keywords:
;   boundary: in, optional, Number
;     Absolute value of the upper latitude boundary to divide into bins, given
;     in degrees. Defaults to the whole map.
;
;-
function binf, data, bin_size, boundary=boundary
  if keyword_set(boundary) then begin
    threshold = latind(boundary)
  endif else begin
    boundary = 90
    threshold = [0,179]
  endelse

  sinlat = findgen(180)/89.5 - 1.0
  deglat = asin(sinlat)*!RADEG
  dim = bindim(bin_size, boundary)
  bin_flux = dblarr(dim(0), dim(1))

  long = indgen(360)

  ; Calculate total flux per each bin
  ; Iterate over bins in latitude
  for i=0, dim(1)-1 do begin
    p = where( (asin(sinlat) gt (-boundary+bin_size*(i-1))*!dtor ) and (asin(sinlat) le (-boundary+bin_size*i)*!dtor) )
    p = [min(p), max(p)]

    ; Iterate over bins in longitude
    for j=0, dim(0)-1 do begin
      q = where( (long ge (j-1)*bin_size) and (long lt j*bin_size) )
      q = [min(q), max(q)]
  
      ;Add flux to bin
      bin_flux(j, i) = total(data(q(0):q(1), p(0):p(1)))
    endfor
  endfor
  return, bin_flux
end
