function latind, limit
  ; Convert latitude argument to radians
  rad = (limit/180.0)*!pi ; Radians

  ; Construct an array of sine latitude corresponding to magnetogram observations
  sine_lat = findgen(180.0)/89.5 - 1.0

  ; Convert the sine latitudes to latitudes and take absolute value
  lat = abs(asin(sine_lat))

  ind = [-1,-1] ; (South, North)

  ; Compute lower and upper bounds for where the absolute value of latitude
  ; does not exceed rad
  for i=0, 179 do begin
    if (lat(i) lt rad) and (ind(0) eq -1) then begin
      ind(0) = i 
    endif
    if (lat(i) lt rad) and (ind(0) ne -1) then begin
      ind(1) = i
    endif
  endfor
  return, ind
end
