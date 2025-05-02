;+
; :Description:
;   Calculates the number of bins in each direction, assuming a synoptic map
;   has dimensions (360,180)
;
; :Returns: Number of x and y bins as an array, [x,y]
;
; :Arguments:
;   bin_size: in, Number
;     Size of the bin in each direction, given in degrees.
;   boundary: in, optional, Number
;     Absolute value of the upper latitude boundary to divide into bins, given
;     in degrees. Defaults to the whole map.
;
;-
function bindim, bin_size, boundary
  sinlat = findgen(180)/89.5 - 1.0 ; sine of latitude
  deglat = asin(sinlat)*!RADEG ; latitude in degrees
  threshold = latind(boundary)
  
  ; Calculate number of longitude bins
  x = 360/bin_size

  ; Calculate number of latitude bins
  y = fix((deglat(threshold(1))-deglat(threshold(0))) / bin_size)+1
  return, [x,y]
end
