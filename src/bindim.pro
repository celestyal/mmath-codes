function bindim, bin_size, lat_lim
  sinlat = findgen(180)/89.5 - 1.0 ;sine of latitudes
  deglat = asin(sinlat)*!RADEG ;latitude in degrees
  threshold = latind(lat_lim)
  
  x = 360/bin_size
  y = fix((deglat(threshold(1))-deglat(threshold(0))) / bin_size)
  return, [x,y]
end
