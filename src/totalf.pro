pro totalf, data, signed=signed, unsigned=unsigned
  ; Calculate the area of each pixel in cm squared
  r = 6.96e10
  area = (4.0*!pi*r^2)/(360.0*180.0)

  if arg_present(signed) then signed = total(data)*area
  if arg_present(unsigned) then unsigned = total(abs(data))*area
end
