pro avgf, data, signed=signed, unsigned=unsigned
  r = 6.96e10 ; radius of the sun
  area = (4.0*!pi*r^2)/(360.0*180.0) ; area per pixel

  ; Dimensions of the array are needed to find correct area
  dim = size(data, /DIMENSIONS)

  if (arg_present(signed)) and (arg_present(unsigned)) then begin
    totalf, data, signed=signed, unsigned=unsigned
  endif else begin
    if arg_present(unsigned) then totalf, data, unsigned=unsigned
    if arg_present(signed) then totalf, data, signed=signed
  endelse

  if arg_present(signed) then signed = signed/(area*dim(0)*dim(1))
  if arg_present(unsigned) then unsigned = unsigned/(area*dim(0)*dim(1))
end
