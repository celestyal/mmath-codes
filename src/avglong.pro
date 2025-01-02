pro avglong, data, signed=signed, unsigned=unsigned
  if arg_present(signed) then begin
    signed = dblarr(180)
    for i=0, 179 do begin
      signed(i) = total(data(*,i))
    endfor
    signed = signed/360
  endif
  if arg_present(unsigned) then begin
    unsigned = dblarr(180)
    for i=0, 179 do begin
      unsigned(i)= total(abs(data(*,i)))
    endfor
    unsigned = unsigned/360
  endif
end
