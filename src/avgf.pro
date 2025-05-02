;+
; :Description:
;   Calculates the physical area (cm) mapped out by the subset of the synoptic
;   map passed to the procedure
;
; :Arguments:
;   data: in, array
;     (360, 180) synoptic map
;   area: out, Number
;     Physical area (cm) of the synoptic map
;
;-
pro avgf, data, area
  r = 6.96e10 ; radius of the sun
  p_area = (4.0*!pi*r^2)/(360.0*180.0) ; area per pixel

  ; Dimensions of the array are needed to find correct area
  dim = size(data, /DIMENSIONS)

  if arg_present(area) then area=p_area*dim(0)*dim(1)
end
