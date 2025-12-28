;+
; :Description:
;   Plots the differential rotation profile of Snodgrass 1983
;
; :Arguments:
;   xsize: in, Number
;     Horizontal length of the image in centimetres
;   aspect: in, Number
;     Ratio of vertical length to horizontal length
;   lambda: in, Number
;     Latitude (absolute value) where the meridional flow falls to zero
;   u0: in, Number
;     Peak meridional flow velocity
;   fontsize: in, Number
;     Size of font in plots
;   filepath: in, String
;     Location to write the image file to
;
;-
pro mflow, xsize, aspect, lambda, u0, fontsize, filepath
    xarr = ((FINDGEN(50001)/50000)*180-90)
    u = -u0*sin(!pi*xarr/(lambda))

    ; Set mflow to 0 at latitudes above lambda
    u(where(xarr lt -abs(lambda))) = 0.0
    u(where(xarr gt abs(lambda))) = 0.0

    !P.FONT = 0
    set_plot, 'ps'
    device, /Times, xsize=xsize, ysize=xsize*aspect, filename=filepath, font_size=fontsize, /encapsulated
        plot, u, xarr, xtitle="Velocity (m/s)", ytitle="Latitude (Deg)", xstyle=0, ystyle=1, thick=2
    device, /close
end
