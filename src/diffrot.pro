;+
; :Description:
;   Plots the differential rotation profile of Snodgrass 1983
;
; :Arguments:
;   xsize: in, Number
;     Horizontal length of the image in centimetres
;   aspect: in, Number
;     Ratio of vertical length to horizontal length
;   fontsize: in, Number
;     Size of font in plots
;   filepath: in, String
;     Location to write the image file to
;
;-
pro diffrot, xsize, aspect, fontsize, filepath
    xarr = ((FINDGEN(50001)/50000)*180-90)
    Omega = 13.38 - 2.30*sin(xarr*!dtor)^2 - 1.62*sin(xarr*!dtor)^4

    !P.FONT = 0
    set_plot, 'ps'

    device, /Times, xsize=xsize, ysize=xsize*aspect, filename=filepath, font_size=fontsize, /encapsulated
        plot, Omega, xarr, xtitle="Angular Velocity (Deg/day)", ytitle="Latitude (Deg)", xstyle=0, ystyle=1, thick=2

        ;+ Plot the Carrington rotation rate too
        xvalue = 13.20
        oplot, fltarr(2) + xvalue, !y.crange, linestyle=1
    device,/close
end
