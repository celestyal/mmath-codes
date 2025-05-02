;+
; :Description:
;   Displays the longitude averaged field and calculates the meridional flow
;   velocity depending on where the cursor is clicked.
;
; :Arguments:
;   data: in, Array
;     Longitude averaged magnetic field
;   scale: in, Number
;     Scale factor to apply to the image before adding the axis
;   low: in, Number
;     Lowest value for the colours to become fully saturated
;   high: in, Number
;     Largest value for the colours to become fully saturated
;   yearmin: in, Number
;     Year that the x axis ticks begin at
;   repeats: in, Number
;     Number of times to repeat the calculation process
;
;-
pro dispbut, data, scale, low, high, repeats
    ; Save the current device type to return to after saving
    prevdev = !d.name

    ; Get dimensions of the image and scale
    xy = size(data, /DIMENSIONS)
    xy = xy*scale

    ; Set plot to X to display on screen
    set_plot, 'x'
    device, decomposed=0

    ; Create appropriately sized window
    window, xsize=xy(0), ysize=xy(1)

    ; Create and save image
    tvimg = bytscl(congrid(data, xy(0), xy(1)), low, high)
    tvlct, r, g, b, /GET
    imageRGB = bytarr(3, xy(0), xy(1), /NOZERO)
    imageRGB[0, *, *] = r[tvimg]
    imageRGB[1, *, *] = g[tvimg]
    imageRGB[2, *, *] = b[tvimg]
    tv, tvimg

    for i=0, repeats do begin
        ; Track where cursor presses and to mark stream edges
        CURSOR, t0, s0, /DEVICE, /down
        CURSOR, t1, s1, /DEVICE, /down

        ;
        ; Scale
        ;
        sc = float(scale)
        t0 = floor(float(t0)/sc)
        t1 = floor(float(t1)/sc)
        s0 = floor(float(s0)/sc)
        s1 = floor(float(s1)/sc)

        ; Sine of latitude
        s0 = (s0-90)*1.0/90.0
        s1 = (s1-90)*1.0/90.0

        ; Convert to latitude
        s0 = asin(s0)
        s1 = asin(s1)

        ; Angular velocity
        omega = (s1-s0)/(t1-t0)

        ; Velocity
        rsun = 6.96e8
        velocity = rsun*omega

        ; Calculate velocity
        print, velocity/(27.2753*86400), " metres per second"
    endfor

    ; Return to previous device
    set_plot, prevdev
end
