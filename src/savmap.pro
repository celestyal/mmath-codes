;+
; :Description:
;   Writes a visualisation of the a synoptic map to the filesystem
;
; :Arguments:
;   data: in, Array
;     Synoptic map array
;   filepath: in, String
;     Path to write the longitude average magnetic field image to
;   scale: in, Number
;     Scale factor to apply to the image before adding the axis
;   low: in, Number
;     Lowest value for the colours to become fully saturated
;   high: in, Number
;     Largest value for the colours to become fully saturated
;
;-
pro savmap, data, filepath, scale, low, high
    ; Save the current device type to return to after saving
    prevdev = !d.name

    ; Get dimensions of the image and scale
    xy = size(data, /DIMENSIONS)
    xy = xy*scale

    ; Set device to postscript
    set_plot, 'ps'

    xsize=12.0
    exportdir = '../examples/kittpeak/'
    device, decomposed=0, xsize=xsize, /encapsulated
        ; Create and save image
        tvimg = bytscl(congrid(data, xy[0], xy[1]), low, high)
        tvlct, r, g, b, /GET
        imageRGB = bytarr(3, xy[0], xy[1], /NOZERO)
        imageRGB[0, *, *] = r[tvimg]
        imageRGB[1, *, *] = g[tvimg]
        imageRGB[2, *, *] = b[tvimg]
        im = image(imageRGB, /buffer, margin=[0.175, 0.075, 0.03, 0.0], image_dimensions=[xy[0], xy[1]], dimensions=xy+64)
        xax = axis('X', location='bottom', TICKDIR=1, MINOR=0, title='Carrington Longitude (Deg)', coord_transform=[0.0, (360.0/fix(xy[0], type=4))], tickvalues=[0, 60, 120, 180, 240, 300, 360], tickfont_name='Times', tickfont_size=10)
        yax = axis('Y', location='left', TICKDIR=1, MINOR=0, coord_transform=[-1.0, (2.0/fix(xy[1], type=4))], title='Sine of Latitude', tickvalues=[-1, -sin(40*!dtor), 0, sin(40*!dtor), 1], tickfont_name='Times', tickfont_size=10)
        im.Save, filepath, /transparent, border=0
    device, /close

    ; Return to previous device
    set_plot, prevdev
end
