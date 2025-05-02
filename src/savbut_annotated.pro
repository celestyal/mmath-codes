;+
; :Description:
;   Writes an annotated visualisation of the longitude averaged magnetic field
;   to the filesystem. The annotations are the ones used in my dissertation to
;   derive an estimation for the order of magnitude of meridional flow.
;   The values are hardcoded into this routine
;
; :Arguments:
;   data: in, Array
;     Longitude averaged magnetic field
;   filepath: in, String
;     Path to write the longitude average magnetic field image to
;   scale: in, Number
;     Scale factor to apply to the image before adding the axis
;   low: in, Number
;     Lowest value for the colours to become fully saturated
;   high: in, Number
;     Largest value for the colours to become fully saturated
;   yearmin: in, Number
;     Year that the x axis ticks begin at
;
;-
pro savbut_annotated, data, filename, scale, low, high, yearmin
    ; Save the current device type to return to after saving
    prevdev = !d.name

    ; Get dimensions of the image and scale
    xy = size(data, /DIMENSIONS)
    xy = xy*scale

    ; Set device to postscript
    set_plot, 'ps'

    ;+ Load greyscale colour table to improve colour visibility
    loadct, 0

    xsize=12.0
    exportdir = '../examples/kittpeak/'
    device, decomposed=0, xsize=xsize, /encapsulated
        ; Create and save image
        tvimg = bytscl(congrid(data, xy(0), xy(1)), low, high)
        tvlct, r, g, b, /GET
        imageRGB = bytarr(3, xy(0), xy(1), /NOZERO)
        imageRGB[0, *, *] = r[tvimg]
        imageRGB[1, *, *] = g[tvimg]
        imageRGB[2, *, *] = b[tvimg]
        im = image(imageRGB, /buffer, margin=[0.175, 0.075, 0.03, 0.0], image_dimensions=[xy(0), xy(1)], dimensions=xy+64)

        ;+ Draw lines and text over image to indicate coordinates

        ;+ Hardcoded coordinates with scaling applied
        a = [53, 129, 57, 147]*scale ; 17.5
        c = [112, 63, 117, 49]*scale ; -9.9
        d = [162, 40, 165, 27]*scale ; 18.3
        e = [196, 46, 200, 29]*scale ; -17.3
        f = [301, 138, 307, 156]*scale ; 12.8
        g = [328, 46, 332, 29]*scale ; -17.3
        h = [339, 133, 344, 151]*scale ; 14.6
        b = [63, 142, 65, 152]*scale ; 21.3

        ;+ Line variables
        color = 'red'
        thick = 2

        ;+ Stream (a)
        streama = polyline([[a(0), a(1)], [a(2), a(3)]], target=im, /DATA, COLOR=color, thick=thick)
        texta = TEXT(27, 129, TARGET=im, '(a)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Stream (b)
        streamb = polyline([[b(0), b(1)], [b(2), b(3)]], target=im, /DATA, COLOR=color, thick=thick)
        textb = TEXT(70, 142, TARGET=im, '(b)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Stream (c)
        streamc = polyline([[c(0), c(1)], [c(2), c(3)]], target=im, /DATA, COLOR=color, thick=thick)
        textc = TEXT(120, 49, TARGET=im, '(c)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Stream (d)
        streamd = polyline([[d(0), d(1)], [d(2), d(3)]], target=im, /DATA, COLOR=color, thick=thick)
        textd = TEXT(145, 43, TARGET=im, '(d)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Stream (e)
        streame = polyline([[e(0), e(1)], [e(2), e(3)]], target=im, /DATA, COLOR=color, thick=thick)
        texte = TEXT(200, 15, TARGET=im, '(e)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Stream (f)
        streamf = polyline([[f(0), f(1)], [f(2), f(3)]], target=im, /DATA, COLOR=color, thick=thick)
        textf = TEXT(282, 138, TARGET=im, '(f)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Stream (g)
        streamg = polyline([[g(0), g(1)], [g(2), g(3)]], target=im, /DATA, COLOR=color, thick=thick)
        textg = TEXT(332, 19, TARGET=im, '(g)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Stream (h)
        streamh = polyline([[h(0), h(1)], [h(2), h(3)]], target=im, /DATA, COLOR=color, thick=thick)
        texth = TEXT(349, 151, TARGET=im, '(h)', /DATA, COLOR=color, FONT_SIZE=9)

        ;+ Draw axis labels
        xax = axis('X', location='bottom', TICKDIR=1, MINOR=0, title='Year', coord_transform=[yearmin, 27.2753/365.25], tickfont_name='Times', tickfont_size=11)
        yax = axis('Y', location='left', TICKDIR=1, MINOR=0, coord_transform=[-1.0, (2.0/fix(xy(1), type=4))], title='Sine of Latitude', tickvalues=[-1, -sin(40*!dtor), 0, sin(40*!dtor), 1], tickfont_name='Times', tickfont_size=11)
        im.Save, filename, /transparent, border=0
    device, /close

    ; Return to previous device
    set_plot, prevdev

    ; Restore colour table
    loadct, 71
end
