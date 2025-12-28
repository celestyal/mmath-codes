;+
; :Description:
;   Processes NSO/Kitt Peak synoptic maps for CR1625-CR2007 (360x180) and
;   generates every plot (plus additional ones) used in the section
;   "Magnetic Field Observations"
;-
pro mmathObservations
    ;+
    ; I/O Parameters
    ;-

    ;+ x-dimension of image
    xsize = 8.0

    ;+ Font size
    fontsize = 10.0

    ;+ Aspect ratios
    aspect1 = 3.0/4.0
    aspect2 = 9.0/16.0
    aspect3 = 1.0

    ;+ Colour table to use
    loadct, 71

    ;+ Default ystyle to use on plots
    ystyle = 2

    ;+ Tick size (and direction)
    tick_orig = !P.TICKLEN
    !P.TICKLEN = -abs(!P.TICKLEN)

    ;+ Directory containing FITS (to read)
    fit_dir = "./data/mag/"

    ;+ Filepath of the Carrington Rotation metadata (to read)
    carr_file = "./data/dates/carrington.rotations"

    ;+ Root export directory
    exportdir = "./examples/magnetic-observations/"

    ;+ Directory to read/write script variables to
    variables_path = exportdir + "processed-field.sav"

    ;+ Directory to save total flux plots to
    flux_dir = exportdir + "flux/"

    ;+ Directory to save map visualisations to
    rotation_dir = exportdir + "maps/"

    ;+ Directory to save butterfly diagram to
    butterfly_dir = exportdir + "butterfly/"

    ; Directories to save active longitude plots
    longitude_dir = exportdir + "active-longitudes/"
    cycle22_dir = longitude_dir + "cycle22/"
    yearly_dir = longitude_dir + "yearly/"
    biyearly_dir = longitude_dir + "biyearly/"

    ;+
    ; End of parameters
    ;-

    ;+ Create output directories if they doesn't exist
    file_mkdir, exportdir
    file_mkdir, rotation_dir
    file_mkdir, longitude_dir
    file_mkdir, cycle22_dir
    file_mkdir, yearly_dir
    file_mkdir, biyearly_dir
    file_mkdir, butterfly_dir
    file_mkdir, flux_dir

    ;+ Extract rotation number from filenames (assumes they're correct)
    filenames = file_basename(file_search(fit_dir + '*.fits'))

    ;+ Recover Carrington rotation number from fit filenames
    rot = stregex(filenames, '[0-9]+', /extract)
    n_rot = n_elements(rot) ;number of carrington rotations

    ;+ Read other Carrington rotation metadata
    carrdata, rot, carr_file, date=date, yy=yy


    ;+ Prepare array to store each synoptic map
    data = dblarr(360, 180, n_rot)

    ; Prepare total flux array
    flux_s = dblarr(n_rot)
    flux_u = dblarr(n_rot)

    ; Prepare positive and negative polarity total flux arrays
    flux_p = dblarr(n_rot)
    flux_n = dblarr(n_rot)

    ; Prepare polar total flux arrays
    flux_pole_n = dblarr(n_rot)
    flux_pole_s = dblarr(n_rot)
    uflux_pole_n = dblarr(n_rot)
    uflux_pole_s = dblarr(n_rot)

    ; Prepare mid latitude total flux arrays
    flux_mid_n = dblarr(n_rot)
    flux_mid_s = dblarr(n_rot)
    uflux_mid_n = dblarr(n_rot)
    uflux_mid_s = dblarr(n_rot)

    ; Prepare active longitude total flux arrays
    flux_active_n = dblarr(n_rot)
    flux_active_s = dblarr(n_rot)
    uflux_active_n = dblarr(n_rot)
    uflux_active_s = dblarr(n_rot)

    ;+ Prepare longitude average field array
    flux_s_long = dblarr(n_rot, 180)

    ;+ Check if the fits have been saved to file already to avoid re-reading them
    filetest = file_test(variables_path)

    ;+ Check if fits have been processed already
    if (filetest eq 1) then begin
        ;+ Restore saved fit data
        restore, variables_path
    endif else begin
        ;+ Read in the fits and save image of each one
        temp = dblarr(360, 180)
        for i=0, n_rot-1 do begin
            fxread, fit_dir+filenames[i], temp
            data[*, *, i] = temp
            savmap, temp, rotation_dir + string(rot[i]) + ".eps", 1.0, -10, 10
        endfor

        ;+ Save processed magnetic data for future reading
        save, data, filename=variables_path
    endelse

    ;+ Take a version of the maps containing positive flux only
    datapos = (data+abs(data)) / 2

    ;+
    ; Consider polar boundary to be 70 degrees
    ; (Use a large region because polar measurements are noisy)
    ; Consider emergence latitudes to be -40 to 40 degrees
    ;+
    polar_boundary = latind(70)
    emerge_boundary = latind(40)

    ;+
    ; Loop to flux calculations
    ; Total flux & longitude averaged flux
    ; Totals include signed and unsigned for the following intervals:
    ; global, 70-90, 40-70 degrees, 0-40 (in both hemispheres)
    ;-
    for i=0, n_rot-1 do begin
        ;+ Global
        totalf, data[*, *, i], signed=a, unsigned=b
        totalf, datapos[*, *, i], unsigned=c
        flux_s[i] = a
        flux_u[i] = b
        flux_p[i] = c
        flux_n[i] = b-c

        ;+ 70-90 degrees
        totalf, data[*, polar_boundary[1]:179, i], signed=a, unsigned=c ; north
        totalf, data[*, 0:polar_boundary[0], i], signed=b, unsigned=d ; south
        flux_pole_n[i] = a
        flux_pole_s[i] = b
        uflux_pole_n[i] = c
        uflux_pole_s[i] = d

        ;+ 40-70 degrees
        totalf, data[*, emerge_boundary[1]:polar_boundary[1], i], signed=a, unsigned=c ; north
        totalf, data[*, polar_boundary[0]:emerge_boundary[0], i], signed=b, unsigned=d ; south
        flux_mid_n[i] = a
        flux_mid_s[i] = b
        uflux_mid_n[i] = c
        uflux_mid_s[i] = d

        ; 0-40 degrees
        totalf, data[*, 90:emerge_boundary[1], i], signed=a, unsigned=c ; north
        totalf, data[*, emerge_boundary[0]:89, i], signed=b, unsigned=d ; south
        flux_active_n[i] = a
        flux_active_s[i] = b
        uflux_active_n[i] = c
        uflux_active_s[i] = d

        ;+ Longitude average (Butterfly diagram)
        avglong, data[*, *, i], signed=a
        flux_s_long[i, *] = a
    endfor

    ;+ Find areas to calculate averages over
    avgf, data[*, *, 0], area
    avgf, data[*, polar_boundary[1]:179, 0], parea
    avgf, data[*, emerge_boundary[1]:polar_boundary[1], 0], marea
    avgf, data[*, emerge_boundary[0]:89, 0], aarea

    ;+
    ; Active Longitudes
    ;
    ; It is known that flux tends to emerge along similar latitudes over
    ; extended periods of time. Calculations used to generate plots to inspect
    ; this property.
    ;-
    
    ;+ Split domain into 5 degree bins in longitude and latitude
    bin_deg = 5
    dim = bindim(bin_deg, 40.0)
    bin_flux = dblarr(dim[0], dim[1], n_rot)

    ;+ Calculate flux in 5 degree square bins per Carrington Rotation
    for i=0, n_rot-1 do begin
        bin_flux[*,*,i] = binf(data[*, *, i], bin_deg, boundary=40)
    endfor

    ;+ Prepare cycle total and cycle average arrays
    cycle_total = dblarr(dim[0], dim[1])
    cycle_avg = dblarr(dim[0], dim[1])

    ;+
    ; Hardcode approx start and end dates for Solar Cycles
    ; (March 1976-September 1986 and September 1986-August 1996)
    ; Dates converted using https://ssd.jpl.nasa.gov/tools/jdc/#/jd
    ;
    ; Dates of CRs chosen to separate each cycle
    ; CR 21: 1976-03-06 23:06:10 (unused)
    ; CR 22: 1986-09-16 12:10:52
    ; CR 23: 1996-08-22 03:07:51
    ;+
    c22 = (where(rot eq 1780))[0]
    c23 = (where(rot eq 1913))[0]

    ;+ Compute the total and average flux (per CR) in each bin across Cycle 22
    cycle_total[*, *] = total(bin_flux[*, *, c22:c23], 3)
    cycle_avg[*, *] =  cycle_total/double(c23-c22+1.0)
 
    ;+ Find number of CRs that occur for each yearly and biyearly interval
    yy_i = yy mod min(yy)
    yy2_i = (yy mod min(yy)) / 2

    ;+ Calculate how many occurances of each year and each two year period occurs in each interval
    annual_n = histogram(yy_i)
    biannual_n = histogram(yy2_i)

    ; Number of years and two year periods exist in each array
    annual_ntimes = n_elements(annual_n)
    biannual_ntimes = n_elements(biannual_n)

    ; Prepare arrays containing annual and biannual flux per bin
    annual_total = dblarr(dim[0], dim[1], n_elements(yy_i))
    annual_avg = dblarr(dim[0], dim[1], n_elements(yy_i))
    biannual_total = dblarr(dim[0], dim[1], n_elements(yy2_i))
    biannual_avg = dblarr(dim[0], dim[1], n_elements(yy2_i))

    ;+ Find average flux across bin in one/two year intervals
    for t=0, n_rot-1 do begin
        ;+ Yearly
        annual_total[*, *, yy_i[t]] = annual_total[*, *, yy_i[t]] + bin_flux[*, *, t]
        annual_avg[*, *, yy_i[t]] = annual_avg[*, *, yy_i[t]] + bin_flux[*, *, t]/annual_n[yy_i[t]]

        ;+ Biyearly
        biannual_total[*, *, yy2_i[t]] = biannual_total[*, *, yy2_i[t]] + bin_flux[*, *, t]
        biannual_avg[*, *, yy2_i[t]] = biannual_avg[*, *, yy2_i[t]] + bin_flux[*, *, t]/biannual_n[yy2_i[t]]
    endfor

    ;+
    ; End of active longitude computations
    ;-

    ;+
    ; Plots
    ;-

    ;+ To use custom fonts
    !P.FONT = 0

    ;+ Plos as postscript image
    set_plot,'ps'

    ; Global Unsigned Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Global Unsigned Flux.eps', /encapsulated
        plot, date, flux_u, xtitle="Year", ytitle="Total Magnetic Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/area, ticklen=0, ytitle="Average Flux Density (G)"
    device,/close

    ; Global Signed Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Global Signed Flux.eps', /encapsulated
        plot, date, flux_s, xtitle="Year", ytitle="Total Magnetic Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/area, ticklen=0, ytitle="Average Flux Density (G)"
    device,/close

    ; Global Signed Flux for each polarity
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Global Signed Flux Separated Polarities.eps', /encapsulated
        plot, date, flux_p, xtitle="Year", ytitle="Total Magnetic Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xtickinterval=10.0
        oplot, date, flux_n, linestyle=1, thick=2
    device,/close

    ; Pole Signed Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Polar Signed flux.eps', /encapsulated
        plot, date, flux_pole_n, xtitle="Year", ytitle="Total Polar Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        oplot, date, flux_pole_s, linestyle=1, thick=2
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/parea, ticklen=0, ytitle="Average Flux Density (G)"
    device,/close

    ; 40-70 degrees signed Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Mid Latitude Signed Flux.eps', /encapsulated
        plot, date, flux_mid_n, xtitle="Year", ytitle="Total Magnetic Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        oplot, date, flux_mid_s, linestyle=1, thick=2
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/marea, ticklen=0, ytitle="Average Flux Density (G)"
    device,/close

    ; 0-40 degrees signed Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Low Latitude Signed Flux.eps', /encapsulated
        plot, date, flux_active_n, xtitle="Year", ytitle="Total Magnetic Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        oplot, date, flux_active_s, linestyle=1, thick=2
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/aarea, ticklen=0, ytitle="Average Flux Density (G)"
    device,/close


    ; Pole unsigned Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Polar Unsigned Flux.eps', /encapsulated
        plot, date, uflux_pole_n, xtitle="Year", ytitle="Total Magnetic flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        oplot, date, uflux_pole_s, linestyle=1, thick=2
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/parea, ticklen=0, ytitle="Average Flux Density (G)"
    device, /close

    ; Mid-latitude unsigned Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Mid Latitude Unsigned Flux.eps', /encapsulated
        plot, date, uflux_mid_n, xtitle="Year", ytitle="Total Magnetic Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        oplot, date, uflux_mid_s, linestyle=1, thick=2
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/marea, ticklen=0, ytitle="Average Flux Density (G)"
    device, /close

    ; Active latitude unsigned Flux
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect1, filename=flux_dir+'Low Latitude Unsigned Flux.eps', /encapsulated
        plot, date, uflux_active_n, xtitle="Year", ytitle="Total Magnetic Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=10.0
        oplot, date, uflux_active_s, linestyle=1, thick=2
        axis, yaxis=1, yrange=(!Y.CRANGE)*1.0/aarea, ticklen=0, ytitle="Average Flux Density (G)"
    device, /close

    ;+ Butterfly diagram
    savbut, flux_s_long, butterfly_dir+"Longitude Averaged.eps", 1, -10, 10, date[0]
    
    ;+ Buttefly diagram which returns mflow estimates (interactive)
    ;dispbut, flux_s_long, 7, -5, 5, 100

    ;+ Annotated version used in dissertation
    savbut_annotated, flux_s_long, butterfly_dir+"Annotated Longitude Averaged.eps", 1, -5, 5, date[0]

    ; Plot total flux and average flux per bin each latitude during Solar Cycle 22
    for i=0, dim[1]-1 do begin
        device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect2, filename=cycle22_dir+strtrim(-abs(40.0)+bin_deg*i, 2)+".eps", /encapsulated
            plot, indgen(72)*bin_deg, cycle_total[*, i], xtitle="Carrington Longitude", ytitle="Total Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=60

            ;+ Plot the average flux per bin at that latitude
            oplot, indgen(72)*bin_deg, cycle_avg[*, i]
        device, /close
    endfor

    ;+ Scale factor to amplify the size of lines for each latitude in the active longitude plots
    scale=2.5

    ;+ Cycle 22 Southern Hemisphere
    fname = cycle22_dir+"Cycle 22 Southern Hemisphere.eps"
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect3, filename=fname, /encapsulated
        plot, indgen(72)*bin_deg, $ 
        -40 + scale*cycle_total[*, 0]/max(abs(cycle_total[*, *, *])), $
        xtitle="Carrington Longitude", ytitle="", xstyle=1,$
        thick=1, ystyle=1, xmargin=[4, 1], xtickinterval=60, $
        yrange=[-40, 0], ytickinterval=5.0, yminor=1
        oplot, !x.crange, fltarr(2)-40, linestyle=2
        for i=1, 7 do begin
            ;+ Plot the Total Flux across the year
            oplot, indgen(72)*bin_deg, (-40+5*i) + scale*cycle_total[*, i]/max(abs(cycle_total[*, *]))
            oplot, !x.crange, fltarr(2) + (-40+5*i), linestyle=2
        endfor
    device, /close

    ;+ Cycle 22 Southern hemisphere
    fname = cycle22_dir+"Cycle 22 Northern Hemisphere.eps"
    device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect3, filename=fname, /encapsulated
        plot, indgen(72)*bin_deg, $
        scale*cycle_total[*, 8]/max(abs(cycle_total[*, *, *])), $
        xtitle="Carrington Longitude", ytitle="", xstyle=1,$
        thick=1, ystyle=1, xmargin=[4, 1], xtickinterval=60, $
        yrange=[-2.5, 37.5], ytickinterval=5.0, yminor=1
        oplot, !x.crange, fltarr(2)-40, linestyle=2
        for i=9, 15 do begin
            ;+ Plot the Total Flux across the year
            oplot, indgen(72)*bin_deg, 5*(i-8) + scale*cycle_total[*, i]/max(abs(cycle_total[*, *]))
            oplot, !x.crange, fltarr(2) + 5*(i-8), linestyle=2
        endfor
    device, /close

    ; Plot individual plots for each latitude and year
    for i=0, dim[1]-1 do begin
        for t=0, annual_ntimes-1 do begin
            ;+ Make directory for that year if it doesn't exist
            new_dir = yearly_dir+strtrim(yy[0]+t, 2)+"/"
            file_mkdir, new_dir
            fname = new_dir+strtrim(-40+bin_deg*i, 2)+".eps"
            device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect2, filename=fname, /encapsulated
                ;+ Plot the Total Flux across the year
                plot, indgen(72)*bin_deg, annual_total[*, i, t], xtitle="Carrington Longitude", ytitle="Total Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=60

                ;+ Plot the average flux per bin at that latitude
                oplot, indgen(72)*bin_deg, annual_avg[*, i, t], linestyle=2
            device, /close
        endfor
    endfor

    ; Plot individual plots for each latitude in two year intervals
    for i=0, dim[1]-1 do begin
        for t=0, biannual_ntimes-1 do begin
            new_dir = biyearly_dir+strtrim(yy[0]+2*t, 2)+"/"
            file_mkdir, new_dir
            fname = new_dir+strtrim(-40+bin_deg*i, 2)+".eps"
            device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect2, filename=fname, /encapsulated
                ;+ Plot the Total Flux across the year
                plot, indgen(72)*bin_deg, biannual_total[*, i, t], xtitle="Carrington Longitude", ytitle="Total Flux (Mx)", xstyle=1, ystyle=ystyle, thick=2, xmargin=[10, 6], xtickinterval=60

                ;+ Plot the average flux per bin at that latitude
                oplot, indgen(72)*bin_deg, biannual_avg[*, i, t], linestyle=2
            device, /close
        endfor
    endfor

    ; Plot with combined latitude across one year intervals
    for a=0, 1 do begin
        label = ["South ", "North "]
        offset = [0, 8]
        ymin = [-40, 0]
        ymax = [0, 40]
        for t=0, annual_ntimes-1 do begin
            fname = yearly_dir + label[a] + strtrim(1975+t, 2) + ".eps"
            device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect3, filename=fname, /encapsulated
                plot, indgen(72)*bin_deg, ymin[a] + scale*annual_total[*, offset, t]/max(abs(annual_total[*, *, *])), $
                xtitle="Carrington Longitude", ytitle="", xstyle=1,$
                thick=1, ystyle=1, xmargin=[4, 1], xtickinterval=60, $
                yrange=[ymin[a]-2.5, ymax[a]-2.5], ytickinterval=5.0, yminor=1
                oplot, !x.crange, fltarr(2) + ymin[0], linestyle=2
                for i=1, 7 do begin
                    ;+ Plot the Total Flux across the year
                    oplot, indgen(72)*bin_deg, (ymin[a] + 5*i) + scale*annual_total[*, i+offset, t]/max(abs(annual_total[*, *, *]))
                    oplot, !x.crange, fltarr(2) + (ymin[a] + 5*i), linestyle=2
                endfor
            device, /close
        endfor
    endfor

    ; Plot with combined latitude across two year intervals
    for a=0, 1 do begin
        label = ["South ", "North "]
        offset = [0, 8]
        ymin = [-40, 0]
        ymax = [0, 40]
        for t=0, biannual_ntimes-1 do begin
            fname = biyearly_dir + label[a] + strtrim(1975 + 2*t, 2) + ".eps"
            device, /Times, font_size=fontsize, xsize=xsize, ysize=xsize*aspect3, filename=fname, /encapsulated
                plot, indgen(72)*bin_deg, ymin[a] + scale*biannual_total[*, offset, t]/max(abs(biannual_total[*, *, *])), $
                xtitle="Carrington Longitude", ytitle="", xstyle=1,$
                thick=1, ystyle=1, xmargin=[4, 1], xtickinterval=60, $
                yrange=[ymin[a]-2.5, ymax[a]-2.5], ytickinterval=5.0, yminor=1
                oplot, !x.crange, fltarr(2) + ymin[0], linestyle=2
                for i=1, 7 do begin
                    ;+ Plot the Total Flux across the year
                    oplot, indgen(72)*bin_deg, (ymin[a] + 5*i) + scale*biannual_total[*, i+offset, t]/max(abs(biannual_total[*, *, *]))
                    oplot, !x.crange, fltarr(2) + (ymin[a] + 5*i), linestyle=2
                endfor
            device, /close
        endfor
    endfor

    ;+ Plot Meridional flow (van Ballegooĳen et al., 1998)
    mflow, xsize, aspect1, 75.0, 11.0, fontsize, exportdir+'Meridional Flow.eps'

    ;+ Differential rotation (Snodgrass, 1983)
    diffrot, xsize, aspect1, fontsize, exportdir+'Differential Rotation.eps'

    ;+ Restore original tick length
    !P.TICKLEN = tick_orig
end
