pro set_params
    common params, fontsize, xsize, aspect
    fontsize = 11.0
    xsize = 8.0 ; centimetres
    aspect = 0.75
    !P.TICKLEN = -abs(!P.TICKLEN)
end

pro read_cft
    ; Define all variables as common so they can be called by other procedures
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin

    ; Selects simulation folder
    folder = DIALOG_PICKFILE(PATH='./', $
        TITLE="Select the folder containing the exported simulation", $
        /DIRECTORY, $
        /MUST_EXIST)

    file = folder+'simulation.bin'

    ; Preamble (we're not interested in them here)
    version = long64(0)
    hlen = long64(0)

    ; Header items
    solv = long64(0)
    bc = long64(0)
    nx = long64(0)
    ny = long64(0)
    sav = long64(0)
    d = double(0.0)
    t = double(0.0)
    lx = double(0.0)
    ly = double(0.0)
    xmin = double(0.0)
    xmax = double(0.0)
    ymin = double(0.0)
    ymax = double(0.0)

    ; Read in file contents in order
    get_lun, u
    openr, u, file
    ; Preamble
    readu, u, version, hlen
    ; Header
    readu, u, solv, bc, nx, ny, sav, d, t, lx, ly, xmin, xmax, ymin, ymax

    ; Allocate arrays using header info
    vx = dblarr(nx+1,ny+1)
    vy = dblarr(nx+1,ny+1)
    times = dblarr(sav)
    ax = dblarr(sav, nx, ny+1)
    ay = dblarr(sav, nx+1, ny)

    ; Velocities
    readu, u, vx, vy

    ; Saves, each of which contains time, ax and ay
    x = double(0.0)
    temp1 = dblarr(nx, ny+1)
    temp2 = dblarr(nx+1, ny)
    for i=0, fix(sav)-1 do begin
            readu, u, x, temp1, temp2
            times(i) = x
            ax(i, *, *) = temp1
            ay(i, *, *) = temp2
    endfor
    close, u

    ; Set up arrays for other quantities of interest that are derived from A
    ;bz = dblarr(sav, nx, ny)
    bz = dblarr(sav, nx+2, ny+2)
    jx = dblarr(sav, nx, ny+1)
    jy = dblarr(sav, nx+1, ny)

    ;
    ; Calculate Bz using finite difference
    ;
    dx = (xmax-xmin)/nx
    dy = (ymax-ymin)/ny

    ; Dummy arrays
    dummy = dblarr(nx, ny)

    ; Read in each Bz
    for i=0, fix(sav)-1 do begin
        dummy = (ay(i, 1:nx, 0:ny-1) - ay(i, 0:nx-1, 0:ny-1))/dx - (ax(i, 0:nx-1, 1:ny)-ax(i, 0:nx-1, 0:ny-1))/dy
        bz(i,1:nx,1:ny) = dummy

        ; Apply BC
        case bc of
            ; Periodic
            1: begin
                bz(i, 0, *) = bz(i, nx, *)
                bz(i, nx+1, *) = bz(i, 1, *)
            end
            ; Infinite
            2: begin
                bz(i, 0, *) = bz(0, 1, *)
                bz(i, nx+1, *) = bz(0, nx, *)
            end
            ; Open
            3: begin
                bz(i, 0, *) = bz(i, 1, *)
                bz(i, nx+1, *) = bz(i, nx, *)
            end
            else: print, "Unsure how to treat this boundary condition"
        endcase
	; Apply closed boundary at horizontal edges
	bz(i, *, 0) = bz(i, *, 1)
	bz(i, *, ny+1) = bz(i, *, ny)
    endfor

    ; Coordinates of cell corners (and vx, vy)
    x = dindgen(nx+1)
    y = dindgen(ny+1)
    x = xmin + x*dx
    y = ymin + y*dy

    ; Coordinates of Bz
    xb = x(0:nx-1)+dx/2.0
    yb = y(0:ny-1)+dy/2.0

    ; Coordinates of Ax and jx
    xax = xb
    yax = y

    ; Coordinates of Ay and jy
    xay = x
    yay = yb
end


pro bz_ycut_multi_t, yslice, t_ind, linestyles, ylim=ylim
    ; Define common variables needed from read_cft
     common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    old_dev = !D.NAME
    set_plot, 'ps'
    !P.FONT = 0

    repeats = size(t_ind, /n_elements)-1

    ; Plot Bz
    device, /Times, xsize=xsize, ysize=xsize*aspect, font_size=fontsize, filename=folder+'bz_cut.eps', /encapsulated
        if keyword_set(ylim) then begin
            plot, xb, bz(t_ind(0), 1:nx, yslice), xtitle="x (non-dimensional)", ytitle="B!Dz!N Strength (G)", xstyle=1, ystyle=3, thick=2, linestyle=linestyles(0), yrange=ylim
            for k=1, repeats do begin
                oplot, xb, bz(t_ind(k), 1:nx, yslice), thick=2, linestyle=linestyles(k)
            endfor
        endif else begin
            plot, xb, bz(t_ind(0), 1:nx, yslice), xtitle="x (non-dimensional)", ytitle="B!Dz!N Strength (G)", xstyle=1, ystyle=3, thick=2, linestyle=linestyles(0)
            for k=1, repeats do begin
                oplot, xb, bz(t_ind(k), 1:nx, yslice), thick=2, linestyle=linestyles(k)
            endfor
        endelse
    device, /close

    print, 'Plotted cut (constant y) of Bz at times (days): ', times(t_ind)/84000.0

    set_plot, old_dev
end


pro current_ycut_multi_t, yslice, t_ind, linestyles
    ; Define common variables needed from read_cft
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    ; jx and jy (curl B)
    jx = dblarr(sav, nx, ny+1)
    jy = dblarr(sav, nx+1, ny)
    jx =  (bz(*, 1:nx, 1:*)-bz(*, 1:nx, 0:ny)) / dy
    jy = -(bz(*, 1:*, 1:ny)-bz(*, 0:nx, 1:ny)) / dx

    old_dev = !D.NAME
    set_plot, 'ps'
    !P.FONT = 0

    repeats = size(t_ind, /n_elements)-1

    device, /Times, xsize=xsize, ysize=xsize*aspect, filename=folder + 'jx.eps', /encapsulated, font_size=fontsize
        plot, xax, jx(t_ind(0), *, yslice) ,xtitle='x (non-dimensional)',ytitle='j!Ix', xstyle=1, ystyle=1, linestyle=5, thick=2
        for k=1, repeats do begin
            oplot, xay, jx(t_ind(k), *, yslice), thick=2
        endfor
    device, /close

    device, /Times, xsize=xsize, ysize=xsize*aspect, filename=folder +'jy.eps', /encapsulated, font_size=fontsize
        plot, xay, jy(t_ind(0), *, yslice), xtitle='x (non-dimensional)', ytitle='j!Ix', xstyle=1, ystyle=1, linestyle=linestyles(0), thick=2, yrange=[-200.0,1.0]
        for k=1, repeats do begin
            oplot, xay, jy(t_ind(k), *, yslice), thick=2, linestyle=linestyles(k)
        endfor
    device, /close

    print, 'Plotted cut (constant y) of jx and jy at times (seconds): ', times(t_ind)

    set_plot, old_dev
end


pro current_ycut_integral, yslice, ylim=ylim
    ; Define common variables needed from read_cft
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    ; Integral of jy along a 1D slice over time
    integral = dblarr(sav)
    for t=0, sav-1 do begin
        integral(t) = 0
        i = 1
        integral(t) = dx*(abs(jy(t, 0, yslice)) + abs(jy(t, nx-1, yslice)))*0.5
        while (i lt nx-2) do begin
            integral(t) = integral(t) + dx*abs(jy(t, i, yslice))
            i = i + 1
        endwhile
    endfor

    old_dev = !D.NAME
    set_plot, 'ps'
    !P.FONT = 0

    device, /Times, xsize=xsize, ysize=xsize*aspect, filename=folder + 'jy_integral.eps', /encapsulated, font_size=fontsize
        if keyword_set(ylim) then begin
            plot, times/86400.0, 20-integral, xtitle='Time (days)', ytitle='Total |j!Dy!N| Absolute Error', xstyle=1, ystyle=1, linestyle=0, thick=2, yrange=ylim
        endif else begin
            plot, times/86400.0, 20-integral, xtitle='Time (days)', ytitle='Total |j!Dy!N| Absolute Error', xstyle=1, ystyle=1, linestyle=0, thick=2
        endelse
    device, /close

    print, 'Plotted cut (constant y) of the integral of jy.'

    set_plot, old_dev
end


; Note: this should be replaced with the IDL code when the updated IDL branch is merged to main
pro savmapc, data, filename, scale, low, high, xmin, xmax, ymin, ymax
    common params, fontsize, xsize, aspect
    ; Save the current device type to return to after saving
    prevdev = !d.name

    ; Get dimensions of the image and scale
    xy = size(data, /DIMENSIONS)
    xy = xy*scale

    ; Set device to postscript
    set_plot, 'ps'

    xsize=8.0
    ; Take negative of the data and load the colour table (so the colour table is inverted)
    loadct, 71, /silent

    device, decomposed=0, xsize=xsize, /encapsulated
        ; Create and save image
        tvimg = bytscl(congrid(data, xy(0), xy(1)), low, high)
        tvlct, r, g, b, /GET
        imageRGB = bytarr(3, xy(0), xy(1), /NOZERO)
        imageRGB[0, *, *] = r[tvimg]
        imageRGB[1, *, *] = g[tvimg]
        imageRGB[2, *, *] = b[tvimg]
        im = image(imageRGB, /buffer, margin=[0.19, 0.19, 0.03, 0.0], image_dimensions=[xy(0), xy(1)], dimensions=xy+64)
        xax = axis('X', location='bottom', TICKDIR=1, MINOR=0, title='x (non-dimensional)', coord_transform=[xmin, (xmax-xmin)/fix(xy(0), type=4)], tickfont_name='Times', tickfont_size=26)
        yax = axis('Y', location='left', TICKDIR=1, MINOR=0, coord_transform=[ymin, (ymax-ymin)/fix(xy(1), type=4)], title='y (non-dimensional)', tickfont_name='Times', tickfont_size=26)
        im.Save, filename, /transparent, border=0
    device, /close

    ; Return to previous device
    set_plot, prevdev
end



pro example01
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    set_params
    read_cft

    bz_ycut_multi_t, floor(double(ny)/2.0), [0, ceil(sav*0.01), ceil(sav*0.2), sav-1], [3, 2, 1, 0]

    current_ycut_multi_t, floor(double(ny)/2.0), [0, ceil(sav*0.2), sav-1], [2, 1, 0]

    current_ycut_integral, floor(double(ny)/2.0)
end


pro example02
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    set_params
    read_cft

    plt_times = [0, ceil(sav*0.01), ceil(sav*0.2), sav-1]
    line_styles = [3, 2, 1, 0]
    yslice = floor(double(ny)/2.0)
    bz_ycut_multi_t, yslice, plt_times, line_styles
    
    ;+ 
    ; Plot analytical solution at same timestamps as simulation
    ;-
    old_dev = !D.NAME
    set_plot, 'ps'
    !P.FONT = 0
    solution = dblarr(size(xb, /n_dimensions))
    device, /Times, xsize=xsize, ysize=xsize*aspect, font_size=fontsize, filename=folder+'bz_analytical.eps', /encapsulated
        ;+ Solution at t=0 (argument to approximate infinity)
        solution = 10*erf(2147483647*xb)
        plot, xb, solution, xtitle="x (non-dimensional)", ytitle="B!Dz!N Strength (G)", xstyle=1, ystyle=3, thick=2, linestyle=line_styles(0)
	    for k=1, size(plt_times, /n_elements)-1 do begin
            ;+ Plots using dimensional of quantities
            solution = 10*erf((xb*lx)/sqrt(4*d*times(plt_times(k))))
            oplot, xb, solution, thick=2, linestyle=line_styles(k)
        endfor
    device, /close
    set_plot, old_dev

    ;+
    ; Plot the error of the two intermediate curves to get a sense of the error.
    ; Plotted separately because the errors vary by an order of magnitude.
    ;-
    old_dev = !D.NAME
    set_plot, 'ps'
    !P.FONT = 0

    ;+ Error for small t
    error = dblarr(size(xb, /n_dimensions))
    device, /Times, xsize=xsize, ysize=xsize*aspect, font_size=fontsize, filename=folder+'bz_error1.eps', /encapsulated
        ;+ Plot penultimate error (higher)
        solution = 10*erf((xb*lx)/sqrt(4*d*times(plt_times(1))))
        plot, xb, abs(solution-bz(plt_times(1), 1:nx, yslice)), xtitle="x (non-dimensional)", ytitle="Absolute Error", xstyle=1, ystyle=3, thick=2, linestyle=line_styles(1)
    device, /close

    ;+ Error for large t
    device, /Times, xsize=xsize, ysize=xsize*aspect, font_size=fontsize, filename=folder+'bz_error2.eps', /encapsulated
        ;+ Plot penultimate error (higher)
        solution = 10*erf((xb*lx)/sqrt(4*d*times(plt_times(2))))
        plot, xb, abs(solution-bz(plt_times(2), 1:nx, yslice)), xtitle="x (non-dimensional)", ytitle="Absolute Error", xstyle=1, ystyle=3, thick=2, linestyle=line_styles(2)
    device, /close

    set_plot, old_dev


    current_ycut_multi_t, floor(double(ny)/2.0), [0, ceil(sav*0.2), sav-1], [2, 1, 0]
    current_ycut_integral, floor(double(ny)/2.0) ;, ylim=[17.0, 23.0]
end


pro example03
    ; Same plots as Example 1
    example01
end


pro example04
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    set_params
    read_cft

    bz_ycut_multi_t, floor(double(ny)/2.0), [sav-1], [0]
end


pro example05
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    set_params
    read_cft

    bz_ycut_multi_t, floor(double(ny)/2.0), [0, sav-1], [1, 0], ylim=[9.5, 10.5]
   
    ; Plot surface at end of simulation to check stability across whole domain
    old_dev = !D.NAME
    set_plot, 'ps'
    !P.FONT = 0
    device, /Times, xsize=xsize, ysize=xsize*aspect, font_size=fontsize, filename=folder +'surface.eps', /encapsulated
            surface, bz(sav-1, *, *), zrange=[9.0,11.0], xstyle=4, ystyle=4, /upper_only
    device, /close
    set_plot, old_dev
end


pro example06
    ; Same plots as Example 5
    example05
end


pro example07
    ; Same plots as Example 5
    example05
end


pro example08
    common sim, ax, ay, bc, bz, d, dx, dy, file, folder, hlen, jx, jy, lx, ly, nx, ny, sav, solv, times, version, vx, vy, x, xax, xay, xb, xmax, xmin, y, ymax, ymin
    common params, fontsize, xsize, aspect

    set_params
    read_cft

    ;t_ind = [0, ceil(sav*0.5), sav-1]
    t_ind = [0, ceil(sav*0.1), sav-1]

    data = dblarr(nx, ny)

    for k=0, 2 do begin
        data(*, *) = bz(t_ind(k), 1:nx, 1:ny)
        savmapc, data, folder+'image-'+string(k)+'.eps', 1.0, -10.0, 10.0, xmin, xmax, ymin, ymax
    endfor

    bz_ycut_multi_t, floor(double(ny)/2.0), t_ind, [2, 1, 0]
end


pro example09
    ; Same plots as Example 8
    example08
end


pro example10
    ; Same plots as Example 8
    example08
end
