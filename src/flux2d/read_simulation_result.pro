; Selects simulation folder
folder = DIALOG_PICKFILE(PATH='../src/flux_2d/examples', /DIRECTORY, $
   TITLE="Choose the root of the directory containing the simulation output.")
array = ''
line = ''

; Read in variables from structured file
; Order: nx, ny, xmax, xmin, ymax, ymin, timestamps(ascending)
get_lun, u
openr, u, (folder + 'simulation')
    while not eof(u) do begin
        readf, u, line
        array = [ array, line ]
    endwhile
close,u

; Convert from strings to floats
array = float(array)

; Give descriptive identifiers to array variables (excluding times)
time_index_displacement = 7
nx = array(1) & ny = array(2)
xmax = array(3) & xmin = array(4)
ymax = array(5) & ymin = array(6)
nt = n_elements(array) - time_index_displacement

; Set up quantity arrays
ax = dblarr(nx, ny+1, nt)
ay = dblarr(nx+1, ny, nt)
bz = dblarr(nx, ny, nt)
jx = dblarr(nx, ny+1, nt)
jy = dblarr(nx+1, ny, nt)

; Populate coordinate arrays
; Calculate dx
dx = (xmax-xmin)/nx
dy = (ymax-ymin)/ny

; Coordinates of Bz
xb = ((findgen(nx)*dx) + xmin + (dx/2.0))
yb = ((findgen(ny)*dy) + ymin + (dy/2.0))

; Coordinates of Ax
xax = ((findgen(nx)*dx) + xmin + (dx/2.0))
yax = ((findgen(ny+1)*dy) + ymin)

; Coordinates of Ay
xay = ((findgen(nx+1)*dx) + xmin)
yay = ((findgen(ny)*dy) + ymin + (dy/2.0))

; Coordinates of j are the same as that of A
xjx = xax
yjx = yax
xjy = xay
yjy = yay

; Read in data arrays
; Bz
openr, u, folder +'/bz.dat', /f77_unformatted
    for t=0, nt-1 do begin
        dummy = dblarr(nx, ny)
        readu, u, dummy
        bz(*,*,t) = dummy
    endfor
close, u

; Ax
openr, u, folder + 'ax.dat', /f77_unformatted
    for t=0, nt-1 do begin
        dummy = dblarr(nx, ny+1)
        readu, u, dummy
        ax(*,*,t) = dummy
    endfor
close, u

; Ay
openr, u, folder +'/ay.dat', /f77_unformatted
    for t=0, nt-1 do begin
        dummy = dblarr(nx+1, ny)
        readu, u, dummy
        ay(*,*,t) = dummy
    endfor
close, u

; jx
openr, u, folder +'/jx.dat', /f77_unformatted
    for t=0, nt-1 do begin
        dummy = dblarr(nx, ny+1)
        readu,u,dummy
        jx(*,*,t) = dummy
    endfor
close, u

; jy
openr, u, folder +'/jy.dat', /f77_unformatted
    for t=0, nt-1 do begin
        dummy = dblarr(nx+1, ny)
        readu,u,dummy
        jy(*,*,t) = dummy
    endfor
close, u
free_lun, u
end
