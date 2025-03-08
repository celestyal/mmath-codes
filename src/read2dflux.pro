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
