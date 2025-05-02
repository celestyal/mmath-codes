;Directory containing the synoptic magnetogram observations
directory = '../data/mag/'

;Recovery carrington information from filenames, assuming they're correct
filenames = file_basename(file_search(directory+'*.fits'))

;Recover carrington rotation number from fit filenames
rot = stregex(filenames, '[0-9]+', /extract)
n_rot = n_elements(rot) ;number of carrington rotations

;Filename of the Carrington rotation metadata
carr_file = '../data/dates/carrington.rotations'

;Read in the other carrington metadata
carrdata, rot, carr_file, date=date, yy=yy

;Prepare arrays
data = dblarr(360,180,n_rot) 
flux_s = dblarr(n_rot)
flux_u = dblarr(n_rot)
flux_s_long = dblarr(n_rot,180)
flux_u_long = dblarr(n_rot,180)
aflux_u = dblarr(n_rot)
aflux_s = dblarr(n_rot)
aflux_npole = dblarr(n_rot)
aflux_spole = dblarr(n_rot)

;Get approx start and end dates for cycles 21 and 22
;(March 1976 to September 1986 and September 1986 to August 1996)
;Julian dates converted using https://ssd.jpl.nasa.gov/tools/jdc/#/jd
c21 = (where(rot eq 1639))(0) ;1976-03-06 23:06:10
c22 = (where(rot eq 1780))(0) ;1986-09-16 12:10:52
c23 = (where(rot eq 1913))(0) ;1996-08-22 03:07:51

;For polar flux, consider |latitudes| greater than 70
boundary = latind(70)

; Read in the fits
for i=0, n_rot-1 do begin
	fxread, directory+filenames(i), dummy
	data(*,*,i) = dummy

	;Total flux
	totalf, data(*,*,i), signed=a, unsigned=b
  flux_s(i) = a
  flux_u(i) = b

	;Polar flux
  totalf, data(*,boundary(1):179,i), signed=a
	totalf, data(*,0:boundary(0),i), signed=b
  aflux_npole(i) = a
  aflux_spole(i) = b

	;Average flux across the domain
	avgf, data(*,*,i), signed=a, unsigned=b
  aflux_s(i) = a
  aflux_u(i) = b

	;Average flux at the poles
	avgf, data(*, boundary(1):179, i), signed=a
	avgf, data(*, 0:boundary(0), i), signed=b
  aflux_npole(i) = a
  aflux_spole(i) = b

	;Longitude average of magnetic field)
	avglong, data(*,*,i), signed=a, unsigned=b
  flux_s_long(i,*) = a
  flux_u_long(i,*) = b
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Active Longitudes
;Consider latitudes between +-40 degrees
boundary = 40.0

; Calculate the flux in square bins at each time
bin_deg = 5
dim = bindim(bin_deg, boundary)
bin_flux = dblarr(dim(0), dim(1), n_rot)
for i=0, n_rot-1 do begin
  bin_flux(*,*,i) = binf(data, bin_deg, boundary=40)
endfor

;To find average flux in yearly, biyearly intervals
yy_i = yy mod min(yy)
yy2_i = (yy mod min(yy)) / 2

;Determine what each entry needs to be divided by to obtain average
annual_n = histogram(yy_i)
biannual_n = histogram(yy_i / 2)

annual = dblarr(dim(0), dim(1), n_elements(annual_n))
biannual = dblarr(dim(0), dim(1), n_elements(biannual_n))

;Average flux across 1 and 2 year
for t=0, n_elements(yy)-1 do begin
    annual(*, *, yy_i(t)) = annual(*, *, yy_i(t)) + bin_flux(*, *, t)/annual_n(yy_i(t))
    biannual(*, *, yy2_i(t)) = biannual(*, *, yy2_i(t)) + bin_flux(*, *, t)/biannual_n(yy2_i(t))
endfor

;Average flux across cycle 21 and 22
c21_avg = dblarr(dim(0), dim(1))
c22_avg = dblarr(dim(0), dim(1))

for t=c21, c22 do begin
    c21_avg(*, *) = c21_avg(*, *) + bin_flux(*, *, t)
endfor
c21_avg = c21_avg / (c22-c21)

;Sum flux across cycle 22
for t=c22, c23 do begin
    c22_avg(*, *) = c22_avg(*,*) + bin_flux(*,*,t)
endfor
c22_avg = c22_avg / (c23-c22)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plots
!P.FONT = 0
set_plot,'ps'
device, /Times, filename='total_unsigned_flux.pse', /encapsulated
plot, date(c21:c23), flux_u(c21:c23),xtitle="Year", ytitle="Magnetic Flux (Mx)", xstyle=1
device,/close

;Total signed Flux
device, /Times, filename='total_signed_flux.pse', /encapsulated
plot, date(c21:c23), flux_s(c21:c23), xtitle="Year", ytitle="Magnetic Flux (Mx)", xstyle=1
device,/close

;Average signed pole flux
device, /Times, filename='avg_polar_flux.pse', /encapsulated
plot, date(c21:c23), aflux_spole(c21:c23), xtitle="Year", ytitle="Average Magnetic Field (G)", xstyle=1
oplot, date(c21:c23), aflux_npole(c21:c23)
device,/close

;Longitude averaged magnetic field butteryfly diagram (Gauss)
scale=2
xy = size(flux_s_long(*,*), /DIMENSIONS)
xy = xy*scale
set_plot, 'z'
device, decompose=0, SET_RESOLUTION=xy
loadct, 0 ;colour table for the plot
;window, 0, xsize=xy(0), ysize=xy(1)
tv, bytscl(congrid(flux_s_long, xy(0), xy(1)), -10, 10)
write_png, "butterfly.png", TVRD(/TRUE)

;Plots for the bins WIP

print,"Saved graphics to IDL's current working directory."
end
