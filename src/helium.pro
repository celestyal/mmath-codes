;--------------------------------------------------------------------
;  Display He I 10830 or magnetograph data:
;--------------------------------------------------------------------
pro display,file,nwin=nwin,bmax=bmax,fil=fil,helium=helium,  $
                 printer=printer,trans=trans
common blk2,file_dat
;
xold=-1.0+2.0*(0.5+findgen(180))/180   ; sine latitude of original data
xnew=sin(0.5*!pi*xold)   ; sine latitude of mapped data (uniform in angle)
index=fltarr(180)
index(0:5)=0.  &  index(174:179)=179.
i=0
for n=6,173 do begin
  while (xold(i+1) lt xnew(n)) do i=i+1
  frac=(xnew(n)-xold(i))/(xold(i+1)-xold(i)) > 0 < 1
  index(n)=float(i)+frac
endfor
;
fxread,file,data
for n=0,359 do data(n,*)=interpolate(reform(data(n,*),180),index)
if not keyword_set(bmax) then bmax=100.
data=(64./bmax)*(data+bmax)
data=byte((data > 0.) < 127.)
if keyword_set(helium) then begin
  data=data+128
  file_dat=strmid(file,0,strlen(file)-4)+'pnts'
endif
;
if keyword_set(printer) then begin
  set_plot,'PS'
  helium_color
  print,'Color copy on '+printer
  device,/landscape,/color,bits=8
  offset=3000
  xsize=18000
  ysize= 9000
  tv,[[[data]],[[data]],[[data]]],offset,offset,/device,  $
       xsize=xsize,ysize=ysize,true=3
endif else begin
  if not keyword_set(nwin) then begin
    if keyword_set(helium) then nwin=8 else nwin=4
  endif
  device,window_state=flag
  if flag(nwin) then wset,nwin else window,nwin,xsize=800,ysize=440
  erase
  offset=50
  xsize=720
  ysize=360
  tv,rebin(data,xsize,ysize),offset,offset
endelse
xtitle='LONGITUDE IN DEGREES'
ytitle='LATITUDE IN DEGREES'
title=string(file,format='("FILE = ",a)')
px=float([offset,offset+xsize])/!d.x_size
py=float([offset,offset+ysize])/!d.y_size
plot,[0,360],[0,0],/nodata,  $
         position=[px(0),py(0),px(1),py(1)],/noerase,       $
         xrange=[0,360],xstyle=1,                           $
         yrange=[-90,90],ystyle=1,                          $
         xticks=4,xtickv=[0,90,180,270,360],xminor=9,       $
         yticks=2,ytickv=[-90,0,90]        ,yminor=9,       $
         xtitle=xtitle,ytitle=ytitle,title=title
if keyword_set(fil) then begin
  if keyword_set(helium) then color=255 else color=0
  filaments,color=color
endif
if keyword_set(printer) then begin
  device,/close
  if keyword_set(trans) then str=' -SInputSlot=Transparency' else str=' '
  spawn,'pslpr -d'+printer+str+' idl.ps'
  set_plot,'X'
  helium_color
endif
end

;--------------------------------------------------------------------
;  Display filaments:
;--------------------------------------------------------------------
pro filaments,color=color
common blk2,file_dat
openr,1,file_dat,error=error
if not keyword_set(color) then color=0
if error eq 0 then begin
  n=0
  while not eof(1) do begin
    readf,1,n
    lon=intarr(n)  &  readf,1,lon
    lat=intarr(n)  &  readf,1,lat
    oplot,lon,lat,color=color
  endwhile
endif
close,1
end

;--------------------------------------------------------------------
;  Manually select filaments:
;--------------------------------------------------------------------
pro select
common blk2,file_dat
print,'Click LEFT to select point, MIDDLE to end filament, RIGHT to exit'
while 1 do begin
  lon=intarr(100)
  lat=intarr(100)
  n=0
  cursor,i,j,3,/device  &  button=!err
  while button eq 1 and n lt 100 do begin
    if n eq 0 then print,'Start filament:'
    lon(n)=(i-50)/2
    lat(n)=(j-50)/2-90
    print,'LON=',strcompress(lon(n)),' LAT=',strcompress(lat(n))
    if n gt 0 then plots,[lon(n-1),lon(n)],[lat(n-1),lat(n)],/data,color=255
    n=n+1
    cursor,i,j,3,/device  &  button=!err
  endwhile
  if n ge 2 then begin
    openw,1,file_dat,/append
    printf,1,n
    printf,1,strcompress(lon(0:n-1))
    printf,1,strcompress(lat(0:n-1))
    close,1
  endif
  if button eq 4 then return
endwhile
end

;--------------------------------------------------------------------
;  Define color table for magnetograms and He I 10830:
;--------------------------------------------------------------------
pro helium_color,gamma=gamma
;
;  BLUE-RED color table for magnetograms (replace first and last
;  entries with BLACK and WHITE):
;
loadct,33,/silent
tvlct,r,g,b,/get
m=!d.n_colors
mh=m/2
r1=rebin(r(0:2*mh-1),mh)
g1=rebin(g(0:2*mh-1),mh)
b1=rebin(b(0:2*mh-1),mh)
r1(0)=0       &  g1(0)=0       &  b1(0)=0        ; black
r1(mh-1)=255  &  g1(mh-1)=255  &  b1(mh-1)=255   ; white
;
;  RED color table for He I 10830:
;
if not keyword_set(gamma) then gamma=1.0
x=(findgen(m-mh)/(m-mh-1))^gamma
r2=byte(255.*(          x/0.66 < 1.))
g2=byte(255.*(1.0+(x-1.0)/0.52 > 0.))
b2=byte(255.*(1.0+(x-1.0)/0.25 > 0.))
;
;  Combine color tables:
;
r=[r1,r2]  &  g=[g1,g2]  &  b=[b1,b2]
tvlct,r,g,b
end

common blk2,file_dat
end

