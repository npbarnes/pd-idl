; read 3d coordinate data 
PRO xneut,nfile,nfrm,t,POSTSCRIPT=postscript
;file='coord.dat'
f_read_coord,'coord.dat',qqx,qqy,qqz,dzg,dzg,nx,ny,nz
close,1
Ni_max=long(0)
nt=0ll
ntout=0ll
frm=0ll

;if keyword_set(postscript) then begin
;set_plot,'ps
;!p.font=0
;device,filename='ve_xe.eps'
;device,/palatino
;endif

file ='xneut_'+strmid(strtrim(nfile,2),0,1)+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

xn=dblarr(Ni_max,3,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xn
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xn
   frmcnt = frmcnt + 1

endwhile
close,1

file ='vneut_'+strmid(strtrim(nfile,2),0,1)+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

vn=dblarr(Ni_max,3,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,vn
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,vn
   frmcnt = frmcnt + 1

endwhile
close,1

xn = xn + t*vn

!p.multi=[0,1,1]

whx = where(xn(*,0) ne 0.0)
why = where(xn(*,1) ne 0.0)
whz = where(xn(*,2) ne 0.0)

szxn = size(xn)
whr = randomu(seed,szxn(1)/10)*szxn(1)

!p.charsize=0.5
window,1,xsize=300,ysize=250
plot,[min(xn(*,0)),max(xn(*,0))],[min(xn(*,2)),max(xn(*,2))],$
  xtitle='x (km)',ytitle='z (km)',xrange=[min(xn(whx,0)),max(xn(whx,0))], $
  yrange=[min(xn(whz,2)),max(xn(whz,2))],/nodata
plots,xn(whr,0),xn(whr,2),psym=3

window,2,xsize=300,ysize=250
plot,[min(xn(*,1)),max(xn(*,1))],[min(xn(*,2)),max(xn(*,2))],$
  xtitle='y (km)',ytitle='z (km)',xrange=[min(xn(why,1)),max(xn(why,1))], $
  yrange=[min(xn(whz,2)),max(xn(whz,2))],/nodata

xnx = xn(whr,0)
xny = xn(whr,1)
xnz = xn(whr,2)
wh = where((xnx gt 1) and (xnx lt 17))
;plots,xn(whr,1),xn(whr,2),psym=3
plots,xny(wh),xnz(wh),psym=3

window,3,xsize=300,ysize=250
plot,[min(xn(*,0)),max(xn(*,0))],[min(xn(*,1)),max(xn(*,1))],$
  xtitle='x (km)',ytitle='y (km)',xrange=[min(xn(whx,0)),max(xn(whx,0))], $
  yrange=[min(xn(why,1)),max(xn(why,1))],/nodata
plots,xn(whr,0),xn(whr,1),psym=3


;!p.multi=[0,1,3]
;plot,xn(*,0),vn(*,0),psym=3
;plot,xn(*,0),vn(*,1),psym=3
;plot,xn(*,0),vn(*,2),psym=3 ;,xrange=[110,120],xstyle=1

plot_xp,nfile,nfrm,xp

xp(*,2) = xp(*,2) - qqz(nz/2)
print,qqz(nz/2)
sx=100
sz=20
nn = fltarr(sx,sz)
np = fltarr(sx,sz)

x = xn(whx,0)
y = xn(why,1)
z = xn(whz,2)

z = z - qqz(nz/2)

maxxp = max(xp(*,0))+1
maxxn = max(x)+1
minxp = min(xp(*,0))-1
minxn = min(x)-1
if (maxxp gt maxxn) then maxx = maxxp
if (maxxp le maxxn) then maxx = maxxn
if (minxp lt minxn) then minx = minxp
if (minxp ge minxn) then minx = minxn
print,minx,maxx

xx = minx + findgen(sx)*(maxx-minx)/sx

maxzp = max(xp(*,2))+1
maxzn = max(z)+1
minzp = min(xp(*,2) > 0.0)-1
minzn = min(z)-1
if (maxzp gt maxzn) then maxz = maxzp
if (maxzp le maxzn) then maxz = maxzn
minz = minzn
;if (minzp lt minzn) then minz = minzp
;if (minzp ge minzn) then minz = minzn
;print,minz,maxz,minzp,minzn

zz = (minz + findgen(sz)*(maxz-minz)/(sz))

dx=xx(1) - xx(0)
dz=zz(1) - zz(0)

qx = fltarr(sx,sz)
qz = fltarr(sx,sz)

for i = 0,sx-1 do begin
  qx(i,*) = xx(i)
endfor
for j = 0,sz-1 do begin
  qz(*,j) = zz(j)
endfor

i = 0l

for i = 0l,(n_elements(whx)-1) do begin
  ii = round((x(i)-minx)/dx)
  jj = round((z(i)-minz)/dz)
;  print,ii,jj,z(i),maxz
  nn(ii,jj) = nn(ii,jj) + 1.0
endfor

qx = qx - 9.4*0.4 - 4.0*0.5


tm = (0.004*20)*(((nfile-1)*3.0) + nfrm)

xc = 9.4*(tm)
;zc = qqz((nz/2))
zc = 0.0
yc = (ny/2)+1
print,tm,xc,zc

r = 1.33*(tm+0.4)


xr = r - findgen(1001)*2*r/1000
zr = sqrt(r^2 - xr^2)

xrr = fltarr(n_elements(xr)*2)
zrr = fltarr(n_elements(xr)*2)
xrr = [xr,reverse(xr)]

zrr = [zr,-zr]

xrr = xrr+xc
zrr = zrr+zc

if keyword_set(postscript) then begin
set_plot,'ps
!p.charsize=1.0
device,filename='disk.ps
@x6x9
@pspl
endif
;window,4,xsize=600,ysize=600
!p.multi=[0,1,2]
!x.range=[0,15]
contour,smooth(nn,2),qx,qz,nlev=10,/isotropic,xtitle='x (km)',ytitle='z (km)',$
  title='Neutral Column Density'
plots,xrr,zrr,/data,thick=2
;,xrange=[minx,maxx],yrange=[minz,maxz],xstyle=1,ystyle=1
;window,5,xsize=800,ysize=200
;image_cont,rebin(nn,sx*5,sz*5)


for i = 0l,(n_elements(xp(*,0))-1) do begin
  ii = round((xp(i,0)-minx)/dx)
  jj = round((xp(i,2)-minz)/dz)
  if ((jj gt 0) and (jj lt sz)) then begin
    np(ii,jj) = np(ii,jj) + 1.0
  endif
endfor

contour,smooth(np,2),qx,qz,nlev=10,/isotropic,xtitle='x (km)',ytitle='z (km)',$
  title='Ion Column Density'
plots,xrr,zrr,/data,thick=2

device,filename='disk2.ps


tm = (0.004*20)*(((nfile-1)*3.0) + nfrm)
tm = tm + t

xc = 9.4*(tm)
;zc = qqz((nz/2))
zc = 0.0
yc = (ny/2)+1
print,tm,xc,zc

r = 1.33*(tm+0.4)


xr = r - findgen(1001)*2*r/1000
zr = sqrt(r^2 - xr^2)

xrr = fltarr(n_elements(xr)*2)
zrr = fltarr(n_elements(xr)*2)
xrr = [xr,reverse(xr)]

zrr = [zr,-zr]

xrr = xrr+xc
zrr = zrr+zc


!p.multi=[0,1,1]
!x.range=[0,max(qx)]
device,/inches,xoffset=1.0,xsize=6.5,yoffset=3.5,ysize=6.0
contour,smooth(nn,2),qx,qz,nlev=10,/isotropic,xtitle='x (km)',ytitle='z (km)',$
  xrange=[0,60],yrange=[-10,10]
plots,xrr,zrr,/data,thick=2


if keyword_set(postscript) then begin
device,/close
endif
set_plot,'x
!p.font=-1



;window,5
;!p.multi=[0,1,1]
;tvscl,rebin(nn,2*sx,2*sz)
write_jpeg,'nn.jpg',rebin(nn,3*sx,3*sz)*250/max(nn)

;if keyword_set(postscript) then begin
;device,/close
;set_plot,'x'
;!p.font=-1
;endif

end









