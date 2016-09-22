; read 3d coordinate data 
PRO plot_xp,nfile,nfrm,xn,POSTSCRIPT=postscript
;file='coord.dat'
f_read_coord,'coord.dat',qqx,qqy,qqz,dzg,dzg,nx,ny,nz
close,1
Ni_max=long(0)
nt=0ll
ntout=0ll
frm=0ll

if keyword_set(postscript) then begin
set_plot,'ps
!p.font=0
device,filename='ve_xe.eps'
device,/palatino
endif

file ='xp_'+strmid(strtrim(nfile,2),0,1)+'.dat'
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

file ='xp_'+strmid(strtrim(nfile,2),0,1)+'.dat'
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

;xn = xn + t*vn

!p.multi=[0,1,1]

whx = where(xn(*,0) ne 0.0)
why = where(xn(*,1) ne 0.0)
whz = where(xn(*,2) ne 0.0)

szxn = size(xn)
whr = randomu(seed,10000)*szxn(1)

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

window,4,xsize=600,ysize=400
;!p.multi=[0,1,3]
;plot,xn(*,0),vn(*,0),psym=3
;plot,xn(*,0),vn(*,1),psym=3
;plot,xn(*,0),vn(*,2),psym=3 ;,xrange=[110,120],xstyle=1

sx=100
sz=20
nn = fltarr(sx,sz)

x = xn(whx,0)
y = xn(why,1)
z = xn(whz,2)

maxx = max(x)+1
minx = min(x)-1
xx = minx + findgen(sx)*(maxx-minx)/sx

maxz = max(z)+1
minz = min(z)-1
zz = minz + findgen(sz)*(maxz-minz)/sz

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
  nn(ii,jj) = nn(ii,jj) + 1.0
endfor

!p.multi=[0,1,1]
contour,smooth(nn,3),qx,qz,nlev=20
;,xrange=[minx,maxx],yrange=[minz,maxz],xstyle=1,ystyle=1
;window,5,xsize=800,ysize=200
;image_cont,rebin(nn,sx*5,sz*5)


if keyword_set(postscript) then begin
device,/close
set_plot,'x'
!p.font=-1
endif

end


