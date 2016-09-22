PRO regrid_xz,ar,xx,zx

@common

;convert irregularly gridded array to regular grid for viewing
;slices in the xz and yz planes.
;arr must be a 2d array

f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz

x=x(minx:maxx)
z=z(minz:maxz)
nx=n_elements(x)
nz=n_elements(z)

sz = size(ar)

xx=fltarr(nx,nz)
;yy=fltarr(ny,nz)
zx=fltarr(nx,nz)
;zy=fltarr(ny,nz)

print,size(xx),size(zx),size(ar)

for i=0,nz-1 do begin 
   xx(*,i)=x(*)
;   yy(*,i)=y(*)
endfor

for i=0,nx-1 do zx(i,*)=z(*)
;for i=0,ny-1 do zy(i,*)=z(*)
   
triangulate,xx,zx,trxz
;triangulate,yy,zy,tryz

gsxz=[(max(x)-min(x))/n_elements(x),(max(z)-min(z))/n_elements(z)]
;gsyz=[(max(y)-min(y))/n_elements(y),(max(z)-min(z))/n_elements(z)]

ar=trigrid(xx,zx,ar,trxz,gsxz)
;ar=congrid(ar,sz(1),sz(2))
sz=size(ar)
ar = ar(0:sz(1)-2,0:sz(2)-2)

print,size(xx),size(zx),size(ar)

return
end





