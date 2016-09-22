PRO regrid,arr

;convert irregularly gridded array to regular grid for viewing
;slices in the xz and yz planes.
;arr must be a 2d array

f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz

sz = size(arr)

xx=fltarr(nx,nz)
yy=fltarr(ny,nz)
zx=fltarr(nx,nz)
zy=fltarr(ny,nz)

for i=0,nz-1 do begin 
   xx(*,i)=x(*)
   yy(*,i)=y(*)
endfor

for i=0,nx-1 do zx(i,*)=z(*)
for i=0,ny-1 do zy(i,*)=z(*)
   
triangulate,xx,zx,trxz
triangulate,yy,zy,tryz

gsxz=[(max(x)-min(x))/n_elements(x),(max(z)-min(z))/n_elements(z)]
gsyz=[(max(y)-min(y))/n_elements(y),(max(z)-min(z))/n_elements(z)]

arr=trigrid(xx,zx,arr,trxz,gsxz)
arr=congrid(arr,sz(1),sz(2))

return
end





