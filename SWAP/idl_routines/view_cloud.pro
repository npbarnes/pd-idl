pro view_cloud,icldxz,icldxy,ncld

file = 'np.dat'

nx=0
ny=0
nz=0

openr,1,file,/f77_unformatted
readu,1,nx
readu,1,ny
readu,1,nz
print,nx,ny,nz

np=fltarr(nx,ny,nz,/nozero)

readu,1,np

close,1


file = 'nn.dat'

nx=0
ny=0
nz=0

openr,1,file,/f77_unformatted
readu,1,nx
readu,1,ny
readu,1,nz
print,nx,ny,nz

nn=fltarr(nx,ny,nz,/nozero)

readu,1,nn

close,1

icldxz=fltarr(nx,nz)
icldxy=fltarr(nx,ny)
ncld=fltarr(nx,nz)

for i=0,nx-1 do begin
   for k=0,nz-1 do begin
      icldxz(i,k) = total(np(i,*,k))
      ncld(i,k) = total(nn(i,*,k))
   endfor
endfor

for i=0,nx-1 do begin
   for j=0,ny-1 do begin
      icldxy(i,j) = total(np(i,j,*))
   endfor
endfor



end