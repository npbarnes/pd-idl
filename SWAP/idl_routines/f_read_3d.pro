; read 3d coordinate data 

PRO f_read_3d,file,x
;file='coord.dat'

nx=0
ny=0
nz=0
nt=0
ntout=0

openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,nx
readu,1,ny
readu,1,nz
print,nt,ntout,nx,ny,nz

x=fltarr(nx,ny,nz,/nozero)

readu,1,m
readu,1,x

close,1

end
