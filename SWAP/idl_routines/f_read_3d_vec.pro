; read 3d coordinate data 

PRO f_read_3d_vec,file,x
;file='coord.dat'

nt=0
nout=0
nx=0
ny=0
nz=0
frm=0

openr,1,file,/f77_unformatted
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz
readu,1,frm
print,nt,nout,nx,ny,nz,frm

x=fltarr(nx,ny,nz,3,/nozero)

readu,1,x

close,1

end
