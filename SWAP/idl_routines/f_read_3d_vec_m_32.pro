; read 3d coordinate data 
PRO f_read_3d_vec_m_32,file,nfrm,x

nt=0l
nout=0l
nx=0l
ny=0l
nz=0l
frm=0l

files = file+'.dat'
openr,1,files,/f77_unformatted
print,'Reading file.....',files
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz
readu,1,frm
print,nt,nout,nx,ny,nz,frm
print,'  image #.....',frm

x=fltarr(nx,ny,nz,3,/nozero)

readu,1,x
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,x
   frmcnt = frmcnt + 1

endwhile

close,1

end

