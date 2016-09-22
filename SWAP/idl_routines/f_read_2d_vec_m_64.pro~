; read 3d coordinate data 
PRO f_read_2d_vec_m,file,nfrm,x

nt=0ll
nout=0ll
nx=0ll
nz=0ll
frm=0ll

files = file+'.dat'
openr,1,files,/f77_unformatted
print,'Reading file.....',files
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,nz
readu,1,frm
print,nt,nout,nx,nz,frm
print,'  image #.....',frm

x=dblarr(nx,nz,3,/nozero)

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






