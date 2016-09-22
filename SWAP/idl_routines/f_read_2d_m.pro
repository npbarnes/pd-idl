; read 2d coordinate data 
PRO f_read_2d_m,file,nfrm,x
;file='coord.dat'

nx=0ll
nz=0ll
nt=0ll
ntout=0ll
frm=0ll

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,nx
readu,1,nz
print,nt,ntout,nx,nz

x=dblarr(nx,nz,/nozero)

readu,1,frm
print,'  image #.....',frm
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
