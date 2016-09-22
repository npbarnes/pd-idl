; read 3d coordinate data 
PRO f_read_2dxz_m_32,file,nfrm,x
;file='coord.dat'

nx=0l
ny=0l
nz=0l
nt=0l
ntout=0l
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,nx
readu,1,ny
readu,1,nz
print,nt,ntout,nx,ny,nz

x=fltarr(nx,nz,/nozero)

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
