; read 3d coordinate data 
PRO vp_dist,file,nfrm,x
;file='coord.dat'

Ni_max=long(0)
nt=0
ntout=0
frm=0

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

x=fltarr(Ni_max,3,/nozero)

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
