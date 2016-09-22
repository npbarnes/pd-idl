; read 3d coordinate data 
PRO plot_orbit,nfile,nfrm,POSTSCRIPT=postscript
;file='coord.dat'
f_read_coord,'coord.dat',qx,qy,qz,dzg,dzg,nx,ny,nz
close,1
;Ni_max=long(0)
nt=0ll
ntout=0ll
Ni_max=0ll
frm=0

if keyword_set(postscript) then begin
set_plot,'ps
!p.font=0
device,filename='ve_xe.eps'
device,/palatino
endif


xfile = 'xp_'+strmid(strtrim(nfile,2),0,1)+'.dat'
print,' reading...',xfile
openr,1,xfile,/f77_unformatted
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

print,x(3,0),x(3,1),x(3,2)


   
if keyword_set(postscript) then begin
device,/close
set_plot,'x'
!p.font=-1
endif

end
