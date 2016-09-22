; read 3d coordinate data 
PRO ve_xe,nfile,nfrm,minj,maxj,POSTSCRIPT=postscript
;file='coord.dat'
f_read_coord,'coord.dat',qx,qy,qz,dzg,dzg,nx,ny,nz
close,1
Ni_max=long(0)
nt=0
ntout=0
frm=0

if keyword_set(postscript) then begin
set_plot,'ps
!p.font=0
device,filename='ve_xe.eps'
device,/palatino
endif

file ='ve_'+strmid(strtrim(nfile,2),0,1)+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
readu,1,nt
readu,1,ntout
readu,1,Ni_max
print,nt,ntout,Ni_max

v=fltarr(Ni_max,3,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,v
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,v
   frmcnt = frmcnt + 1

endwhile
close,1


xfile = 'xe_'+strmid(strtrim(nfile,2),0,1)+'.dat'
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

minv = min(v)
maxv = max(v)
minx = min(x)
maxx = max(x)

;plot,[minx,maxx],[minv,maxv],/nodata

wh = where( (x(*,1) gt qy(minj)) and (x(*,1) lt qy(maxj)))
plot,x(wh,2),v(wh,2),psym=3,xtitle='z (km)',ytitle='v!dz!n (km/s)',$
xrange=[min(x(wh,2)),max(x(wh,2))],xstyle=1
;for i = 0,n_elements(x)-1 do begin
;   plots,x(i),v(i)
;endfor
   
if keyword_set(postscript) then begin
device,/close
set_plot,'x'
!p.font=-1
endif

end
