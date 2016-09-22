PRO pdp_diag_vec,file,nfiles,nfrm,ri,rj,rk,dt,POSTSCRIPT=postscript

;wavemovie.pro
;reads sequential data file and generates movie
;xyz specifies the component to view...0,1,2

ax = 0
ay = 0
az = 0
 
;t = findgen(nfiles*nfrm)*dt
t = 0

close,1

nt=0
nout=0
nx=0
ny=0
nz=0
frm=0

openr,1,file+'1.dat',/f77_unformatted

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

;rj = 17
;rk = 74

print,nt,nout,nx,ny,nz

arr = fltarr(nx,nz,nfiles*nfrm) 
temparr = fltarr(nx,ny,nz,3,/nozero)

close,1

frmcount=0
for m=1,nfiles do begin

   files = file+strtrim(string(m),1)+'.dat'
   openr,1,files,/f77_unformatted
   print,'Reading file.....',files

   readu,1,nt
   readu,1,nout
   readu,1,nx
   readu,1,ny
   readu,1,nz

   while (not(eof(1))) do begin
      readu,1,frm
      t = [t,frm*dt]
      print,'  image #.....',frm,frm*dt
      readu,1,temparr
      ax = [ax,temparr(ri,rj,rk,0)]
      ay = [ay,temparr(ri,rj,rk,1)]
      az = [az,temparr(ri,rj,rk,2)]
 
;      arr(*,*,frmcount) = reform(temparr(*,ny/2,*,xyz))
;;      surface,reform(temparr(*,ny/2,*,xyz))
;;      print,median(temparr(*,*,*,xyz))
;;      wait,.2
;;      print,median(temparr(*,*,*,xyz))
      frmcount = frmcount + 1
   endwhile

close,1

endfor

ax = ax(1:*)
ay = ay(1:*)
az = az(1:*)
t = t(1:*)
t0=0.01
t = t + t0
print,t

;window,1,xsize=500,ysize=700
if keyword_set(postscript) then begin
   set_plot,'ps
   !p.font=0
   device,/inches,xoffset=1.0,xsize=6.5,yoffset=1.0,ysize=9.0
;   device,/encapsulated
   device,filename='profile.ps
   device,/palatino
endif

!p.multi=[0,1,3]
!p.charsize=2.0
plot,t,ax,title='E!dx!n',xtitle='t (s)',ytitle='B!dx!n (s!u-1!n)'
plot,t,ay,title='E!dy!n',xtitle='t (s)',ytitle='B!dy!n (s!u-1!n)'
plot,t,az,title='E!dz!n',xtitle='t (s)',ytitle='B!dz!n (s!u-1!n)'

;print,ax,ay,az,t
!p.multi=[0,1,1]

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
endif

end




