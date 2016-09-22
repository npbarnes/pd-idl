; read 3d coordinate data 
PRO ve_dist,path,nfile,nfrm,x,xx,OP=op,cl
;file='coord.dat'
f_read_coord,'coord.dat',qx,qy,qz
close,1
Ni_max=long(0)
nt=0
ntout=0
frm=0

pth = '/pf/d5/pad/pad/hybrid/'
pth = pth+path

file = pth+'ve_'+strmid(strtrim(nfile,2),0,1)+'.dat'
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

vv = sqrt(v(*,0)^2 + v(*,1)^2 + v(*,2)^2)

!p.multi=[0,1,1]

print,'maximum velocity...',max(vv)

wh = where((x(*,2) gt max(qz)*0.60) and (x(*,1) gt 0.4*max(qy)) $
            and (x(*,1) lt 0.6*max(qy)))
;wh = where(x(*,2) gt min(qz))

;binsize=0.01
;nv = fix((max(v(wh,0)) - min(v(wh,0)))/binsize)
;dv = (max(v(wh,0)) - min(v(wh,0)))/nv
;vel = min(v(wh,0)) + dv*findgen(nv)
;plot,vel,histogram(v(wh,0),min=min(v(wh,0)),max=max(v(wh,0))/2,bin=0.01)

;binsize=0.01
;nv = fix((max(v(wh,1)) - min(v(wh,1)))/binsize)
;dv = (max(v(wh,1)) - min(v(wh,1)))/nv
;vel = min(v(wh,1)) + dv*findgen(nv)
;plot,vel,histogram(v(wh,1),min=min(v(wh,1)),max=max(v(wh,1))/2,bin=0.01)

binsize=0.2
nv = fix((max(v(wh,2)) - min(v(wh,2)))/binsize)
dv = (max(v(wh,2)) - min(v(wh,2)))/nv
vel = min(v(wh,2)) + dv*findgen(nv)
v2 = sqrt(v(wh,2)^2 + v(wh,1)^2 + v(wh,0)^2)
vel2 = min(v2) + dv*findgen(nv)
h = histogram(v(wh,2),min=min(v(wh,2)),max=max(v(wh,2)),bin=binsize)
h2 = histogram(v2,min=min(v2),max=max(v2),bin=binsize)
if keyword_set(OP) then begin
;   oplot,vel,h,color=cl
   oplot,vel2,h2,color=cl
endif else begin
;   plot_io,vel,h>1.0,xrange=[-10,10],xtitle='v (km/s)'   
    plot,vel2,h2>1.0,xrange=[0,10],xtitle='v!dz!n (km/s)',ytitle='f(v)'
endelse

;;binsize=0.01
;nv = fix((max(vv) - min(vv))/binsize)
;dv = (max(vv) - min(vv))/nv
;vel = min(vv) + dv*findgen(nv)
;h = histogram(vv,min=min(vv),max=max(vv),bin=binsize)
;if keyword_set(OP) then begin
;   oplot,vel,h,color=100
;endif else begin
;   plot_io,vel,h>1.0
;endelse

close,1

end
