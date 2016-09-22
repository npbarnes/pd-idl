PRO mom_m,nfiles,dt,POSTSCRIPT=postscript

;mom_m.pro
;plots momentum as a function of time for hybrid code


;file=''
close,1
openr,1,'momentum_1.dat',/f77_unformatted

nt=0ll
nout=0ll
ts=0ll
pup=dblarr(3,/nozero)
puf=dblarr(3,/nozero)
peb=dblarr(3,/nozero)
input_p=dblarr(3,/nozero)

readu,1,nt
readu,1,nout

print,nt,nout,nt/nout
readu,1,ts
tstep = ts
print,'time step #.....',ts
readu,1,pup,puf,peb,input_p

pupx=pup(0)
pupy=pup(1)
pupz=pup(2)

pufx=puf(0)
pufy=puf(1)
pufz=puf(2)

input_px = input_p(0)
input_py = input_p(1)
input_pz = input_p(2)

while (not(eof(1))) do begin
   readu,1,ts
   tstep = [tstep,ts]
;   print,'time step #.....',ts

   readu,1,pup,puf,peb,input_p

   pupx = [pupx,pup(0)]
   pupy = [pupy,pup(1)]
   pupz = [pupz,pup(2)]

   pufx = [pufx,puf(0)]
   pufy = [pufy,puf(1)]
   pufz = [pufz,puf(2)]

   input_px = [input_px, input_p(0)]
   input_py = [input_py, input_p(1)]
   input_pz = [input_pz, input_p(2)]

endwhile

close,1

nfile = 2
while (nfile le nfiles) do begin
file = 'momentum_'+strtrim(string(nfile),1)+'.dat'
openr,1,file,/f77_unformatted

;nt=0
;nout=0
;ts=0
;pup=dblarr(3,/nozero)
;puf=dblarr(3,/nozero)
;input_p=dblarr(3,/nozero)

readu,1,nt
readu,1,nout

print,nt,nout,nt/nout
;readu,1,ts
;tstep = ts
;print,'time step #.....',ts
;readu,1,pup,puf,input_p

;pupx=pup(0)
;pupy=pup(1)
;pupz=pup(2)

;pufx=puf(0)
;pufy=puf(1)
;pufz=puf(2)

;input_px = input_p(0)
;input_py = input_p(1)
;input_pz = input_p(2)

while (not(eof(1))) do begin
   readu,1,ts
   tstep = [tstep,ts]
;   print,'time step #.....',ts

   readu,1,pup,puf,peb,input_p

   pupx = [pupx,pup(0)]
   pupy = [pupy,pup(1)]
   pupz = [pupz,pup(2)]

   pufx = [pufx,puf(0)]
   pufy = [pufy,puf(1)]
   pufz = [pufz,puf(2)]

   input_px = [input_px, input_p(0)]
   input_py = [input_py, input_p(1)]
   input_pz = [input_pz, input_p(2)]

endwhile

close,1
nfile = nfile+1

endwhile

p_cons_m,nfiles,dt,tot_surf

if keyword_set(postscript) then begin
  set_plot,'ps
  !p.font=7
  device,/times
  device,filename='mom.eps'
  device,/encapsulated
  device,xsize=6.0,ysize=8.0,xoffset=1.0,yoffset=1.0,/inches
end


!p.multi=[0,3,5]
!y.title='p (kg m/s)'
!x.title='t (s)'
!x.ticks = 2
!x.style = 1
;!x.range = [0,2.5]
tstep = tstep
print,tstep
plot,tstep,pupx,title='(u!dp!n)!dx!n',charsize=1.5
plot,tstep,pupy,title='(u!dp!n)!dy!n',charsize=1.5
plot,tstep,pupz,title='(u!dp!n)!dz!n',charsize=1.5
plot,tstep,pufx,title='(u!df!n)!dx!n',charsize=1.5
plot,tstep,pufy,title='(u!df!n)!dy!n',charsize=1.5
plot,tstep,pufz,title='(u!df!n)!dz!n',charsize=1.5
plot,tstep,tot_surf(*,0),title='x surface',charsize=1.5
plot,tstep,tot_surf(*,1),title='y surface',charsize=1.5
plot,tstep,tot_surf(*,2),title='z surface',charsize=1.5
plot,tstep,input_px,title='input p!dx!n ',charsize=1.5
plot,tstep,input_py,title='input p!dy!n',charsize=1.5
plot,tstep,input_pz,title='input p!dz!n',charsize=1.5
;plot,tstep,(pufx+tot_surf(*,0))/input_py,title='x normalized',charsize=1.5
;plot,tstep,(pufy+tot_surf(*,1))/input_py,title='y normalized',charsize=1.5
;plot,tstep,(pufz+tot_surf(*,2))/input_py,title='z normalized',charsize=1.5


!y.title = ' '
normpx = (pupx+pufx-tot_surf(*,0))/input_py
plot,tstep,normpx,title='normalized p!dx!n',charsize=1.5
;	dx = !x.crange(1) - !x.crange(0)
;	dy = !y.crange(1) - !y.crange(0)
;	xyouts,0.4*dx+!x.crange(0),0.1*dy+!y.crange(0), $
; 'ave: '+ strmid(strtrim(string(total(normpx/n_elements(normpx))),1),0,3),/data

normpy = (pupy+pufy-tot_surf(*,1))/input_py
;normpy = (pupy+pufy)
plot,tstep,normpy,title='normalized p!dy!n',charsize=1.5
;	dx = !x.crange(1) - !x.crange(0)
;	dy = !y.crange(1) - !y.crange(0)
;	xyouts,0.4*dx+!x.crange(0),0.1*dy+!y.crange(0), $
; 'ave: '+ strmid(strtrim(string(total(normpy/n_elements(normpy))),1),0,3),/data
normpz = (pupz+pufz-tot_surf(*,2))/input_py
plot,tstep,normpz,title='normalized p!dz!n',charsize=1.5
;	dx = !x.crange(1) - !x.crange(0)
;	dy = !y.crange(1) - !y.crange(0)
;	xyouts,0.4*dx+!x.crange(0),0.1*dy+!y.crange(0), $
; 'ave: '+ strmid(strtrim(string(total(normpz/n_elements(normpz))),1),0,3),/data
;dx = !x.crange(1) - !x.crange(0)
;dy = !y.crange(1) - !y.crange(0)
;plots,[0.2*dx,0.4*dx]+!x.crange(0),[0.2*dy,0.2*dy]+!y.crange(0),linestyle=2, $
;      /data
;xyouts,0.42*dx+!x.crange(0),0.2*dy+!y.crange(0),'no b1z',/data

;xyouts,!d.x_size*0.6,!d.y_size*1.05,systime(0),/device

!p.multi=[0,1,1]

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
endif

return
end






