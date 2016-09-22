PRO mom,file,POSTSCRIPT=postscript

;mom.pro
;plots momentum as a function of time for hybrid code

if keyword_set(postscript) then begin
  set_plot,'ps
  !p.font=7
  device,/times
  device,filename='mom.ps'
  device,xsize=6.5,ysize=9.0,xoffset=1.0,yoffset=1.0,/inches
end

;file=''
close,1
openr,1,file,/f77_unformatted

nt=0
nout=0
ts=0
pup=fltarr(3,/nozero)
puf=fltarr(3,/nozero)
input_p=fltarr(3,/nozero)

readu,1,nt
readu,1,nout

print,nt,nout,nt/nout
readu,1,ts
tstep = ts
print,'time step #.....',ts
readu,1,pup,puf,input_p

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
   print,'time step #.....',ts

   readu,1,pup,puf,input_p

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

!p.multi=[0,3,4]
plot,tstep,pupx,title='pupx',charsize=1.5
plot,tstep,pupy,title='pupy',charsize=1.5
plot,tstep,pupz,title='pupz',charsize=1.5
plot,tstep,pufx,title='pufx',charsize=1.5
plot,tstep,pufy,title='pufy',charsize=1.5
plot,tstep,pufz,title='pufz',charsize=1.5
plot,tstep,input_px,title='input_px ',charsize=1.5
plot,tstep,input_py,title='input_py',charsize=1.5
plot,tstep,input_pz,title='input_pz',charsize=1.5
normpx = (pupx+pufx)/input_px
plot,tstep,normpx,title='normalized px',charsize=1.5
	dx = !x.crange(1) - !x.crange(0)
	dy = !y.crange(1) - !y.crange(0)
	xyouts,0.4*dx+!x.crange(0),0.1*dy+!y.crange(0), $
 'ave: '+ strmid(strtrim(string(total(normpx/n_elements(normpx))),1),0,3),/data

normpy = (pupy+pufy)/input_py
plot,tstep,normpy,title='normalized py',charsize=1.5
	dx = !x.crange(1) - !x.crange(0)
	dy = !y.crange(1) - !y.crange(0)
	xyouts,0.4*dx+!x.crange(0),0.1*dy+!y.crange(0), $
 'ave: '+ strmid(strtrim(string(total(normpy/n_elements(normpy))),1),0,3),/data
normpz = (pupz+pufz)/input_pz
plot,tstep,(pupz+pufz)/input_pz,title='normalized pz',charsize=1.5
	dx = !x.crange(1) - !x.crange(0)
	dy = !y.crange(1) - !y.crange(0)
	xyouts,0.4*dx+!x.crange(0),0.1*dy+!y.crange(0), $
 'ave: '+ strmid(strtrim(string(total(normpz/n_elements(normpz))),1),0,3),/data
;dx = !x.crange(1) - !x.crange(0)
;dy = !y.crange(1) - !y.crange(0)
;plots,[0.2*dx,0.4*dx]+!x.crange(0),[0.2*dy,0.2*dy]+!y.crange(0),linestyle=2, $
;      /data
;xyouts,0.42*dx+!x.crange(0),0.2*dy+!y.crange(0),'no b1z',/data

xyouts,!d.x_size*0.6,!d.y_size*1.05,systime(0),/device

!p.multi=[0,1,1]

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
endif

return
end






