PRO energy,file

;energy.pro
;plots energy as a function of time for hybrid code

;file=''
;file='energy.dat'
close,1
openr,1,file,/f77_unformatted

nt=0
nout=0
ts=0
readu,1,nt
readu,1,nout

print,nt,nout,nt/nout
readu,1,ts
tstep = ts
print,'time step #.....',ts
readu,1,d0,d05,d1,d2,d3,d4,d5,d6,d7,d8
print,d0,d05,d1,d2,d3,d4,d5,d6,d7,d8
input_E = d0
input_EeP = d05
Evp = d1
Euf = d2
EB1 = d3
EB1x = d4
EB1y = d5
EB1z = d6
EE = d7
EeP = d8

while (not(eof(1))) do begin
   readu,1,ts
   tstep = [tstep,ts]
   print,'time step #.....',ts
   readu,1,d0,d1,d2,d3,d4,d5,d6,d7
   input_E = [input_E,d0]
   Evp = [Evp,d1]
   Euf = [Euf,d2]
   EB1 = [EB1,d3]
   EB1x = [EB1x,d4]
   EB1y = [EB1y,d5]
   EB1z = [EB1z,d6]
   EE = [EE,d7]
endwhile

close,1

total_E = Evp+Euf+EB1+EE
norm_E = total_E/input_E
print,'normalized energy...',tstep,norm_E

!p.multi=[0,3,3]
plot,tstep,Evp,title='vp',charsize=2.0
plot,tstep,Euf,title='uf',charsize=2.0
plot,tstep,EB1,title='b1',charsize=2.0
plot,tstep,EB1x,title='b1x',charsize=2.0
plot,tstep,EB1y,title='b1y',charsize=2.0
plot,tstep,EB1z,title='b1z',charsize=2.0
plot,tstep,EE,title='E',charsize=2.0
plot,tstep,total_E,title='total',charsize=2.0
plot,tstep,norm_E,title='normalized',charsize=2.0
oplot,tstep,(Evp+Euf+EB1x+EB1y+EE)/input_E,linestyle=2
dx = !x.crange(1) - !x.crange(0)
dy = !y.crange(1) - !y.crange(0)
plots,[0.2*dx,0.4*dx]+!x.crange(0),[0.2*dy,0.2*dy]+!y.crange(0),linestyle=2, $
      /data
xyouts,0.42*dx+!x.crange(0),0.2*dy+!y.crange(0),'no b1z',/data

!p.multi=[0,1,1]

return
end






