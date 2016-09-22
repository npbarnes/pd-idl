PRO energy_m,nfiles,POSTSCRIPT=postscript

;energy.pro
;plots energy as a function of time for hybrid code

if keyword_set(postscript) then begin
set_plot,'ps
device,filename='energy.ps'
device,/palatino
device,/portrait
device,/encapsulated
device,/inches,xsize=6.0,xoffset=1.0,ysize=8.0,yoffset=1.0
end

file=''
file='energy_'
close,1

nt=0
nout=0
ts=0
tstep=0
input_E = 0.0
input_EeP = 0.0
Evp = 0.0
Euf = 0.0
EB1 = 0.0
EB1x = 0.0
EB1y = 0.0
EB1z = 0.0
EE = 0.0
EeP = 0.0
input_chex = 0.0
input_bill = 0.0

for m=1,nfiles do begin

   files = file+strtrim(string(m),1)+'.dat'

   openr,1,files,/f77_unformatted

   readu,1,nt
   readu,1,nout

   print,nt,nout,nt/nout

   while (not(eof(1))) do begin
      readu,1,ts
      tstep = [tstep,ts]
      print,'time step #.....',ts
      readu,1,d0,d05,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
      input_E = [input_E,d0]
      input_EeP = [input_EeP,d05]
      Evp = [Evp,d1]
      Euf = [Euf,d2]
      EB1 = [EB1,d3]
      EB1x = [EB1x,d4]
      EB1y = [EB1y,d5]
      EB1z = [EB1z,d6]
      EE = [EE,d7]
      EeP = [EeP,d8]
      input_chex = [input_chex,d9]
      input_bill = [input_bill,d10]
   endwhile

   close,1
endfor

tstep=tstep(1:*)*0.005
input_E=input_E(1:*)
input_EeP = input_EeP(1:*)
Evp=Evp(1:*)
Euf=Euf(1:*)
EB1=EB1(1:*)
EB1x=EB1x(1:*)
EB1y=EB1y(1:*)
EB1z=EB1z(1:*)
EE=EE(1:*)
EeP = EeP(1:*)
input_chex = input_chex(1:*)
input_bill = input_bill(1:*)

;total_E = Evp+Euf+EB1+EE+EeP
total_E = Evp+Euf+EB1
norm_E = total_E/(input_E)
print,'normalized energy...',tstep,norm_E

!p.multi=[0,2,3]
!p.font=0
!x.title = 't (s)'
!y.title = 'Energy (J)'
plot,tstep,Evp,title='v!dp!n',charsize=2.0
plot,tstep,Euf,title='u!de!n',charsize=2.0  ;electron energy
plot,tstep,EB1,title='b!d1!n',charsize=2.0
;plot,tstep,EB1x,title='b1x',charsize=2.0
;plot,tstep,EB1y,title='b1y',charsize=2.0
;plot,tstep,EB1z,title='b1z',charsize=2.0
;plot,tstep,EE,title='E',charsize=2.0
;plot,tstep,EeP,title='EeP',charsize=2.0
plot,tstep,total_E,title='total',charsize=2.0
;plot,tstep,total_E/(input_E+input_EeP),title='normalized (EeP)',charsize=2.0
plot,tstep,input_E,title = 'Input E',charsize=2.0
!y.title = ' '
plot,norm_E,title='normalized',charsize=2.0
;oplot,tstep,(Evp+Euf+EB1x+EB1y+EE)/input_E,linestyle=2
;dx = !x.crange(1) - !x.crange(0)
;dy = !y.crange(1) - !y.crange(0)
;plots,[0.2*dx,0.4*dx]+!x.crange(0),[0.2*dy,0.2*dy]+!y.crange(0),linestyle=2, $
;      /data
;xyouts,0.42*dx+!x.crange(0),0.2*dy+!y.crange(0),'no b1z',/data
;plot,tstep,input_chex/(input_E - input_chex - input_bill),$
;     title = 'Input Chex',charsize = 2.0 
;plot,tstep,input_bill/(input_E - input_chex - input_bill),$
;     title = 'Input Billiard',charsize = 2.0 

!p.multi=[0,1,1]
if keyword_set(postscript) then device,/close

return
end






