vo=1.33e3
valf = 200e3
vsat = 9.6e3
nel = 1e12

t = 0.01

mba = 2.3e-25
mO = 26.0e-27
vol_cloud = (t*vo)^3
No = 1e25
tau = 23.0 ;s
Ni = No*(1-exp(-t/tau))
input_p = Ni*mba*vsat

d = 2*vo*t
z=valf*t

vol = 2*z*!pi*(d/2.0)^2

mom_trans = nel*vol*mO*vsat
print,'Momentum transfered....',t,nel*vol*mO*vsat,input_p

dt = 0.1
t = fltarr(21)
for i=1,20 do begin
   t(i)=t(i-1)+dt
   Ni = No*(1-exp(-t(i)/tau))
   input_p = [input_p,Ni*mba*vsat]
   d = 2*vo*t(i)
   z=valf*t(i)
   vol = 2*z*!pi*(d/2.0)^2
   mom_trans = [mom_trans,nel*vol*mO*vsat]
   print,'Momentum transfered....',t(i),mom_trans(i),input_p(i)
endfor

plot,t,input_p,linestyle = 1,ytitle='p!dx!n (kg m/s)',xtitle = 't (s)'
oplot,t,mom_trans,linestyle = 0
legend,['Input p!dx!n','Transfered p!dx!n'],linestyle=[1,0]


end
