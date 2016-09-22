n=0
dt=0.004

openr,1,'chex.dat',/f77_unformatted
readu,1,n
print,n

chex_rate=0
t=0.0
time=0.0

while (not(eof(1))) do begin
   readu,1,y 
   chex_rate=[chex_rate,y]
   t=t+dt
   time=[time,t]
endwhile

plot,time,chex_rate,ytitle = 'dNi/dt', xtitle = 'release time (sec)',title='Charge exchange rates
xyouts,0.5,0.6,'% released material = ' + string(100.*total(chex_rate*dt)/2.6e24),/normal

print,100.*total(chex_rate*dt)/2.6e24
close,1


No=2.4e24
tau=28.0

dnidt = (No/tau)*exp(-time/tau)
oplot,time,dnidt,linestyle=1

end