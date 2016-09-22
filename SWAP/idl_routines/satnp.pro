n=0
dt=0.004

openr,1,'satnp.dat',/f77_unformatted
readu,1,n
print,n

satnp_arr=0
t=0.2
time=0.0

while (not(eof(1))) do begin
   readu,1,y 
   satnp_arr=[satnp_arr,y]
   t=t+dt
   time=[time,t]
endwhile

plot,time,satnp_arr,ytitle = 'np', $
     xtitle = 'release time (sec)',title='Barium density at satellite position'

close,1

end