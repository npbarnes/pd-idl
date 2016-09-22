vo=1.33e3
valf = 300e3
vsat = 9.6e3
nel = 1e12

t = 0.1

d = 2*vo*t
dt = d/vsat
z=valf*dt

vol = 2*z*!pi*(d/2.0)^2

num_e = nel*vol

print,'number of electrons....',num_e

for i=0,9 do begin
   t=t+dt
   d = 2*vo*t
   dt = d/vsat
   z=valf*dt
   vol = 2*z*!pi*(d/2.0)^2
   num_e = nel*vol
   print,'number of electrons....',t,num_e
endfor

end
