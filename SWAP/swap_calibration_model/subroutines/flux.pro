function  flux,  phi, theta, vr,theta1,phi1
close, 110
openr, 110,'inputs.txt'
readf, 110,  format='(3F30.8)', n, vo, t1
close, 110
;print, n, vo, t1
mp=1.67262158D-27
kb=1.380658D-23
vth=sqrt((2.0D)*kb*t1/mp)*100.d

amp=n/(vth^3)/(!dpi)^1.5

vz=vr*sin(theta)
vzo=vo*sin(theta1)

vy=vr*cos(theta)*sin(phi)
vyo=vo*cos(theta1)*sin(phi1)

vx=vr*cos(theta)*cos(phi)
vxo=vo*cos(theta1)*cos(phi1)
f=amp*exp(-( (vx-vxo)^2.0D +(vy-vyo)^2.0 + (vz-vzo)^2.0 )/(vth^2) )
fprime=(f)*cos(theta)*(vr)^3.0

return, fprime
end 
