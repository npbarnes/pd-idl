n = 10.0e20
dxr = 0.5
alpha=2.0e-21
b0r = 3.0e-5
q=1.6e-19
mO=16e-27

b0r = q*b0r/mO
kr = 2/dxr

a1 = kr^2*b0r/(alpha*n)
a2 = kr^2*b0r^2/(alpha*n)

omegar = 0.5*(a1 + sqrt(a1^2 + 4*a2))
vphi = omegar/kr
print,b0r,omegar,vphi

qr = 400.0
;--------------------------------------------
z=findgen(201)

dz = 0.5 + 0.001*(z-100)^2
plot,z,dz,psym=1

qz = fltarr(201)

qz(100) = qr
for i=101,200 do qz(i) = qz(i-1) + dz(i)
for i=0,99 do qz(99-i) = qz(100-i) - dz(99-i)

k=2/dz

b0 = 0.5*(-omegar + (omegar/k)*sqrt(k^2 + 4*alpha*n))

plot,b0

end