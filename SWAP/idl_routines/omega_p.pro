;  calcuate dispersion relations for electromagnetic waves parallel to B

q = double(1.6e-19)    ;C
me = double(9.1e-31)   ;kg
ep = double(8.85e-12)  
q_over_me = q/me
c=double(3.0e8)        ;m/s


no = 1e10 + (1e12 - 1e10)*dindgen(20) ;1/m^3
B = 1e-5 + (5e-5 - 1e-5)*dindgen(20)     ;T

wpe = sqrt(no*q*q/(ep*me))
wce = q_over_me*B

plot,no,wpe
plot,b,wce


end
