;  calcuate dispersion relations for electromagnetic waves parallel to B

q = double(1.6e-19)    ;C
me = double(9.1e-31)   ;kg
ep = double(8.85e-12)  
mu = double(!pi*4e-7)
mO = 16*1.67e-27
q_over_me = q/me
c=double(3.0e8)        ;m/s

no = double(1e11)   ;1/m^3
B = double(3e-5)    ;T

va=B/sqrt(mu*mO*no)
print,va/1e3


end
