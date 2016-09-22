;  calcuate dispersion relations for electromagnetic waves parallel to B

q = double(1.6e-19)    ;C
me = double(9.1e-31)   ;kg
ep = double(8.85e-12)  
q_over_me = q/me
c=double(3.0e8)        ;m/s

no = double(1e11)   ;1/m^3
B = double(3e-5)    ;T

wpe = sqrt(no*q*q/(ep*me))
wce = q_over_me*B

wr = 0.5*(wce + sqrt(wce^2 + 4*wpe^2))
wl = 0.5*(-wce + sqrt(wce^2 + 4*wpe^2))
wh = sqrt(wpe^2 + wce^2)
 
w=wpe + (wh-wpe)*(0.01+dindgen(101)/102)
fw=1.0/(1-(wpe^2/w^2)*((w^2 - wpe^2)/(w^2 - wh^2)))

plot,w,fw,xtitle='!7x',ytitle='!6(v!i!7u!n!6/c)!e2', $
          title='x mode (!7x!i!6p!n < !7x!6 < !7x!6!ih!n)
xyouts,0.6,0.7,'!7x!6!ip!n = ' + string(wpe),/normal
xyouts,0.6,0.65,'!7x!6!ic!n = ' + string(wce),/normal
xyouts,0.6,0.6,'!7x!6!ih!n = ' + string(wh),/normal
xyouts,0.6,0.75,'!6B = ' + string(B) + ' (t)',/normal
xyouts,0.6,0.8,'!6no = ' + string(no) + ' (m^-3)',/normal
print,wpe,wce,wr

end
