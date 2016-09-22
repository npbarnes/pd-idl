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

w=wce*(0.01+dindgen(101)/102)
fw=1/(1-((wpe^2/w^2)/(1-wce/w)))

;plot,w,fw,xtitle='!7x',ytitle='!6(v!i!7u!n!6/c)!e2',title='whistler mode
;xyouts,0.4,0.5,'!6v!7!ix!6!n =' + string(sqrt(max(fw))*c/1e3)+' km/s',/normal
;xyouts,0.4,0.45,'B = '+string(B)+' (T)',/normal
;xyouts,0.4,0.4,'no = '+string(no)+' (m^-3)',/normal

vphi = sqrt(fw)*c
k=w/vphi
plot,w,k,ytitle='k (m!e-1!n)',xtitle='!7x!6 (rad/sec)',title='whistler mode
xyouts,0.4,0.65,'B = '+string(B)+' (T)',/normal
xyouts,0.4,0.7,'no = '+string(no)+' (m^-3)',/normal

print,wpe,wce,wr

end
