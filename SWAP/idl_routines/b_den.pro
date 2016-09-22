;calculate density profile for simulation 
;!p.multi=[0,2,3]
 
re=6.38e6    ;earth radius in m 
ro=2.0e6     ;m 
M=8e22       ;Amps m^2
mu=!pi*4e-7  ;N/amp^2
mO=16*1.67e-27 ;kg

;window,1,xsize=450,ysize=680,title='win 1 
 
lambda=(!pi/4.0) - findgen(1000)*!pi/2000.0 
rearth = re*sqrt(cos(lambda)^2 + sin(lambda)^2) 
r=(re+ro)*cos(lambda)^2 
wh=where((r gt re+5e4) and (lambda ge 0)) 
!x.style=1 
!y.style=1 
;plot,rearth/1e3,lambda,/polar,xrange=[4000,9000],yrange=[-4000,4000],$
;     xtitle='km',ytitle='km 
;oplot,r/1e3,lambda,/polar,linestyle=1 

print,acos(sqrt(re/(re+ro)))*!radeg 
dl=lambda(0)-lambda(1) 
intgl = total(dl*(ro+re)*cos(lambda(wh))^2) 
print,intgl 
 
alt=[0,50,100,200,300,400,500,600,700,800,900,1000,1500,2000]*1e3  ;m 
nO=[0,2e3,7e4,1.3e5,5e5,1e6,5e5,2e5,6.9e4,3.2e4,2e4,1.5e4,1.0e4,0.9e4]*(1e2)^3 ;m 
t=(50+findgen(1000)*2)*1e3  ;m 
tt=r(wh)-re
tt= min(tt)+findgen(350)*(max(tt)-min(tt))/350.0                 ;m
nno=spline(alt,no,tt) 
plot_oo,nno,tt,xtitle='nO (m-3)',ytitle='alt (m) 
oplot,no,alt,psym=1 
 
nwh=n_elements(wh) 
n=findgen(nwh) 
a=findgen(nwh) 

for i=0,nwh-1 do begin 
   wh1=where(r(wh(i))-re gt tt) 
   nwh1=n_elements(wh1) 
   n(i)=nno(nwh1-1) 
   intgl = total(dl*r(wh(0):wh(i))*cos(lambda(wh(0):wh(i))^2)) 
   a(i)=intgl
endfor 

B = mu*M*sqrt(1+3*sin(lambda(wh))^2)/(4*!pi*r(wh)^3)
Va=B/sqrt(mu*mO*n)
;Vphi=sqrt(Va)

;plot_io,a/1e3,n,xtitle='Dist along B (km)',ytitle='nO (m-3)' 
;plot,a/1e3,b,ytitle='B (T)',xtitle='Dist along B (km)'
;plot,a/1e3,va/1e3,yrange=[0,1500],ytitle='Va (km/s)',xtitle='Dist along B (km)
;plot,a/1e3,vphi/1e3,ytitle='Vphi (km/s)',xtitle='Dist along B (km), $
;     yrange=[0,1.4]

openw,1,'oden.dat'
for i=0,n_elements(tt)-1 do printf,1,format='(e12.4,e12.4)',tt(i),nno(i)
close,1


end














