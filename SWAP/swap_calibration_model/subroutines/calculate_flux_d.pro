function calculate_flux_d ,vc,thetac,phic, w, s4,eff
aeff=0.033
aeff=aeff/0.0882d
dect_eff=eff
aeff=dect_eff*aeff
mp=1.67262158D-27
kb=1.380658D-23
e=1.6E-19 
ee=.5*mp/e*(vc*1000.)^2 
degrad=!dpi/180.D
szp=14.;15.
szt=20.
szv=5;10.;;;;12.
phi=findgen(szp)*(30)/(szp-1) -15.d
theta=findgen(szt)*10./(szt-1)  -5.d
phi=phi*degrad
theta=theta*degrad

f=dblarr(szp,szt, szv)
do_plot='n'
temp=dblarr(szp,szt)
;_______________
;___Read phi
;____________

wp=interpol(w.w, w.phi,phi/degrad)
;______________
;___Read Weenie
;________________

junk=min(abs(ee-s4.ecen),iee)
;print, ee, iee
temp1=dblarr(szp, szt)
temp2=dblarr(szp)
;*******Note: using -s4.x below to convert azimuth to theta angle)
for i=0, szp-1 do begin
  for j=0, szt-1 do begin
  
   er1=s4.y(iee,*)*ee
  
   ;print, format='(30F20.5)', er1
   junk=min(abs(theta(j)/degrad -(-s4.x)),ix)
   tt=s4.arr(*,ix,iee)
   ipos=where(tt gt 0)
   enew=min(er1(ipos)) + findgen(szv)*(max(er1(ipos))-min(er1(ipos)))/(szv-1.)
   
   ttnew=interpol(tt(ipos),er1(ipos), enew)
   vr1=sqrt(2.d*enew*e/mp)*100.d
   k=findgen(szv)
   ;print, format='(20F20.5)',ttnew, enew
   temp=FLUX(phi(i), theta(j), vr1,thetac,phic)*Aeff*wp(i)*ttnew 
   ;print, format='(20F20.5)', temp, enew
   f(i,j,k)=temp
   
    res1=int_tabulated(vr1,f(i,j,*),/double,/sort)  ; speed integral
    temp1(i,j)=res1
  endfor
  res2=int_tabulated(theta,temp1(i,*),/double,/sort) ;theta integral
  temp2(i)=res2
endfor
res=int_tabulated(phi,temp2,/double,/sort) ;phi integral

return, res
end
