set_plot,'ps
!p.font=0
!x.thick=2
!y.thick=2
!p.thick=2
device,filename='dptmp.ps



nT=100
nrh=120

tempF = 50.+ findgen(nT)/2
rh = 20 + findgen(nrh)/2
tempC = (tempF-32)*5./9.
svp = 6.11*10^(7.5*tempC/(237.7+tempC))

vp = fltarr(nT,nrh)
dptC = fltarr(nT,nrh)
for i = 0,nT-1 do begin
   for j = 0,nrh-1 do begin
     vp(i,j) = svp(i)*rh(j)/100.
     dptC(i,j) =(-430.22+237.7*alog(vp(i,j)))/(-alog(vp(i,j))+19.08)

   endfor 
endfor


dptF = (9./5.)*dptC + 32.

contour,dptF,tempF,rh,nlev=10,c_label=[1,1,1,1,1,1,1,1,1,1],yrange=[20,80],$
  levels=[10,20,30,40,50,60,70,80,90],ytitle='Relative Humidity (%)',$
  xtitle='Temp (F)',c_thick=[1,1,1,1,1,1,3,1]

plots,[85,85],[!y.crange(0),!y.crange(1)],/data,linestyle=1

device,/close


end
