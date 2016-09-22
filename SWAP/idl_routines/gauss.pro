vo=1.33
vth=0.29

tau = 500.0
t=5.0

r=findgen(1000)/100
dr=r(1)-r(0)

ni = exp(-t/tau)*exp(-(r-vo*t)^2/(vth*t)^2)/(4*!pi*r^2*t*vth*sqrt(!pi))

niDH = exp(-t/tau)*exp(-(r-vo*t)^2/(vth*t)^2)/ $
           (4*!pi*vth*t^3*(sqrt(!pi)*(vth^2)/4.0 + vo*vth + sqrt(!pi)*(vo^2)/2.0))


plot,r,ni
oplot,r,niDH
print,total(4*!pi*r^2*ni*dr)
print,total(4*!pi*r^2*niDH*dr)


end