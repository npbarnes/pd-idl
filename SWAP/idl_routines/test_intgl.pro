pro test_intgl,m,x

f_read_2d_m,'t1_1',m,t1
f_read_2d_m,'t2_1',m,t3

sz = size(t1)
nx = sz(1)
nz = sz(2)

a1 = t1(x,*)
a3 = t3(x,*)

window,0
a1 = smooth(a1,2)
plot,a1,yrange=[-2*abs(max(a1)),+2*abs(max(a1))]
;oplot,a3,linestyle=1


sum = fltarr(nz)
for i = 1,nz-1 do begin
   sum(i) = total(a1(0:i))/2.0
   print,sum(i)
endfor

oplot,sum,linestyle=2

window,1
plot,a1-a3,linestyle=1,yrange=[-2*abs(max(a3)),+2*abs(max(a3))]

sum = fltarr(nz)
for i = 1,nz-1 do begin
   sum(i) = total(a1(0:i)-a3(0:i))/2.0
;   print,sum(i)
endfor

oplot,sum,linestyle=3
print,total(sum)
end