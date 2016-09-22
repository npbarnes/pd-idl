pro integrate_t1,j

f_read_3d_m,'t1_1',1,t1

;t1 = smooth(t1,3)

a = t1(1,j,*)


sz  = size(t1)

b = fltarr(sz(3))

b(0) = a(0)

for i = 1,sz(3)-1 do begin
   b(i) = total(a(0:i))
endfor

!p.multi=[0,1,2]

window,0
plot,a
plot,b

end