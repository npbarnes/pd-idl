;-------------------------------------------------------------------
pro btrace,file,b
;btrace.pro  Traces magnetic field lines.
;-------------------------------------------------------------------

f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz
f_read_3d_vec,file,b
print,nx,ny,nz

dx = max(x)/nx
dz = max(z)/nz
print,dx,dz

bx = reform(b(*,18,*,0))
bz = reform(b(*,18,*,2))
regrid,bx
regrid,bz

;surface,bx
;surface,bz
;wait,10

contour,bx,/nodata,xrange=[0,nx*dx],yrange=[0,nz*dz],xstyle=1,ystyle=1
;plots,5,100,/data,psym=1


x0 = 0.0

repeat begin

x0 = x0 + 0.5
x1 = x0
z1 = 0.0
flg = 0

while (flg ne 1) do begin
  i = round(nx*x1/max(x))
  k = round(nz*z1/max(z))

  if (i gt nx-1) then goto, BAIL
  if (k gt nz-1) then goto, BAIL


  bx1 = bx(i,k)
  bz1 = bz(i,k)
;  print,'bx1,bz1...',bx1,bz1,i,k
  z2 = z1+0.5*dz
;  x2 = x1+dx
  x2 = x1+dz*bx1/bz1
;  print,'x2,z2...',x2,z2	
  plots,[x1,x2],[z1,z2],/data
  x1 = x2
  z1 = z2

  if (x1 lt 0) then flg = 1
  if (x1 gt max(x)) then flg = 1

  if (z1 lt 0) then flg = 1
  if (z1 gt max(z)) then flg = 1

endwhile

BAIL:

endrep until(x1 gt max(x))


;contour,bx

return
end
;-------------------------------------------------------------------