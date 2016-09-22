f_read_3d_vec_m,'Eall_1',10,b1
f_read_3d_m,'npall_1',10,np
f_read_coord,'coord.dat',xx,yy,zz,dzc,dzg,nx,ny,nz

set_plot,'ps
!p.font=0
device,/palatino
device,/encapsulated
device,filename='exy_cont.eps'


x1=4
x2 = 36
y1 = 3
y2 = 28

a = b1(x1:x2,y1:y2,*,*)
b = np(x1:x2,y1:y2,*)
x = xx(x1:x2)
y = yy(y1:y2)



zs = 91

velovect,smooth(a(*,*,zs,0),2),smooth(a(*,*,zs,1),2),x,y,$
   xrange=[fix(xx(x1)),round(xx(x2))],yrange=[round(yy(y1)),round(yy(y2))],$
   xtitle='x (km)',ytitle='y (km)',title = 'E (max = 1.0 V/m)'
   print,max(sqrt(a(*,*,zs,0)^2 + a(*,*,zs,1)^2))*1e3*16*1.27e-27/1.6e-19
  
contour,b(*,*,zs)/1e15,x,y,levels=[1e6,1e7,1e8,2e8,4e8,6e8],$
   c_thick=[2,1,1,1,1,1],$
   xstyle=1,ystyle=1,/noerase,$
   xrange=[fix(xx(x1)),round(xx(x2))],yrange=[round(yy(y1)),round(yy(y2))]

device,/close


end