;------------------------------------------------------------
pro cross,a,b,c
;------------------------------------------------------------

   c(0) = a(1)*b(2) - a(2)*b(1)
   c(1) = a(2)*b(0) - a(0)*b(2)
   c(2) = a(0)*b(1) - a(1)*b(0)

return
end
;------------------------------------------------------------


;-------------------------------------------------------------
pro dot,a,b,c
   c = a(0)*b(0) + a(1)*b(1) + a(2)*b(2)
return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
PRO Initialize,f,p,g,u,e
;-------------------------------------------------------------

e.initial_E = 0.0
vol = g.dx(0)*g.dx(1)*g.dx(2)
f.B(*,*,*,2) = 300.0
for i=0,g.ni-1 do begin
   p(i).v(0) = 10.0 + 2.0*(0.5-randomu(seed))
   p(i).v(1) = 4.0*(0.5-randomu(seed))
   p(i).v(2) = 4.0*(0.5-randomu(seed))
   p(i).x(0) = g.ri*g.dx(0)+g.dx(0)*(0.5 - randomu(seed))
   p(i).x(1) = g.rj*g.dx(1)+g.dx(1)*(0.5 - randomu(seed))
   p(i).x(2) = g.rk*g.dx(2)+g.dx(2)*(0.5 - randomu(seed))
   p(i).ijkp(0) = fix(p(i).x(0)/g.dx(0))
   p(i).ijkp(1) = fix(p(i).x(1)/g.dx(1))   
   p(i).ijkp(2) = fix(P(i).x(2)/g.dx(2))
   ii = p(i).ijkp(0)
   jj = p(i).ijkp(1)
   kk = p(i).ijkp(2)
   u.np(ii,jj,kk) = u.np(ii,jj,kk) + 1.0/(g.dx(0)*g.dx(1)*g.dx(2))
   p(i).v1 = p(i).v
endfor

;for i = 0,g.ni-1 do begin
;   ii = p(i).ijkp(0)
;   jj = p(i).ijkp(1)
;   kk = p(i).ijkp(2)
;   u.up1(ii,jj,kk,*) = u.up1(ii,jj,kk,*) + p(i).v(*)/(u.np(ii,jj,kk)*vol)
;;   print,i,p(i).v(0),u.up1(ii,jj,kk,0)
;;   print,i,p(i).v(1),u.up1(ii,jj,kk,1)
;;   print,i,p(i).v(2),u.up1(ii,jj,kk,2)
;endfor

u.up = u.up1
e.initial_E = total(p.v(0)^2 + p.v(1)^2 + p.v(2)^2)

return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
pro get_energy,p,e
;-------------------------------------------------------------

e.ve = total(p.v(0)^2 + p.v(1)^2 + p.v(2)^2)
print,'Normalized enery....',e.ve/e.initial_E

return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
pro extrapol_up,u,p,g
;-------------------------------------------------------------

;u.up = 1.5*u.up1 - 0.5*u.up
;u.up3 = u.up1
v_at_n = fltarr(3)  ;extrapolated value for v at time level n

vol = g.dx(0)*g.dx(1)*g.dx(2)
u.up=0.0
for i=0,g.ni-1 do begin
   ii=p(i).ijkp(0)
   jj=p(i).ijkp(1)
   kk=p(i).ijkp(2)
   v_at_n = 1.5*p(i).v - 0.5*p(i).v1
   u.up(ii,jj,kk,*) = u.up(ii,jj,kk,*) + v_at_n/(u.np(ii,jj,kk)*vol)   

endfor


return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
pro get_Ep,f,u
;-------------------------------------------------------------


f.Ep(*,*,*,0) = u.up(*,*,*,1)*f.B(*,*,*,2) - u.up(*,*,*,2)*f.B(*,*,*,1)
f.Ep(*,*,*,1) = u.up(*,*,*,2)*f.B(*,*,*,0) - u.up(*,*,*,0)*f.B(*,*,*,2)
f.Ep(*,*,*,2) = u.up(*,*,*,0)*f.B(*,*,*,1) - u.up(*,*,*,1)*f.B(*,*,*,0)

f.Ep = -f.Ep
;print,f.Ep(where

return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
pro improve_up,p,f,u,g
;-------------------------------------------------------------

u.up = 0.0

B = reform(f.B(0,0,0,*))
B2=0.0
dot,B,B,B2

a_d = 1 + (B2*g.dt*g.dt/4.0)
a1 = (1-(B2*g.dt*g.dt/4.0))/a_d
a2 = g.dt/a_d
a3 = 0.5*g.dt*g.dt/a_d
print,'coeff...',a_d,a1,a2,a3

vm=fltarr(3)
vp=fltarr(3)
vmxB=fltarr(3)
vmdotB=0.0

vol = g.dx(0)*g.dx(1)*g.dx(2)
u.up = 0.0
for i=0,g.ni-1 do begin

   ii = p(i).ijkp(0)
   jj = p(i).ijkp(1)
   kk = p(i).ijkp(2)

   vm = p(i).v + 0.5*g.dt*f.Ep(ii,jj,kk,*)
   cross,vm, B, vmxB
   dot,vm, B, vmdotB
   vp = a1*vm + a2*vmxB + a3*vmdotB*B
   v_at_n = 0.5*(vm + vp)
   p(i).v1 = v_at_n
   u.up(ii,jj,kk,*) = u.up(ii,jj,kk,*) + v_at_n(*)/(u.np(ii,jj,kk)*vol)
 
endfor

return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
pro move_ions,g,p,f,u
;-------------------------------------------------------------
xplus = fltarr(g.ni,3)             ;half time step position
nplus = fltarr(g.nx,g.ny,g.nz)     ;half time step density
ijk = fltarr(g.ni,3)               ;half time step indices

B = reform(f.B(0,0,0,*))
B2=0.0
dot,B,B,B2

a_d = 1 + (B2*g.dt*g.dt/4.0)
a1 = (1-(B2*g.dt*g.dt/4.0))/a_d
a2 = g.dt/a_d
a3 = 0.5*g.dt*g.dt/a_d

vm=fltarr(3)
vp=fltarr(3)
vmxB=fltarr(3)
vmdotB=0.0

vol = g.dx(0)*g.dx(1)*g.dx(2)

for i=0,g.ni-1 do begin

;   p(i).v1 = p(i).v

   ii = p(i).ijkp(0)
   jj = p(i).ijkp(1)
   kk = p(i).ijkp(2)

   vm = p(i).v + 0.5*g.dt*f.Ep(ii,jj,kk,*)
   cross,vm, B, vmxB
   dot,vm, B, vmdotB
   vp = a1*vm + a2*vmxB + a3*vmdotB*B

   p(i).v = vp + 0.5*g.dt*f.Ep(ii,jj,kk,*)    ;half time step move
;   xplus(i,*) = p(i).x + 0.5*g.dt*p(i).v
;   ijk(i,*) = fix(xplus(i,*)/g.dx(*))
;   nplus(ijk(i,0),ijk(i,1),ijk(i,2)) = nplus(ijk(i,0),ijk(i,1),ijk(i,2)) + $
;                                       1.0/vol   

   p(i).x = p(i).x + g.dt*p(i).v            ;full time step move
   p(i).ijkp = fix(p(i).x/g.dx) 
endfor

;u.up1 = 0.0
;for i=0,g.ni-1 do begin       ;update up1 at half time step postion
;   ii = ijk(i,0)
;   jj = ijk(i,1)
;   kk = ijk(i,2)
;   u.up1(ii,jj,kk,*) = u.up1(ii,jj,kk,*) + p(i).v/(nplus(ii,jj,kk)*vol)
;endfor

u.np=0                              ;update densities at full time step
for i=0,g.ni-1 do begin
   ii = p(i).ijkp(0)
   jj = p(i).ijkp(1)
   kk = p(i).ijkp(2)
   u.np(ii,jj,kk) = u.np(ii,jj,kk) + 1.0/vol 
endfor

return
end
;-------------------------------------------------------------



;------------------------------------------------------------- 
PRO Move_neutral,p,u,g
;-------------------------------------------------------------

u.np = 0.0
for i=0,g.ni-1 do begin
   for m=0,2 do begin   
      p(i).x(m) = p(i).x(m) + g.dt*p(i).v(m)
      p(i).ijkp(m) = p(i).x(m)/g.dx(m)
   endfor
   ii = p(i).ijkp(0)
   jj = p(i).ijkp(1)
   kk = p(i).ijkp(2)
   u.np(ii,jj,kk) = 1.0/(g.dx(0)*g.dx(1)*g.dx(2))
endfor

return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
pro make_plot,p,g
;-------------------------------------------------------------

xmax = max(p.x(0))
zmax = max(p.x(2))
;erase
loadct,2
;surface,dist(5),xrange=[0,g.nx*g.dx(0)],$
;                yrange=[0,g.ny*g.dx(1)],$
;                zrange=[0,g.nz*g.dx(2)],/save,/nodata

plot,[0,fix(g.nx*g.dx(0))],[0,fix(g.nz*g.dx(2))],/nodata
plots,p.x(0),p.x(2),psym=3,color=30
;for i=0,g.ni-1 do begin
;   plots,p(i).x(0),p(i).x(1),p(i).x(2),psym=3,color=i*20,/T3D  
;endfor
plots,g.ri*g.dx(0), g.rk*g.dx(2), psym=2,color=80
return
end
;-------------------------------------------------------------


;-------------------------------------------------------------
; main
;-------------------------------------------------------------
ni = 500
nx = 20
ny = 6
nz = 6
dt = 0.05
nt = 1
ri=2
rj=ny/2
rk=nz/2

g = {info, dx: fltarr(3), nx: 0, ny: 0, nz: 0, ni: 0, dt: 0.0, $
           ri: 0, rj: 0, rk: 0}
p = replicate({particles, v: fltarr(3), x: fltarr(3), $
                          v1: fltarr(3), ijkp: fltarr(3)},ni)
f = {fields, B: fltarr(nx,ny,nz,3), Ep: fltarr(nx,ny,nz,3)}
u = {fluid, up: fltarr(nx,ny,nz,3), up1: fltarr(nx,ny,nz,3), $
            up3: fltarr(nx,ny,nz,3), np: fltarr(nx,ny,nz)}
e = {energy, initial_E: 0.0, ve: 0.0}


g.dx(0) = .5
g.dx(1) = .5
g.dx(2) = .5
g.nx = nx
g.ny = ny
g.nz = nz
g.ni = ni
g.dt = dt

g.ri=ri
g.rj=rj
g.rk=rk


Initialize,f,p,g,u,e

for n = 0,nt-1 do begin
   time = g.dt*n
   print,'time.....',n,time
;   Move_neutral,p,u,g
   extrapol_up,u,p,g
   get_Ep,f,u
   improve_up,p,f,u,g
   get_Ep,f,u
   move_ions,g,p,f,u

   get_energy,p,e
   make_plot,p,g
endfor



end
;-------------------------------------------------------------



