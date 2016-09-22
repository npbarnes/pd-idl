; read 3d coordinate data 

PRO f_read_coord,file,x,y,z,dz_grid,dz_cell,nx,ny,nz
;file='coord.dat'

nx=0l
ny=0l
nz=0l

openr,11,file,/f77_unformatted
readu,11,nx
readu,11,ny
readu,11,nz
print,nx,ny,nz

x=fltarr(nx,/nozero)
y=fltarr(ny,/nozero)
z=fltarr(nz,/nozero)
dz_grid=fltarr(nz,/nozero)
dz_cell=fltarr(nz,/nozero)
readu,11,x
readu,11,y
readu,11,z
readu,11,dz_grid
readu,11,dz_cell

;q=fltarr(nx,ny,nz,3)

;for i=0,nx-1 do begin
;   for j=0,ny-1 do begin
;      for k=0,nz-1 do begin
;         for m=0,2 do begin
;	    q(i,j,k,0)=x(i)
;	    q(i,j,k,1)=y(j)
;            q(i,j,k,2)=z(k)
;         endfor
;      endfor
;   endfor
;endfor


;window,1,xsize=300,ysize=500,title='grid'
;plots,q(*,0,*,0)/max(x),q(*,0,*,2)/max(z),psym=3,/normal
;;plot,x(*,0,*),z(*,0,*),psym=3,/normal

;window,2,xsize=300,ysize=400,title='dz'
;dz=z(0,0,*)-shift(z(0,0,*),1)
;plot,dz,yrange=[0,max(dz<5000)]

close,11

return
end
