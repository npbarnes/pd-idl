; read 2d coordinate data 

PRO f_read_coord_2d_32,file,x,z,dz_grid,dz_cell,nx,nz

nx=0l
nz=0l

openr,11,file,/f77_unformatted
readu,11,nx
readu,11,nz
print,nx,nz

x=fltarr(nx,/nozero)
z=fltarr(nz,/nozero)
dz_grid=fltarr(nz,/nozero)
dz_cell=fltarr(nz,/nozero)
readu,11,x
readu,11,z
readu,11,dz_grid
readu,11,dz_cell

close,11

return
end
