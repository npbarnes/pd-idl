file = 'f77test1.dat'
close,1
openr,1,file,/f77_unformatted

nx=0.
ny=0.
nz=0.

readu,1,nx
readu,1,ny
readu,1,nz

print,nx,ny,nz

end



