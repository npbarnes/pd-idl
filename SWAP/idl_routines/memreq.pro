pro memreq,nx,ny,nz,ni
;nx,ny,nz are the x,y,z dimensions of the arrays
;ni is the maximum number of particles

f3dvector = 13.*nx*ny*nz*3.
pvector = 8.*ni*3.
pscalar = 3*ni
f3dscalar = 3.*nx*ny*nz
f1dscalar = 1.*nz
p8d = 1.*ni*8.

print,f3dvector,pvector,f3dscalar,f1dscalar,pscalar,p8d


tot = (f3dvector+pvector+f3dscalar+f1dscalar+pscalar+p8d)

print,'f3dvector....',f3dvector/tot
print,'pvector......',pvector/tot
print,'f3dscalar....',f3dscalar/tot
print,'f1dscalar....',f1dscalar/tot
print,'pscalar......',pscalar/tot
print,'p8d..........',p8d/tot


print,''
print,'Total Memory Requirements:'
print,''
print,'     64 bit...',tot*8/1e6,'(MB)',tot,'(Mw)
print,'     32 bit...',tot*4/1e6,'(MB)'



return 
end
