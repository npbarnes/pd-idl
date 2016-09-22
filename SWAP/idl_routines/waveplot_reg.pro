PRO waveplot_reg,file,xyz,arr,nplts

;wavemovie.pro
;reads sequential data file and generates movie
;nplts is free parameter controling how many plots to stack.
;xyz specifies plane to slice through

close,1
openr,1,file,/f77_unformatted
nt=0
nout=0
readu,1,nt
readu,1,nout
close,1

f_read_coord,'coord.dat',x,y,z,dz_grid,dz_cell,nx,ny,nz

wavemovie,file,xyz,arr

sz = size(arr)
nfrms = sz(3)

arr1 = fltarr(nz,nt/nout) 

for i=0,nfrms-1 do begin
   a1=reform(arr(*,*,i))
   regrid,a1,xyz,sz,x,y,z,dz,nx,ny,nz
   arr1(*,i) = reform(a1(0,*))
endfor

for i=0,nz-1 do begin
   for j=0,(nt/nout)-1 do begin
      arr1(i,j) = arr1(i,j) + 0.012*j
   endfor
endfor

plot,arr1(*,0),/noclip,xstyle=9,ystyle=8,xtitle='z',yrange=[0,2*max(arr)]

rng = fix(256./((nt/nout)/nplts))

for i=1,(nt/nout)/nplts do begin
   oplot,arr1(*,(i*nplts-1)),/noclip,color=255-i*rng
endfor

return
end







