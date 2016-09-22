PRO waveplot,file,arr,nplts

;wavemovie.pro
;reads sequential data file and generates movie
;nplts is free parameter controling how many plots to stack.

close,1
openr,1,file,/f77_unformatted

nt=0
nout=0
nx=0
ny=0
nz=0
frm=0
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

arr = fltarr(nz,nt/nout) 
arr1 = fltarr(nz,nt/nout) 
temparr = fltarr(nx,ny,nz,3,/nozero)

;for i = 0,nt-1 do begin
i=0
while (not(eof(1))) do begin
   readu,1,frm
   print,'image #.....',frm
   readu,1,temparr
   arr(*,i) = reform(temparr(0,0,*,1))
   i=i+1
endwhile


for i=0,nz-1 do begin
   for j=0,(nt/nout)-1 do begin
      arr1(i,j) = arr(i,j) + 0.012*j
   endfor
endfor

plot,arr1(*,0),/noclip,xstyle=9,ystyle=8,xtitle='z',yrange=[0,2*max(arr)]

rng = fix(256./((nt/nout)/nplts))

for i=1,(nt/nout)/nplts do begin
   oplot,arr1(*,(i*nplts-1)),/noclip,color=255-i*rng
endfor

return
end







