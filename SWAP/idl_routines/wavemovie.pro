PRO wavemovie,file,xyz,arr,temparr

;wavemovie.pro
;reads sequential data file and generates movie
;xyz specifies the component to view...0,1,2

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

print,nt,nout,nx,ny,nz

print,nt,nout,nt/nout
arr = fltarr(nx,nz,nt/nout) 
temparr = fltarr(nx,ny,nz,3,/nozero)

i=0
while (not(eof(1))) do begin
   readu,1,frm
   print,'image #.....',frm
   readu,1,temparr
   arr(*,*,i) = reform(temparr(*,ny/2,*,xyz))
;   surface,reform(temparr(*,ny/2,*,xyz))
;   print,median(temparr(*,*,*,xyz))
   wait,.2
;   print,median(temparr(*,*,*,xyz))
   i = i + 1
endwhile

;movie,arr
close,1

end




