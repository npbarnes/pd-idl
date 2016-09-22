PRO wavemovie_m,file,xyz,arr,temparr,nfiles,nfrm

;wavemovie.pro
;reads sequential data file and generates movie
;xyz specifies the component to view...0,1,2

close,1

nt=0
nout=0
nx=0
ny=0
nz=0
frm=0

openr,1,file+'1.dat',/f77_unformatted

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

print,nt,nout,nx,ny,nz

arr = fltarr(nx,nz,nfiles*nfrm) 
temparr = fltarr(nx,ny,nz,3,/nozero)

close,1

frmcount=0
for m=1,nfiles do begin

   files = file+strtrim(string(m),1)+'.dat'
   openr,1,files,/f77_unformatted
   print,'Reading file.....',files

   readu,1,nt
   readu,1,nout
   readu,1,nx
   readu,1,ny
   readu,1,nz

   while (not(eof(1))) do begin
      readu,1,frm
      print,'  image #.....',frm
      readu,1,temparr
      arr(*,*,frmcount) = reform(temparr(*,ny/2,*,xyz))
;      surface,reform(temparr(*,ny/2,*,xyz))
;      print,median(temparr(*,*,*,xyz))
;      wait,.2
;      print,median(temparr(*,*,*,xyz))
      frmcount = frmcount + 1
   endwhile

close,1

endfor

end




