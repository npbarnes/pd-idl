PRO partmovie_xy_m,nfrm,nfiles

;partmovie.pro
;reads sequential data file and generates movie

zm = 6  ;image enlargement factor

close,1
openr,1,'npall_1.dat',/f77_unformatted

nt=0ll
nout=0ll
nx=0ll
ny=0ll
nz=0ll
frm=0ll
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

close,1

print,nt,nout,nx,ny,nz

arr = dblarr(nx,ny,nfrm*nfiles) 
icld = dblarr(nx,ny,nz,/nozero)
ncld = dblarr(nx,ny,nz,/nozero) 

frmcount=0
for m=1,nfiles do begin

   file = 'npall_'+strtrim(string(m),1)+'.dat'
   print,'Reading....',file

   openr,1,file,/f77_unformatted

   readu,1,nt
   readu,1,nout
   readu,1,nx
   readu,1,ny
   readu,1,nz

   for n = 1,nfrm do begin

      readu,1,frm
      print,'   image #.....',frm
      readu,1,icld
;      print,'max,min icld...',max(icld),min(icld)
;      readu,1,ncld

      for i=0,nx-1 do begin
         for j=0,ny-1 do begin
;            arr(i,j,frmcount) = total(icld(i,j,*)+ncld(i,j,*))
            arr(i,j,frmcount) = total(icld(i,j,*))
         endfor
      endfor
      frmcount=frmcount+1
   endfor

   close,1
endfor

;movie,arr
sz=size(arr)
nfrms=sz(3)

zm = 512.0/nx
arr2=congrid(bytscl(arr),nx*zm,ny*zm,nfiles*nfrm)

;for i = 0,nfiles*nfrm-1 do begin
;   file = '/usr/users/pad/apex/model_'+strtrim(string(i),2)+'.gif'
;   write_gif,file,arr2(*,*,i)
;endfor

xinteranimate,set=[nx*zm,ny*zm,nfiles*nfrm]

for i=0,(nfiles*nfrm)-1 do begin
   xinteranimate,image=arr2(*,*,i),frame=i
endfor

xinteranimate

end



