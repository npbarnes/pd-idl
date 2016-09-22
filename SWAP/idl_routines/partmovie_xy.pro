PRO partmovie_xy,file,nfrm

;partmovie.pro
;reads sequential data file and generates movie

zm = 8  ;image enlargement factor

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

arr = fltarr(nx,ny,nfrm) 
icld = fltarr(nx,ny,nz,/nozero)
ncld = fltarr(nx,ny,nz,/nozero) 

for n = 0,nfrm-1 do begin
   readu,1,frm
   print,'image #.....',frm
   readu,1,icld
   print,'max np....',max(icld)
   print,'Ni tot....',total(icld)*0.5^3


;   readu,1,ncld
   for i=0,nx-1 do begin
      for j=0,ny-1 do begin
;         arr(i,j,n) = total(icld(i,j,*)+ncld(i,j,*))
         arr(i,j,n) = total(icld(i,j,*))
      endfor
   endfor
endfor

close,1


;movie,arr
sz=size(arr)
nfrms=sz(3)

arr2=congrid(bytscl(arr),nx*zm,ny*zm,nfrm)

xinteranimate,set=[nx*zm,ny*zm,nfrm]

for i=0,nfrm-1 do begin
   xinteranimate,image=arr2(*,*,i),frame=i
endfor

xinteranimate

end



