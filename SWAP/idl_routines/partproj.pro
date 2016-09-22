PRO partproj,file,rk,ry,rx

;partproj.pro
;reads sequential data file and generates movie
;xyz specifies the component to view...0,1,2
;projection plane will be the x-z plane...choose appropriate projections.

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

print,nt,nout,nt/nout
arr = fltarr(nx,nz,nt/nout) 
icld = fltarr(nx,ny,nz,/nozero)
ricld = fltarr(nx,ny,nz,/nozero)
;ncld = fltarr(nx,ny,nz,/nozero) 

nfrms=0
while (not(eof(1))) do begin
   readu,1,frm
   print,'image #.....',frm
   readu,1,icld
;   readu,1,ncld

   for i=0,nx-1 do begin                     ;rotation about x
      ricld(i,*,*) = rot(reform(icld(i,*,*)),0)
   endfor

   for k=0,nz-1 do begin                     ;rotation about z
      ricld(*,*,k) = rot(reform(ricld(*,*,k)),rk)
   endfor

   for j=0,ny-1 do begin                     ;rotation about y
      ricld(*,j,*) = rot(reform(ricld(*,j,*)),ry)
   endfor

   for i=0,nx-1 do begin                     ;rotation about x
      ricld(i,*,*) = rot(reform(ricld(i,*,*)),rx)
   endfor

   for i=0,nx-1 do begin
      for k=0,nz-1 do begin
;         arr(i,k,nfrms) = total(icld(i,*,k)+ncld(i,*,k))
         arr(i,k,nfrms) = total(ricld(i,*,k))
      endfor
   endfor
   nfrms = nfrms + 1
endwhile

close,1

;movie,arr
sz=size(arr)
nfrms=sz(3)

arr2=congrid(bytscl(arr),nx*1,nz*1,nfrms)

xinteranimate,set=[nx*1,nz*1,nfrms]

for i=0,nfrms-1 do begin
   xinteranimate,image=arr2(*,*,i),frame=i
endfor

xinteranimate

xinteranimate,/close

end


