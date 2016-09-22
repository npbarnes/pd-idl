PRO partmovie_yz,file,nfrm

;partmovie.pro
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

;print,nt,nout,nfrm
arr = fltarr(nx,nz,nfrm) 
icld = fltarr(nx,ny,nz,/nozero)
ncld = fltarr(nx,ny,nz,/nozero) 

for n=0,nfrm-1 do begin
   readu,1,frm
   print,'image #.....',frm
   readu,1,icld
;   readu,1,ncld
   for i=0,nx-1 do begin
      for k=0,nz-1 do begin
;         arr(i,k,n) = total(icld(i,*,k)+ncld(i,*,k))
         arr(i,k,n) = total(icld(i,*,k))
      endfor
   endfor
endfor

close,1

arr2=congrid(bytscl(arr),nx*8,nz*8,nfrm)

xinteranimate,set=[nx*8,nz*8,nfrm]

for i=0,nfrm-1 do begin
   xinteranimate,image=arr2(*,*,i),frame=i
endfor

xinteranimate

end


