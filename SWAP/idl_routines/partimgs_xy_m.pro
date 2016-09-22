PRO partimgs_xy_m,nfrm,nfiles,POSTSCRIPT=postscript

;partmovie.pro
;reads sequential data file and generates movie

WINDOW, 0, XSIZE=500, YSIZE=740, TITLE='XY Images'

close,1
openr,1,'npall_1.dat',/f77_unformatted

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

close,1

reread = 1
if (reread eq 0) then begin


print,nt,nout,nx,ny,nz

arrxy = fltarr(nx,ny,nfrm*nfiles) 
arrxz = fltarr(nx,nz,nfrm*nfiles)
icld = fltarr(nx,ny,nz,/nozero)
ncld = fltarr(nx,ny,nz,/nozero) 

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
      for i=0,nx-1 do begin
         for j=0,ny-1 do begin
            arrxy(i,j,frmcount) = total(icld(i,j,*))
         endfor
      endfor
      for i=0,nx-1 do begin
         for k=0,nz-1 do begin
            arrxz(i,k,frmcount) = total(icld(i,*,k))
         endfor
      endfor
      frmcount=frmcount+1
   endfor

   close,1
endfor

;movie,arr
;sz=size(arrxy)
;nfrms=sz(3)

save,arrxy,filename='xyimage.sav'
save,arrxz,filename='xzimage.sav'
endif

restore,filename='/pf/d3/pad/hybrid/cont/xyimage.sav'
restore,filename='/pf/d3/pad/hybrid/cont/xzimage.sav'


zm = 3  ;image enlargement factor

maxcdenxy = max(arrxy)
maxcdenxz = max(arrxz)
if (maxcdenxy gt maxcdenxz) then maxcden = maxcdenxy
if (maxcdenxy le maxcdenxz) then maxcden = maxcdenxz

print,maxcden,maxcdenxy,maxcdenxz
arrxy2=congrid(bytscl(arrxy/maxcden),nx*zm,ny*zm,nfiles*nfrm)
arrxz2=congrid(bytscl(arrxz/maxcden),nx*zm,nz*zm,nfiles*nfrm)

;loadct,5
;stretch,10,200
psn1 = fltarr(nfiles)
psn2 = fltarr(nfiles)
psn1 = [0,2,4,6,8,10,12,14,16,18,20,22]
psn2 = psn1+1
for i=0,(nfiles*nfrm)-1 do begin
   tvscl,arrxy2(0:zm*75,zm*20:zm*40,i),psn1(i)
   print,'xy...',max(arrxy2(0:zm*75,zm*20:zm*40,i))
   tvscl,arrxz2(0:zm*75,zm*137:zm*157,i),psn2(i)
   print,'xz...',max(arrxz2(0:zm*75,zm*137:zm*157,i))
;   xinteranimate,image=arr2(*,*,i),frame=i
endfor

xloadct

print,max(arrxy)
print,max(arrxz)
print,total(arrxy(*,*,nfiles-1))
print,total(arrxz(*,*,nfiles-1))

cimg=tvrd()
tvscl,cimg
sz = size(cimg)

if keyword_set(postscript) then begin
set_plot,'ps
!p.font=0
device,filename='xyimages.eps'
device,bits=8
device,/palatino,/color
device,/portrait
device,/encapsulated
device,/inches,xsize=6.0,xoffset=1.0,ysize=6.0*sz(2)/sz(1),yoffset=1.0
end

tv,cimg*(!d.n_colors-35)/max(cimg)


if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
endif

WINDOW, 0, XSIZE=25, YSIZE=200, TITLE='Color Bar'

cbar = bytarr(25,!d.n_colors)

for j = 0,!d.n_colors-1 do begin
   cbar(*,j) =  j
endfor

tvscl,bytscl(cbar)
sz = size(cbar)

if keyword_set(postscript) then begin
set_plot,'ps
!p.font=0
device,filename='cbar.eps'
device,bits=8
device,/palatino,/color
device,/portrait
device,/encapsulated
device,/inches,xsize=0.5,xoffset=1.0,ysize=0.5*sz(2)/sz(1),yoffset=1.0
end

cbar = bytarr(25,!d.n_colors)
for j = 0,!d.n_colors-1 do begin
   cbar(*,j) =  j
endfor
tv,cbar*(!d.n_colors-2)/max(cbar)

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
endif


end



