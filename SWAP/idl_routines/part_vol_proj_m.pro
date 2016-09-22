PRO part_vol_proj_m,nfrm,nfiles,POSTSCRIPT=postscript

;partmovie.pro
;reads sequential data file and generates movie

zm = 8.0  ;image enlargement factor
;!p.multi=[0,4,4]
if keyword_set(postscript) then begin
   set_plot,'ps
   device,bits=4
   device,filename='part_img.ps
   device,/inches,xsize=6.5,xoffset=1.0,ysize=9.0,yoffset=1.0
end

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

print,nt,nout,nx,ny,nz

arr = fltarr(nx,ny,nfrm*nfiles) 
icld = fltarr(nx,ny,nz,/nozero)
ncld = fltarr(nx,ny,nz,/nozero) 

frmcount=0

T3D, /RESET
T3D, TRANSLATE=[-0.5, -0.5, -0.5]
T3D, SCALE=[0.7, 0.7, 0.7]
T3D, ROTATE=[30, -30, 60]
T3D, TRANSLATE=[0.5, 0.5, 0.5]




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
;      readu,1,ncld
;      for i=0,nx-1 do begin
;         for j=0,ny-1 do begin
;;            arr(i,j,frmcount) = total(icld(i,j,*)+ncld(i,j,*))
;            arr(i,j,frmcount) = total(icld(i,j,*))
;         endfor
;      endfor
      arr1 = PROJECT_VOL(icld, 64,64,64,DEPTH_Q=0.7, $
                OPAQUE=opaque, TRANS=(!P.T))
;      tvscl,smooth(congrid(arr(*,*,frmcount),nx*zm,ny*zm),4),frmcount
;;         title='frame #....'+strtrim(string(frmcount),1)
      tvscl,congrid(arr1,200,200),frmcount
      frmcount=frmcount+1
   endfor

   close,1
endfor

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
end


end



