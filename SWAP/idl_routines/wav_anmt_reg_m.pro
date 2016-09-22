PRO wav_anmt_reg_m,file,xyz,nfiles,nfrm

;xyz specifies which axis to slice across...i.e. 1 means slice across y
;to get the xz plane.

f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz

wavemovie_m,file,xyz,arr,tmparr,nfiles,nfrm

sz=size(arr)
nfrms=sz(3)

for i=0,nfrms-1 do begin
   a1=reform(arr(*,*,i)) 
;   regrid,a1,xyz,sz,x,y,z,dzg,nx,ny,nz
   regrid,a1
   arr(*,*,i)=a1
endfor

arr1=bytscl(arr)
arr2=congrid(arr1,100,200,nfrms)

xinteranimate,set=[100,200,nfrms]

for i=0,nfrms-1 do begin
   xinteranimate,image=arr2(*,*,i),frame=i
endfor

xinteranimate

xinteranimate,/close

return
end