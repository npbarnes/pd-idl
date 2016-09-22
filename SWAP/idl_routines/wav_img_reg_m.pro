PRO wav_img_reg_m,file,xyz,nfiles,nfrm,POSTSCRIPT=postscript

;xyz specifies which component to look at

if keyword_set(postscript) then begin
   set_plot,'ps
   device,bits=8
   device,filename='wave.ps'
   device,/inches,xsize=6.5,xoffset=1.0,ysize=9,yoffset=1.0
end

wavemovie_m,file,xyz,arr,tmparr,nfiles,nfrm

sz=size(arr)
nfrms=sz(3)

for i=0,nfrms-1 do begin
   a1=reform(arr(*,*,i)) 
   regrid,a1
   arr(*,*,i)=a1
   image_cont,a1
   wait,1
endfor

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x'
end

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