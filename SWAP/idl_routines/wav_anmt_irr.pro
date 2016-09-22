PRO wav_anmt_irr,file,xyz

;xyz specifies which axis to slice across...i.e. 1 means slice across y
;to get the xz plane.

wavemovie,file,xyz,arr

sz=size(arr)
nfrms=sz(3)

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