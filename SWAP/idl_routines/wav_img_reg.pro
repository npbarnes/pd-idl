PRO wav_img_reg,file,xyz,POSTSCRIPT=postscript

;xyz specifies which component to look at

f_read_coord,'coord.dat',x,y,z,dz_grid,dz_cell,nx,ny,nz
xarr = min(x) + findgen(n_elements(x))*max(x)/(n_elements(x)-1)
yarr = min(y) + findgen(n_elements(y))*max(y)/(n_elements(y)-1)
zarr = min(z) + findgen(n_elements(z))*max(z)/(n_elements(z)-1)

close,1
openr,1,file,/f77_unformatted
nt=0
nout=0
nx=0
ny=0
nz=0
dt=0.004
tbegin=0.2
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

temparr = fltarr(nx,ny,nz,3,/nozero)
frm=0
fr=0
while (not(eof(1))) do begin
   readu,1,fr
   frm=[frm,fr]
   readu,1,temparr
endwhile
frm = frm(1:*)

tit = strmid(strtrim(string(frm*dt + tbegin),1),0,4) + ' s'

if keyword_set(postscript) then begin
   set_plot,'ps
   loadct,0
;   stretch,0,240
   !p.font=7
   !p.multi=[0,3,4]
    device,bits=4
   device,/color,/times,font_size=12
   device,filename='wave.ps'
   device,/inches,xsize=6.5,xoffset=1.0,ysize=9,yoffset=1.0
end

wavemovie,file,xyz,arr,tmparr

sz=size(arr)
nfrms=sz(3)

maxarr = fltarr(nfrms)
minarr = fltarr(nfrms)
for i=0,nfrms-1 do begin
   maxarr(i) = max(arr(*,*,i))
   minarr(i) = min(arr(*,*,i))
endfor

arr=bytscl(arr)
for i=0,nfrms-1 do begin
   a1=reform(arr(*,*,i)) 
   regrid,a1
   arr(*,*,i)=a1
   for j=0,sz(2)-1 do begin
      arr(sz(1)-fix(sz(1)*0.1):*,j,i) = j*!d.n_colors/(sz(2)-1)
   endfor
   print,'tit...',tit(i)
   image_cont_xy,arr(*,*,i),xarr,zarr,tit(i) ;,/window_scale,/aspect
   y1vals = [strmid(strtrim(string(minarr(i)),1),0,4), $
             strmid(strtrim(string(maxarr(i)),1),0,4)]
   axis, yaxis=1, yrange = [minarr(i),maxarr(i)],yticks=1,ystyle=1, $
         ytickv = y1vals, charsize = 1.1
   wait,0.1
endfor

if keyword_set(postscript) then begin
   xyouts,0.6*!d.x_size,1.05*!d.y_size,systime(0),color=0,/device,charsize=1.1
   xyouts,0.6*!d.x_size,1.02*!d.y_size,file + string(xyz),color=0,/device,charsize=1.1

   device,/close
   set_plot,'x'
   !p.multi=[0,1,1]
   !p.font=-1
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

;--------------------------------------------------------------------

;---------------------------------------------------------------------;
;  NAME:  Peter                                                       ;
;  DATE:  7/20/93                                                     ;
;  FILE:  image_cont_xy.pro                                           ;
;  SUMY:  image cont EXPAND.C simulation results...                   ;
;---------------------------------------------------------------------;

pro image_cont_xy, a, ax, ay, tit,  WINDOW_SCALE = window_scale, $
                ASPECT = aspect, INTERP = interp
;+
; NAME:
;	IMAGE_CONT
; PURPOSE:
;	Overlay an image and a contour plot.
; CATEGORY:
;	General graphics.
; CALLING SEQUENCE:
;	IMAGE_CONT, A
; INPUTS:
;	A = 2 dimensional array to display.
; KEYWORD PARAMETERS:
;	/WINDOW_SCALE = set to scale the window size to the image size,
;		otherwise the image size is scaled to the window size.
;		Ignored when outputting to devices with scalable pixels.
;	/ASPECT = set to retain image's aspect ratio.  Assumes square
;		pixels.  If /WINDOW_SCALE is set, the aspect ratio is
;		retained.
;	/INTERP = set to bi-linear interpolate if image is resampled.
; OUTPUTS:
;	No explicit outputs.
; COMMON BLOCKS:
;	none.
; SIDE EFFECTS:
;	The currently selected display is affected.
; RESTRICTIONS:
;	None that are obvious.
; PROCEDURE:
;	If the device has scalable pixels then the image is written over
;	the plot window.
; MODIFICATION HISTORY:
;	DMS, May, 1988.
;-

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour
contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px(1)-px(0)		;Size in x in device units
swy = py(1)-py(0)		;Size in Y
six = float(sz(1))		;Image sizes
siy = float(sz(2))
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
	endif

  tvscl,a,px(0),py(0),xsize = swx, ysize = swy, /device

endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tvscl,a,px(0),py(0)	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
	tvscl,poly_2d(bytscl(a),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px(0),py(0)
	endelse			;window_scale
  endelse			;scalable pixels

mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
m = max(a)
contour,a,ax,ay,/noerase,/nodata,/yst,$	;Do the contour
	  pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
          c_color =  colors, charsize=1.1, $
          xtitle='x (km)', ytitle='z (km)', title = tit, $
;          levels=[0.02*m,0.04*m,0.06*m,0.2*m,0.4*m,0.5*m,0.7*m,0.9*m],
	  xrange=[min(ax),max(ax)], xstyle=1
;dx = !x.crange(1) - !x.crange(0)
;dy = !y.crange(1) - !y.crange(0)
;xyouts,0.05*dx+!x.crange(0),0.85*dy+!y.crange(0),tit,/data

return
end
;---------------------------------------------------------------------------;

