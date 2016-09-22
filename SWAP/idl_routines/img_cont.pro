;---------------------------------------------------------------------;
;  NAME:  Peter                                                       ;
;  DATE:  7/20/93                                                     ;
;  FILE:  xpd_plot.pro                                                ;
;  SUMY:  image cont EXPAND.C simulation results...                   ;
;---------------------------------------------------------------------;

pro img_cont_xy, a, WINDOW_SCALE = window_scale, $
                ASPECT = aspect, INTERP = interp,POSTSCRIPT=postscript
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

nclr = !d.n_colors

if keyword_set(postscript) then begin
set_plot,'ps
device,filename='img.eps
device,/encapsulated
!p.font=0
device,/palatino
device,bits=8
device,/color
endif

f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz
;ax = y
;az = z

maxa=max(a)
mina=min(a)
print,!d.n_colors,nclr
a = (not(bytscl(a)))*nclr/256
print,max(a)
sz=size(a)
;ab = bytscl(a)
;a=rebin(a,4*sz(1),4*sz(2))
;sz=size(a)
 
!x.range=[min(y),max(y)]
!y.range=[min(z),max(z)]
;!x.title = 'y (10!u3!n km)' 
;!y.title = 'z (10!u3!n km)'
!x.title = 'y (km)' 
!y.title = 'z (km)'

!p.charsize=1.4
!y.margin = [4,4]
!x.margin = [10,10]
clrb = findgen(sz(2))*nclr/max(findgen(sz(2)))
;clrb = findgen(sz(2))*255/max(findgen(sz(2)))
for i=sz(1)-15,sz(1)-1 do begin
   a(i,*) = max(clrb)-clrb
;   a(i,*) = clrb
endfor

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour
contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1

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
  print,max(a)
  tv,a,px(0),py(0),xsize = swx, ysize = swy, /device

endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tv,a,px(0),py(0)	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
        print,max(a)
	tv,poly_2d((a),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px(0),py(0)
	endelse			;window_scale
  endelse			;scalable pixels

     axis,xaxis=2,xstyle=1
     axis,yaxis=2,ystyle=1
     axis,yaxis=1,ystyle=1,$
;     yticks=2,$
;     ytickv=[min(y),(min(y)+max(y))/2,max(y)], $
;     ytickname=[string(min(a)),string((min(a)+max(a))/2),string(max(a))], $
;     ytitle='E!dz!n (mV/m)',$
;     yrange = [mina*(2.3e-25/1.6e-19)*1e6,maxa*(2.3e-25/1.6e-19*1e6)]
     ytitle='(u!de!n)!dz!n (km/s)',$
     yrange = [mina,maxa]


mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
m = max(a)
;contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1
;contour,a,/noerase,/nodata,/yst,$	;Do the contour
;	  pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
;          c_color =  colors, charsize=1.5, $
;          xtitle='x (km)', ytitle='y (km)', $
;          levels=[0.02*m,0.04*m,0.06*m,0.2*m,0.4*m,0.5*m,0.7*m,0.9*m]
;;	  xrange=[min(ax),max(ax)], xstyle=1
;dx = !x.crange(1) - !x.crange(0)
;dy = !y.crange(1) - !y.crange(0)
;xyouts,0.05*dx+!x.crange(0),0.85*dy+!y.crange(0),tit,/data

if keyword_set(postscript) then begin
device,/close
set_plot,'x
!p.font=-1
endif

return
end
;---------------------------------------------------------------------------;


;---------------------------------------------------------------------------;
PRO img_cont_xy_m,nfrm,nfiles,t0,dt,toffset
;t0 is the ionization turn on time
;dt is the simulation time step
;toffset is the time used to locate the initial burst on the simulation grid
;---------------------------------------------------------------------------;
loadct,5
vsat=9.6
header=systime(0)
;!p.multi=[0,2,2]

f_read_coord,'coord.dat',x,y,z,dz_grid,dz_cell,nx,ny,nz

x=x+vsat*(toffset)
print,x
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
maxarr = fltarr(nfrm*nfiles)
icld = fltarr(nx,ny,nz,/nozero)
ncld = fltarr(nx,ny,nz,/nozero) 

frmcount=0
frame=0
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
      frame = [frame,frm]
      print,'   image #.....',frm
      readu,1,icld
;      readu,1,ncld
      for i=0,nx-1 do begin
         for j=0,ny-1 do begin
;            arr(i,j,frmcount) = total(icld(i,j,*)+ncld(i,j,*))
            arr(i,j,frmcount) = total(icld(i,j,*))
         endfor
      endfor
      maxarr(frmcount) = max(arr(*,*,frmcount))
      arr(*,*,frmcount) = bytscl(arr(*,*,frmcount))
      frmcount=frmcount+1
   endfor

   close,1
endfor

frame = frame*dt + t0

;arr = bytscl(arr)
arr = 255b - arr
;wh=where(arr eq 0)
;arr(wh) = 255b

for i=0,frmcount-1 do begin
      tit = strmid(strtrim(string(frame(i+1)),1),0,3)+' s'  
      con_arr = smooth(congrid(arr(*,*,i),nx*3,ny*3),4)
;      wh=where(con_arr eq 0)
;      con_arr(wh) = 255b
      conx = congrid(x,nx*3)
      cony = congrid(y,ny*3)
      image_cont_xy,con_arr,conx,cony,tit,/window_scale
;	surface,con_arr
endfor

set_plot,'ps
device,/close
file='np.ps'
print,file
device,filename=file 
device,/color,/times,font_size=12
device,/portrait,/inches,xsize=6.5,ysize=9.0,xoffset=1.0,yoffset=1.0
loadct,5
stretch,0,240
flg=8
filecnt=0

!p.multi=[0,2,4]
!p.font=7
for i=0,frmcount-1 do begin
;   if (flg eq 8) then begin
;      device,/close
;      flg=0
;      file='np'+strtrim(string(filecnt),1)+'.ps'
;      print,file
;      device,filename=file 
;      device,/color,/times,font_size=10
;      device,/portrait,/inches,xsize=6.0,ysize=9.0,xoffset=0.75,yoffset=1.0
;      filecnt=filecnt+1
;   endif
   tit = strmid(strtrim(string(frame(i+1)),1),0,3)+' s'  
   con_arr = smooth(congrid(arr(*,*,i),nx*4,ny*4),4)
;   wh=where(con_arr eq 0)
;   con_arr(wh) = 255b
   sz = size(con_arr)
   for j=0,sz(2)-1 do begin
      con_arr(sz(1)-sz(1)/40:*,j) = !d.n_colors - j*!d.n_colors/(sz(2)-1)
   endfor
   conx = congrid(x,nx*4)
   cony = congrid(y,ny*4)

   image_cont_xy,con_arr,conx,cony,tit,/window_scale,/aspect
   y1vals = [string(0.0),strmid(strtrim(string(maxarr(i)/max(maxarr)),1),0,3)]
   axis, yaxis=1, yrange = [0,maxarr(i)/max(maxarr)],yticks=1,ystyle=1, $
         ytickv = y1vals, charsize = 1.5

   xyouts,0.6*!d.x_size,1.05*!d.y_size,header,color=0,/device,charsize=1.0
;   flg = flg+1
endfor

device,/close
set_plot,'x

end
;--------------------------------------------------------------------------;









