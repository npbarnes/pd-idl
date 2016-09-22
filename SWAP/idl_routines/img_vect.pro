;-------------------------------------------------------------------
PRO img_velovec,a,U,V,X,Y, Missing = Missing, Length = length, $ 
                Dots = dots, Color=color, _EXTRA = extra
;-------------------------------------------------------------------
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

        on_error,2                      ;Return to caller if an error occurs
        s = size(u)
        t = size(v)
        if s(0) ne 2 then begin 
baduv:   message, 'U and V parameters must be 2D and same size.'
                endif
        if total(abs(s(0:2)-t(0:2))) ne 0 then goto,baduv
;
        if n_params(0) lt 3 then x = findgen(s(1)) else $
                if n_elements(x) ne s(1) then begin
badxy:                  message, 'X and Y arrays have incorrect size.'
                        endif
        if n_params(1) lt 4 then y = findgen(s(2)) else $
                if n_elements(y) ne s(2) then goto,badxy
;
        if n_elements(missing) le 0 then missing = 1.0e30
        if n_elements(length) le 0 then length = 1.0

        mag = sqrt(u^2+v^2)             ;magnitude.
                ;Subscripts of good elements
        nbad = 0                        ;# of missing points
        if n_elements(missing) gt 0 then begin
                good = where(mag lt missing) 
                if keyword_set(dots) then bad = where(mag ge missing, nbad)
        endif else begin
                good = lindgen(n_elements(mag))
        endelse

        ugood = u(good)
        vgood = v(good)
        x0 = min(x)                     ;get scaling
        x1 = max(x)
        y0 = min(y)
        y1 = max(y)
	x_step=float(x1-x0)/float(s(1))   ; Convert to float. Integer math
	y_step=float(y1-y0)/float(s(2))   ; could result in divide by 0

	maxmag=max([abs(max(ugood/x_step)),abs(max(vgood/y_step))])
	sina = length * (ugood/maxmag)
	cosa = length * (vgood/maxmag)
;
        if n_elements(title) le 0 then title = ''
        ;--------------  plot to get axes  ---------------
        if n_elements(color) eq 0 then color = !p.color
        x_b0=x0-x_step
	x_b1=x1+x_step
	y_b0=y0-y_step
	y_b1=y1+y_step
        if n_elements(position) eq 0 then begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/noerase,/nodata,/xst,/yst, $
            color=color, _EXTRA = extra, $
            pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev
        endif else begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/noerase,/nodata,/xst,/yst, $
            color=color, _EXTRA = extra, $
            pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev
        endelse
;
        r = .3                          ;len of arrow head
        angle = 22.5 * !dtor            ;Angle of arrowhead
        st = r * sin(angle)             ;sin 22.5 degs * length of head
        ct = r * cos(angle)
;
        for i=0,n_elements(good)-1 do begin     ;Each point
                x0 = x(good(i) mod s(1))        ;get coords of start & end
                dx = sina(i)
                x1 = x0 + dx
                y0 = y(good(i) / s(1))
                dy = cosa(i)
                y1 = y0 + dy
		xd=x_step
		yd=y_step
                plots,[x0,x1,x1-(ct*dx/xd-st*dy/yd)*xd, $
			x1,x1-(ct*dx/xd+st*dy/yd)*xd], $
                      [y0,y1,y1-(ct*dy/yd+st*dx/xd)*yd, $
			y1,y1-(ct*dy/yd-st*dx/xd)*yd], $
                      color=color
                endfor
        if nbad gt 0 then $             ;Dots for missing?
                oplot, x(bad mod s(1)), y(bad / s(1)), psym=3, color=color
end
;-------------------------------------------------------------------










