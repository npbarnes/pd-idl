;-------------------------------------------------------------------
pro read_np
;-------------------------------------------------------------------
@common

close,1
ntypefile = ''
print,'ntype...',ntype
if (ntype eq 0) then begin
   ntypefile = 'nfall_'
endif 
if (ntype eq 1) then begin
   ntypefile = 'npall_'
endif
if (ntype eq 2) then begin
   ntypefile = 'pf_'
endif
if (ntype eq 3) then begin
   ntypefile = 't1_'
endif
if (ntype eq 4) then begin
   ntypefile = 't2_'
endif
if (ntype eq 5) then begin
   ntypefile = 't3_'
endif
openr,1,ntypefile+'1.dat',/f77_unformatted

frame=0l
nt=0l
nout=0l
nx=0l
ny=0l
nz=0l
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

close,1

icld = fltarr(nx,ny,nz,/nozero)
ncld = fltarr(nx,ny,nz,/nozero) 

mfile = ((frm-1)/(nfrm/nfiles)) + 1
;print,mfile,frm
frmcount=1
frmcount=1+((frm-1) mod (nfrm/nfiles))
openr,1,ntypefile+strtrim(string(mfile),1)+'.dat',/f77_unformatted
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

cnt=1
while not(eof(1)) do begin
   readu,1,frame
   print, ntypefile+strtrim(string(mfile),1)+' image #.....',frame
   readu,1,icld
   if (cnt eq frmcount) then goto, BAIL
   cnt=cnt+1
endwhile

close,1

;frmcount=0
;for m=1,nfiles do begin

;   npfile = 'npall_'+strtrim(string(m),1)+'.dat'
;   print,'Reading....',npfile

;   openr,1,npfile,/f77_unformatted

;   readu,1,nt
;   readu,1,nout
;   readu,1,nx
;   readu,1,ny
;   readu,1,nz

;   for n = 1,nfrm/nfiles do begin
;      readu,1,frame
;      print,'   image #.....',frame
;      readu,1,icld
;;      readu,1,ncld
;      if (frmcount eq frm-1) then goto, BAIL    
;      frmcount=frmcount+1
;   endfor

;   close,1
;endfor

BAIL:  

tot_icld = icld
icld = icld(minx:maxx,miny:maxy,minz:maxz)

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro read_file,flg
;-------------------------------------------------------------------
@common

close,1

frame=0l
nt=0l
nout=0l
nx=0l
ny=0l
nz=0l

mfile = ((frm-1)/(nfrm/nfiles)) + 1
print,mfile,frm
frmcount=1
frmcount=1+((frm-1) mod (nfrm/nfiles))
openr,1,file(0)+strtrim(string(mfile),1)+'.dat',/f77_unformatted
readu,1,nt
print,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz
print,nt,nout,nx,ny,nz

temparr = fltarr(nx,ny,nz,3,/nozero)
cnt=1
while not(eof(1)) do begin
   readu,1,frame
   print,file(0)+strtrim(string(mfile),1)+' image #.....',frame
   readu,1,temparr
   if (cnt eq frmcount) then goto, BAIL
   cnt=cnt+1
endwhile

close,1


;for m=1,nfiles do begin
   
;   openr,1,file(0)+strtrim(string(m),1)+'.dat',/f77_unformatted

;   readu,1,nt
;   readu,1,nout
;   readu,1,nx
;   readu,1,ny
;   readu,1,nz

;   temparr = fltarr(nx,ny,nz,3,/nozero)

;;   for i = 0,frm-1 do begin
;;      readu,1,frame
;;      print,'image #.....',frame
;;      readu,1,temparr
;;   endfor

;   while not(eof(1)) do begin
;      readu,1,frame
;      print,'image #.....',frame
;      readu,1,temparr
;      if (frmcount eq frm-1) then goto, BAIL
;      frmcount = frmcount+1
;   endwhile

;   close,1

;endfor

BAIL: 

tot_arr=temparr
if (flg eq 0) then begin
    arr = temparr 
endif else begin
   arr = temparr(minx:maxx,miny:maxy,minz:maxz,*)
endelse   

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_surface
;-------------------------------------------------------------------
@common
;!p.font=-1
case pln of

   'xy': begin
	 !x.title='x'
	 !y.title='y'
         surface,reform(arr(*,*,slcxy,vcomp)),charsize=1.0
    	 end;'xy'

   'xz': begin
	 !x.title='x'
	 !y.title='z'
         surface,reform(arr(*,slcxz,*,vcomp)),charsize=1.0
	 !x.title='x'
	 !y.title='y'
	 end;'xz'

   'yz': begin
	 !x.title='y'
	 !y.title='z'
         surface,reform(arr(slcyz,*,*,vcomp)),charsize=1.0
	 !x.title='x'
	 !y.title='y'
	 end;'yz'
endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_contour
;-------------------------------------------------------------------
@common

f_read_coord,'coord.dat',x,y,z,dzc,dzg,nx,ny,nz
x = x(minx:maxx)
y = y(miny:maxy)
z = z(minz:maxz)

zbin = (max(z)-min(z))/n_elements(z)
xbin = (max(x)-min(x))/n_elements(x)
ybin = (max(y)-min(y))/n_elements(y)

case pln of

   'xy': begin
	 !x.title = 'x (km)'
	 !y.title = 'y (km)'
         ar = reform(arr(*,*,slcxy,vcomp))
         contour,ar, $
                 charsize=1.0,nlevels=nlev; $
;                 xrange=[min(x),max(x)],xsty=1,yrange=[min(y),max(y)],ysty=1
;                 xtitle = !x.title + '('+strtrim(string(xbin),2)+'  km)', $
;                 ytitle = !y.title + '('+strtrim(string(ybin),2)+'  km)'

         end;'xy'
    
   'xz': begin
         !x.title='x (km)'
	 !y.title='z (km)'
         ar=reform(arr(*,slcxz,*,vcomp))
         regrid_xz,ar,xx,zx
         print,'array sizes...',size(arr)
         print,size(ar)
         print,size(xx)
	 print,size(zx)
         contour,ar,xx,zx,$
                 charsize=1.0,nlevels=nlev, $
                 xrange=[min(x),max(x)],xsty=1,yrange=[min(z),max(z)],ysty=1
;                 xtitle = !x.title + ' (km)', $
;                 ytitle = !y.title + ' (km)', $

	 !x.title='x'
	 !y.title='y'
         end;'xz'

   'yz': begin
	 !x.title='y (km)'
	 !y.title='z (km)'
         ar=reform(arr(slcyz,*,*,vcomp))
         regrid_yz,ar,yy,zy
         contour,ar,yy,zy,$
                 charsize=1.0,nlevels=nlev, $
                 xrange=[min(y),max(y)],xsty=1,yrange=[min(z),max(z)],ysty=1
;                 xtitle = !x.title + ' (km)', $
;                 ytitle = !y.title + ' (km)', $

	 !x.title='x'
	 !y.title='y'
	 end;'yz'

endcase

return
end
;-------------------------------------------------------------------

;-------------------------------------------------------------------
PRO img_velovec,img,U,V,X,Y, Missing = Missing, Length = length, $ 
                Dots = dots, Color=color, _EXTRA = extra, $
                WINDOW_SCALE = window_scale, ASPECT = aspect, INTERP = interp
;-------------------------------------------------------------------
@common

;a = reverse(not(bytscl(img)))
;;a = reverse(a,2)
;U = reverse(smooth(-U,2))
;V = reverse(smooth(V,2))

a = (bytscl(img))
;a = reverse(a,2)
U = (smooth(U,2))
V = (smooth(V,2))


;U = reverse(U,2)
;V = reverse(-V,2)

sz=size(a)
ab = bytscl(a)
a=rebin(a,6*sz(1),6*sz(2))
sz=size(a)
 

wpi = 2200.0
X = (X - X(15))/wpi
Y = (Y - Y(58/2))/wpi


;!x.range=[min(X),max(X)]
;!y.range=[min(Y),max(Y)]
;!x.margin=[10,4]
;!y.margin=[10,4]
;!x.title = 'y (10!u3!n km)' 
;!y.title = 'z (10!u3!n km)'
;!x.title = 'y' 
;!y.title = 'z'

case pln of
   'xy': begin
         !x.title = 'x (km)' 
         !y.title='y (km)'
         end
   'xz': begin
         !x.title = 'x (c/wpi)'
         !y.title='z (c/wpi)'
         end
   'yz': begin
         !x.title = 'y (km)'
         !y.title='z (km)'
         end
endcase

;!y.margin = [4,9]
;!x.margin = [6,9]

clrb = findgen(sz(2))*255/max(findgen(sz(2)))
for i=sz(1)-3,sz(1)-1 do begin
   a(i,*) = clrb
endfor

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour
contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4,xmargin=[10,10],$
            xrange=[min(X),max(X)],yrange=[min(Y),max(Y)]

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

  tv,a*numclr/max(a),px(0),py(0),xsize = swx, ysize = swy, /device

endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tv,a*numclr/max(a),px(0),py(0)	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
	tv,poly_2d(a*numclr/max(a),$	;Have to resample image
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
        if n_params(0) lt 4 then x = findgen(s(1)) else $
                if n_elements(x) ne s(1) then begin
badxy:                  message, 'X and Y arrays have incorrect size.'
                        endif
        if n_params(1) lt 5 then y = findgen(s(2)) else $
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

	maxmag=max([abs(max(abs(ugood/x_step))),abs(max(abs(vgood/y_step)))])
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
          plot,[x_b0,x_b1],[y_b1,y_b0],/noerase,/nodata,/xst,yst=9, $
            color=color, _EXTRA = extra,$
            pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev
        endif else begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/noerase,/nodata,/xst,yst=9, $
            color=color, _EXTRA = extra,$
            pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev
        endelse
     
     axis,yaxis=1,ystyle=1,$
;     yticks=2,$
;     ytickv=[min(y),(min(y)+max(y))/2,max(y)], $
;     ytickname=[string(min(a)),string((min(a)+max(a))/2),string(max(a))], $
;     ytitle='Density (ions/cm!u3!n)',$
     ytitle='Density (km^-3)',$
     yrange = [min(img),max(img)]



;
        r = .5                          ;len of arrow head
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
;                      color=(ab(x0,y0) + vclr ) mod !d.n_colors 
;                       color = abs(vclr-ab(x0,y0))
                       color = 0
                endfor
        if nbad gt 0 then $             ;Dots for missing?
                oplot, x(bad mod s(1)), y(bad / s(1)), psym=3, color=vclr
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO VELOVEC,U,V,X,Y, Missing = Missing, Length = length, Dots = dots,  $
        Color=color, _EXTRA = extra
;-------------------------------------------------------------------

U = smooth(U,2,/edge_truncate)
V = smooth(V,2,/edge_truncate)

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

	maxmag=max([abs(max(abs(ugood/x_step))),abs(max(abs(vgood/y_step)))])
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
          plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
            color=color, _EXTRA = extra
        endif else begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
            color=color, _EXTRA = extra
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


;--------------------------------------------------------------------------
PRO make_vector
;--------------------------------------------------------------------------
@common

case pln of

   'xy': begin

      vxy = reform(arr(*,*,slcxy,0))
      vyx = reform(arr(*,*,slcxy,1))

      if not((abs(max(vxy)) eq 0) and (abs(max(vyx)) eq 0)) then begin
         !x.title='x'
         !y.title='y'
         !p.title=!p.title + ' (max = '+strtrim(string(max(sqrt(vxy^2 + vyx^2))),2)+')'
         velovec,vxy,vyx,length=1.0, charsize=1.0
;                 title='max ='+string(max(sqrt(vxy^2 + vyx^2)))
         !p.title=''
      endif
      end;xy

   'xz': begin
 
      vxz = reform(arr(*,slcxz,*,0))
      vzx = reform(arr(*,slcxz,*,2))

      if not((abs(max(vxz)) eq 0) and (abs(max(vzx)) eq 0)) then begin
         !x.title='x'
         !y.title='z'
         !p.title=!p.title + ' (max = '+strtrim(string(max(sqrt(vxz^2 + vzx^2))),2)+')'
         velovec,vxz,vzx,length=1.0,charsize=1.0
;                 title='max ='+string(max(sqrt(vxz^2 + vzx^2))), $
         !x.title='x'
         !y.title='y'
         !p.title=''
      endif

      end;'xz'

   'yz': begin

      vyz = reform(arr(slcyz,*,*,1))
      vzy = reform(arr(slcyz,*,*,2))

      if not((abs(max(vyz)) eq 0) and (abs(max(vzy)) eq 0)) then begin
         !x.title='y'
         !y.title='z'
         !p.title=!p.title + ' (max = '+strtrim(string(max(sqrt(vyz^2 + vzy^2))),2)+')'
         velovec,vyz,vzy,length=1.0, charsize=1.0
;                 title='max ='+string(max(sqrt(vyz^2 + vzy^2))), $
         !x.title='x'
         !y.title='y'
         !p.title=''
      endif
      end;yz

endcase

return
end
;------------------------------------------------------------------------


;--------------------------------------------------------------------------
PRO make_img_vector
;--------------------------------------------------------------------------
@common

f_read_coord,'coord.dat',x,y,z,dzc,dzg,nx,ny,nz
x = x(minx:maxx)
y = y(miny:maxy)
z = z(minz:maxz)

case pln of

   'xy': begin
      erase
      img = reform(icld(*,*,slcxy))
      vxy = reform(arr(*,*,slcxy,0))
      vyx = reform(arr(*,*,slcxy,1))

      if not((abs(max(vxy)) eq 0) and (abs(max(vyx)) eq 0)) then begin
         !x.title='x'
         !y.title='y'  
         ptit = !p.title
         !p.title=''
;         !p.title=!p.title + ' (max = '+strtrim(string(max(sqrt(vxy^2 + vyx^2))),2)+')'
         img_velovec,img,vxy,vyx,x,y,length=1.0,charsize=1.0,/aspect,$
              title=ptit+' (max = '+strtrim(string(max(sqrt(vxy^2 + vyx^2))),2)+')'
      endif
      end;xy

   'xz': begin
      erase
      img = reform(icld(*,slcxz,*))
      vxz = reform(arr(*,slcxz,*,0))
      vzx = reform(arr(*,slcxz,*,2))

      if not((abs(max(vxz)) eq 0) and (abs(max(vzx)) eq 0)) then begin
         !x.title='x'
         !y.title='z'
         ptit = !p.title
         !p.title=''
;         !p.title=!p.title + ' (max = '+strtrim(string(max(sqrt(vxz^2 + vzx^2))),2)+')'
         img_velovec,img,vxz,vzx,x,z,length=1.0,charsize=1.0,/aspect, $
              title=ptit+' (max = '+strtrim(string(max(sqrt(vxz^2 + vzx^2))),2)+')'
         !x.title='x'
         !y.title='y'
      endif
      end;'xz'

   'yz': begin
      erase
      img = reform(icld(slcyz,*,*))
      vyz = reform(arr(slcyz,*,*,1))
      vzy = reform(arr(slcyz,*,*,2))

      if not((abs(max(vyz)) eq 0) and (abs(max(vzy)) eq 0)) then begin
         !x.title='y'
         !y.title='z'
	 ptit = !p.title
	 !p.title=''
;         !p.title=!p.title + ' (max = '+strtrim(string(max(sqrt(vyz^2 + vzy^2))),2)+')'
         img_velovec,img,vyz,vzy,y,z,length=1.0,charsize=1.0,/aspect, $
              title=ptit+' (max = '+strmid(strtrim(string(max(sqrt(vyz^2 + vzy^2))),2),0,4)+')'
         !x.title='x'
         !y.title='y'
      endif
      end;yz

endcase

return
end
;------------------------------------------------------------------------


;--------------------------------------------------------------------------
PRO make_profile
;--------------------------------------------------------------------------
@common

case pln of

   'xy': begin
      erase
      img = reform(icld(*,*,slcxy))
      sz = size(img)
      img=rebin(img,10*sz(1),10*sz(2))
      tvscl,img
      profiles,img
      end;xy

   'xz': begin
      erase
      img = reform(icld(*,slcxz,*))
      sz = size(img)
      img=rebin(img,10*sz(1),10*sz(2))
      tvscl,img
      profiles,img
      end;'xz'

   'yz': begin
      erase
      img = reform(icld(slcyz,*,*))
      sz = size(img)
      img=rebin(img,10*sz(1),10*sz(2))
      tvscl,img
      profiles,img
      end;yz

endcase

return
end
;------------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_btrace
;-------------------------------------------------------------------
@common 

if (file ne 'b1all_') then begin
   file = 'b1all_'
   read_file,1
endif

f_read_coord,'coord.dat',x,y,z,dzg,dzc,nnx,nny,nnz

case pln of

   'xy': begin
      erase
      img = reform(icld(*,*,slcxy))
      bx = reform(arr(*,*,slcxy,0))
      bz = reform(arr(*,*,slcxy,1))
      nxx = nx
      nzz = ny
      xx = x
      zz = y
      !x.title='x'
      !y.title='y'
      end;xy

   'xz': begin
      erase
      img = reform(icld(*,slcxz,*))
      bx = reform(arr(*,slcxz,*,0))
      bz = reform(arr(*,slcxz,*,2))
      nxx = nx
      nzz = nz
      xx = x
      zz = z
      !x.title='x'
      !y.title='z'
      end;'xz'

   'yz': begin
      erase
      img = reform(icld(slcyz,*,*))
      bx = reform(arr(slcyz,*,*,1))
      bz = reform(arr(slcyz,*,*,2))
      nxx = ny
      nzz = nz
      xx = y
      zz = z
      !x.title='y'
      !y.title='z'
      end;yz

endcase

dx = max(xx)/nxx
dz = max(zz)/nzz

print,dx,dz
;regrid,bx
;regrid,bz

;image_cont,img
contour,bx,/nodata,xrange=[0,nxx*dx],yrange=[0,nzz*dz],xstyle=1,ystyle=1,$
        xtitle = 'x (km)', ytitle = 'z (km)'

x0 = 0.0

repeat begin

x0 = x0 + 20.0
x1 = x0
z1 = 0.0
flg = 0

while (flg ne 1) do begin
  i = round(nxx*x1/max(xx))
  k = round(nzz*z1/max(zz))

  if (i gt nxx-1) then goto, BAIL
  if (i lt 0) then goto, BAIL
  if (k gt nzz-1) then goto, BAIL
  if (k lt 0) then goto, BAIL

  bx1 = bx(i,k)
  bz1 = bz(i,k) + 2e-5*1.6e-19/2.3e-25
  z2 = z1+dz
  x2 = x1+dz*bx1/abs(bz1)
  if ((x2-x1) gt dx) then begin
	x2 = x1 + dx
	z2 = z1 + dz*abs(bz1)/bx1
	endif
  plots,[x1,x2],[z1,z2],/data
;  print,x1,x2,z1,z2
  x1 = x2
  z1 = z2

  if (x1 lt 0) then flg = 1
  if (x1 gt max(xx)) then flg = 1

  if (z1 lt 0) then flg = 1
  if (z1 gt max(zz)) then flg = 1



endwhile

BAIL:

endrep until(x1 gt max(xx))

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_movie_frames
;-------------------------------------------------------------------
@common

close,1
files = file+strtrim(string(1),1)+'.dat'
openr,1,files,/f77_unformatted

nt=0l
nout=0l
nx=0l
ny=0l
nz=0l

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

close,1

case pln of
   'xy': anmt_arr = fltarr(nx,ny,nfrm)
   'xz': anmt_arr = fltarr(nx,nz,nfrm)
   'yz': anmt_arr = fltarr(ny,nz,nfrm)
endcase

temparr = fltarr(nx,ny,nz,3,/nozero)

i=0
for m=1,nfiles do begin

   files = file+strtrim(string(m),1)+'.dat'
   openr,1,files,/f77_unformatted
   print,'Reading file.....',files

   readu,1,nt
   readu,1,nout
   readu,1,nx
   readu,1,ny
   readu,1,nz

   while (not(eof(1))) do begin
      readu,1,frame
      readu,1,temparr
      case pln of
         'xy': anmt_arr(*,*,i) = reform(temparr(*,*,slcxy,vcomp))
         'xz': begin
               ar = reform(temparr(*,slcxz,*,vcomp))
               regrid_xz,ar,xx,xz
               anmt_arr(*,*,i) = ar
               end
         'yz': anmt_arr(*,*,i) = reform(temparr(slcyz,*,*,vcomp))
      endcase
      i = i + 1
   endwhile
   close,1

endfor

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_animate
;-------------------------------------------------------------------
@common

case pln of

   'xy': begin
      get_movie_frames
      arr1=bytscl(anmt_arr)
      arr2=congrid(arr1,500,500,nfrm)
      xinteranimate,set=[500,500,nfrm]
      for i=0,nfrm-1 do begin
         xinteranimate,image=arr2(*,*,i),frame=i
      endfor
      xinteranimate
      end;xy
    
   'xz': begin
      get_movie_frames
      arr1=bytscl(anmt_arr)
      arr2=congrid(arr1,500,500,nfrm)
      xinteranimate,set=[500,500,nfrm]
      for i=0,nfrm-1 do begin
         xinteranimate,image=arr2(*,*,i),frame=i
      endfor
      xinteranimate
      end;xz
        
   'yz': begin
      get_movie_frames
      arr1=bytscl(anmt_arr)
      arr2=congrid(arr1,500,500,nfrm)
      xinteranimate,set=[500,500,nfrm]
      for i=0,nfrm-1 do begin
         xinteranimate,image=arr2(*,*,i),frame=i
      endfor
      xinteranimate
      end;yz

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_graphics
;-------------------------------------------------------------------
@common

case graph_type of

   'surf': begin
      make_surface
      end;surf

   'cont': begin
      make_contour
      end;cont

   'anim': begin
      make_animate
      end;anim

   'vect': begin
      make_vector
      end;vect

   'img_vect': begin
      make_img_vector
      end;img_vect

   'profile': begin
      make_profile
      end;profile

   'btrc': begin
      make_btrace
      end;btrc

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro view_event,ev
;-------------------------------------------------------------------
@common

widget_control, ev.id, get_uvalue = uvalue

case uvalue of

   'X': begin
      vcomp=0
      if not(ev.select eq 0) then get_graphics
      end;X

   'Y': begin
      vcomp=1
      if not(ev.select eq 0) then get_graphics
      end;Y

   'Z': begin
      vcomp=2
      if not(ev.select eq 0) then get_graphics
      end;Z

   'XY': begin
      widget_control, ev.id, get_value = plane
      if (slcxy gt maxz-minz) then slcxy = maxz-minz-1
      pln=plane
      widget_control,wxy,sensitive=1
      widget_control,wxz,sensitive=0
      widget_control,wyz,sensitive=0
      if not(ev.select eq 0) then get_graphics
      end;XY

   'XZ': begin
      widget_control, ev.id, get_value = plane
      if (slcxz gt maxy-miny) then slcxz = maxy-miny-1
      pln=plane
      widget_control,wxy,sensitive=0
      widget_control,wxz,sensitive=1
      widget_control,wyz,sensitive=0
      if not(ev.select eq 0) then get_graphics
      end;XZ

   'YZ': begin
      widget_control, ev.id, get_value = plane
      if (slcyz gt maxx-minx) then slcyz = maxx-minx-1
      pln=plane
      widget_control,wxy,sensitive=0
      widget_control,wxz,sensitive=0
      widget_control,wyz,sensitive=1
      if not(ev.select eq 0) then get_graphics
      end;YZ

   'SLICEXY': begin
      widget_control, ev.id, get_value = slice
      slcxy=slice-minz
      get_graphics
      end;SLICE

   'SLICEXZ': begin
      widget_control, ev.id, get_value = slice
      slcxz=slice-miny
      get_graphics
      end;SLICE

   'SLICEYZ': begin
      widget_control, ev.id, get_value = slice
      slcyz=slice-minx
      get_graphics
      end;SLICE

   'FRAME': begin
      widget_control, ev.id, get_value = frame
      frm=frame
      print,'vecscl...',vecscl
      if (vecscl eq 1) then begin
         print,'reading vector...'
         read_file,1
         read_np
      endif
      if (vecscl eq 0) then begin
	 print,'reading scalar...'
         read_np
         arr = icld
      endif
      get_graphics
      end;FRAME

   'SURFACE': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=0
      graph_type = 'surf'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'surf'
      end;SURFACE

   'CONTOUR': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=1
      graph_type = 'cont'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'cont'
      end;CONTOUR

   'VECTOR': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=0
      graph_type = 'vect'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'vect'
      end;VECTOR

   'IMG_VECTOR': begin
      widget_control,wclr,sensitive=1
      widget_control,wnlev,sensitive=0
      graph_type = 'img_vect'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'img_vect'
      end;IMG_VECTOR

   'PROFILES': begin
      widget_control,wclr,sensitive=1
      widget_control,wnlev,sensitive=0
      graph_type = 'profile'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'profile'
      end;PROFILES

   'BTRACE': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=0
      graph_type = 'btrc'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'btrc'
      end;BTRACE


   'ANIMATE': begin 
      graph_type = 'anim'
      get_graphics
      graph_type = graph_type_prev
      end;ANIMATE
  
   'VIEWXY': partmovie_xy_m,nfrm/nfiles,nfiles
   'VIEWXZ': partmovie_xz_m,nfrm/nfiles,nfiles
   'VIEWYZ': partmovie_yz_m,nfrm/nfiles,nfiles
   
   'CTABLES':  xloadct

   'VEC_COLOR': begin
      widget_control, ev.id, get_value = vclr
      graph_type = 'img_vect'
      get_graphics
      graph_type_prev = 'img_vect'
      end;VEC_COLOR

   'NLEV': begin
      widget_control, ev.id, get_value = nlev
      graph_type = 'cont'
      get_graphics
      graph_type_prev = 'cont'
      end;NLEV

   'B1': begin
      vecscl = 1
      file = 'b1all_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;B1

   'UF': begin 
      vecscl = 1
      file = 'ufall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UF

   'AJ': begin 
      vecscl = 1
      file = 'ajall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;AJ

   'UGRADU': begin
      vecscl = 1 
      file = 'ugraduall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;AJ
      
   'E' : begin
      vecscl = 1
      file = 'Eall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;E      

   'EF' : begin
      vecscl = 1
      file = 'Efall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;EF      

   'UP' : begin
      vecscl = 1
      file = 'upall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UP      

   'UI' : begin
      vecscl = 1
      file = 'uiall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UI      

   'UF2' : begin
      vecscl = 1
      file = 'uf2all_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UF2      

   'UFP2' : begin
      vecscl = 1
      file = 'ufp2all_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UFP2    
  
   'NP' : begin
      vecscl = 1
      ntype = 1
      read_np
      if not(ev.select eq 0) then get_graphics
      end;NP

   'NP1' : begin
      vecscl = 0
      ntype = 1
      read_np
      arr = icld
      if not(ev.select eq 0) then get_graphics
      end;NP1

   'NF' : begin
      vecscl = 1
      ntype = 0
      read_np
      if not(ev.select eq 0) then get_graphics
      end;NF

   'NF1' : begin
      vecscl = 0
      ntype = 0
      read_np
      arr = icld
      if not(ev.select eq 0) then get_graphics
      end;NF
   
   'PF' : begin
      vecscl = 1
      ntype = 0
      read_np
      nfarr = icld
      ntype = 2
      read_np
      icld = icld/(1.38e-29*nfarr*11600.)
      if not(ev.select eq 0) then get_graphics
      end;PF

   'PF1' : begin
      vecscl = 0
      ntype = 2
      read_np
      arr = icld
      if not(ev.select eq 0) then get_graphics
      end;PF

   'T1' : begin
      ntype = 3
      read_np
      if not(ev.select eq 0) then get_graphics
      end;T1

   'T2' : begin
      ntype = 4
      read_np
      if not(ev.select eq 0) then get_graphics
      end;T2

   'T3' : begin
      ntype = 5
      read_np
      if not(ev.select eq 0) then get_graphics
      end;T3

   'PATCH' : begin
      get_sub_arrays
      end;PATCH
   
   'CHOPPER' : begin
        COMMON VOLUME_DATA, A
	A = alog(icld)       
	SLICER
      end;CHOPPER

   'PS' : begin
      get_ps_options
      end;PS

   'DONE': widget_control,/destroy, ev.top

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_view_widget
;-------------------------------------------------------------------
@common

base = widget_base(title='CRRES Hybrid Code 3-D Viewer',/row)
lc = widget_base(base,/column,space=5)
mc = widget_base(base,/frame,/column)
mct = widget_base(mc,/row,space=5)
mctl = widget_base(mct,/row,/frame,space=5)
mctr = widget_base(mct,/row,/frame,space=5)
rc = widget_base(base,/column,space=5)
rct = widget_base(rc,/column,space=5)
rcanim = widget_base(rc,/frame,/column,space=5)
rcct = widget_base(rc,/column,space=5)

xmenu,['x','y','z'],lc,/column,/exclusive,$
      title='vector component',space=10,uvalue=['X','Y','Z']
xmenu,['xy','xz','yz'],lc,/column,/exclusive,$
      title='plane',space=10,uvalue=['XY','XZ','YZ']
xmenu,['surface','contour','vector','img_vector','np profile','B trace'],$
      rct,/column, $
      /exclusive,title='graphics options',space=10,$
      uvalue=['SURFACE','CONTOUR','VECTOR','IMG_VECTOR','PROFILES','BTRACE']

w1 = widget_label(rcanim,value = 'Animations')
w1 = widget_button(rcanim,value = 'fluid',uvalue = 'ANIMATE')
w1 = widget_button(rcanim,value = 'cloud xy',uvalue='VIEWXY')
w1 = widget_button(rcanim,value = 'cloud xz',uvalue='VIEWXZ')
w1 = widget_button(rcanim,value = 'cloud yz',uvalue='VIEWYZ')

w2 = widget_button(rcct,value = 'color tables',uvalue='CTABLES')

wclr = widget_slider(rct,title='vector color offset', uvalue='VEC_COLOR', $
                       maximum = 255, minimum=0)
widget_control,wclr,sensitive=0



wnlev = widget_slider(rct,title='# contour levels', uvalue='NLEV', $
                    maximum = 30, minimum = 1)
widget_control,wnlev,sensitive=0

wxy = widget_slider(lc,title='z slice',uvalue='SLICEXY', $
                         maximum=maxz, minimum = minz)
widget_control,wxy,sensitive=0
wxz = widget_slider(lc,title='y slice',uvalue='SLICEXZ', $
                         maximum=maxy, minimum=miny)
widget_control,wxz,sensitive=0
wyz = widget_slider(lc,title='x slice',uvalue='SLICEYZ', $
                         maximum=maxx, minimum=minx)
widget_control,wyz,sensitive=0
wfrm = widget_slider(lc,title='frame #',uvalue='FRAME', $
                       maximum=nfrm, minimum = 1)

xmenu,['b1','uf','aj','E','ugradu','up','np','nf','pf'],$
      mctl,/row,/exclusive,$
      space=1,uvalue=['B1','UF','AJ',$
                       'E','UGRADU','UP','NP1','NF1','PF1']

xmenu,['np','nf','pf','t1','t2','t3'],mctr,/row,/exclusive,$
      space=1,uvalue=['NP','NF','PF','T1','T2','T3']

w1 = widget_button(rcct,value = ' Extract Subarray ',uvalue='PATCH')
w1 = widget_button(rcct,value = ' Slicer', uvalue='CHOPPER')
w1 = widget_button(rcct,value = ' PostScript ',uvalue = 'PS')

draw = widget_draw(mc,xsize=700,ysize=600)

w1 = widget_button(mc,value = 'done',uvalue='DONE')

widget_control,/realize,base,/hourglass
widget_control,get_value = window,draw
wset,window

xmanager,'view',base

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO size_event,ev
;-------------------------------------------------------------------
@common

widget_control, ev.id, get_uvalue = uvalue
widget_control,ev.id,get_value = minmax

case uvalue of

   'RESET': begin
      minx=0 & maxx=nx-1
      miny=0 & maxy=ny-1
      minz=0 & maxz=nz-1
      widget_control,wmnx,/realize,set_slider_min=0
      widget_control,wmxx,/realize,set_slider_max=nx-1
      widget_control,wmny,/realize,set_slider_min=0
      widget_control,wmxy,/realize,set_slider_max=ny-1
      widget_control,wmnz,/realize,set_slider_min=0
      widget_control,wmxz,/realize,set_slider_max=nz-1
      widget_control,wxy,/realize,set_slider_min=0,set_slider_max=nz-1
      widget_control,wxz,/realize,set_slider_min=0,set_slider_max=ny-1
      widget_control,wyz,/realize,set_slider_min=0,set_slider_max=nx-1
      icld=tot_icld
      arr=tot_arr
      get_graphics
      end;RESET

   'XMIN': begin
      minx = minmax
      widget_control,wyz,/realize,set_slider_min=minx,set_slider_max=maxx
      if (slcyz lt minx) then slcyz = 0
      if (minx ge nx-3) then begin
         widget_control,wmxx,/realize,set_slider_min=nx-2
         minx=nx-3
         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
         get_graphics
      endif else begin
         widget_control,wmxx,/realize,set_slider_min=minx+1
         print,minx,maxx,miny,maxy,minz,maxz
         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
         get_graphics
      endelse
      end;XMIN
   'XMAX': begin 
      maxx = minmax
      widget_control,wyz,/realize,set_slider_min=minx,set_slider_max=maxx
      if (slcyz gt maxx-minx) then slcyz = maxx-minx-1
      widget_control,wmnx,set_slider_max=maxx-1
      icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
      arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
      get_graphics
      end;XMAX

   'YMIN': begin
      miny = minmax
      widget_control,wxz,/realize,set_slider_min=miny,set_slider_max=maxy
      if (slcxz lt miny) then slcxz = 0
      if (miny ge ny-3) then begin
         widget_control,wmxy,/realize,set_slider_min=ny-2
         miny=ny-3
         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
         get_graphics
      endif else begin
         widget_control,wmxy,/realize,set_slider_min=miny+1
         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
         get_graphics
      endelse
      end;YMIN
   'YMAX': begin
      maxy = minmax
      widget_control,wxz,/realize,set_slider_min=miny,set_slider_max=maxy
      if (slcxz gt maxy-miny) then slcxz = maxy-miny-1
      widget_control,wmny,set_slider_max=maxy-1
      icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
      arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
      get_graphics
      end;YMAX

   'ZMIN': begin
      minz = minmax
      widget_control,wxy,/realize,set_slider_min=minz,set_slider_max=maxz
      if (slcxy lt minz) then slcxy = 0
      if (minz ge nz-3) then begin
         widget_control,wmxz,/realize,set_slider_min=nz-2
         minz=nz-3
         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
         get_graphics
      endif else begin
         widget_control,wmxz,/realize,set_slider_min=minz+1
         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
         get_graphics
      endelse
      end;ZMIN
   'ZMAX': begin
      maxz = minmax
      widget_control,wxy,/realize,set_slider_min=minz,set_slider_max=maxz
      if (slcxy gt maxz-minz) then slcxy = maxz-minz-1
      widget_control,wmnz,set_slider_max=maxz-1
      icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
      arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
      get_graphics
      end;ZMAX

   'DONE': begin 
      widget_control,/destroy, ev.top
;      icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
;      arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
;      get_graphics
      end;DONE

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO get_sub_arrays
;-------------------------------------------------------------------
@common

wsizer = widget_base(title='Sizer',/column)

w9 = widget_button(wsizer,value = 'reset',uvalue='RESET')

wmnx = widget_slider(wsizer,title='min x',uvalue='XMIN', $
                         maximum=nx-1, minimum=0)
wmxx = widget_slider(wsizer,title='max x',uvalue='XMAX', $
                         maximum=nx-1, minimum=0)
wmny = widget_slider(wsizer,title='min y',uvalue='YMIN', $
                         maximum=ny-1, minimum=0)
wmxy = widget_slider(wsizer,title='max y',uvalue='YMAX', $
                         maximum=ny-1, minimum=0)
wmnz = widget_slider(wsizer,title='min z',uvalue='ZMIN', $
                         maximum=nz-1, minimum=0)
wmxz = widget_slider(wsizer,title='max z',uvalue='ZMAX', $
                         maximum=nz-1, minimum=0)
w8 = widget_button(wsizer,value = 'done',uvalue='DONE')

widget_control,/realize,wsizer,/hourglass
xmanager,'size',wsizer

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO post_event,ev
;-------------------------------------------------------------------
@common

widget_control, ev.id, get_uvalue = uvalue
widget_control,ev.id,get_value = label
label=label(0)
name=uvalue
val=label

case uvalue of

   'FILE': begin
       psfile = label 
       print,label
       end;FILE
   'XTIT': !x.title = label
   'YTIT': !y.title = label
   'ZTIT': !z.title = label
   'PTIT': !p.title = label
   'PORT': begin
       prt = 1
       lnd = 0
       end;PORT 
   'LAND': begin
       lnd = 1
       prt = 0
       end;LAND
   'YSIZE': begin
       ysz = label
       end;YSIZE
   'ENCAP': begin
       encp = 1
       noencp = 0
       end;ENCAP
   'NONENCAP': begin
       noencp = 1
       encp = 0
       end;NONENCAP
   'COLOR': begin
       clr = 1
       nclr = 0
       end;COLOR
   'NOCOLOR': begin
       clr = 0
       nclr = 1
       end;COLOR
   'GO': begin 
      numclr = !d.n_colors-1
      set_plot,'ps
      !p.font=0
      !p.charsize=1.0
      device,filename=psfile,/palatino,font_size=10,/color,bits=8
      if (prt eq 1) then device,/portrait
      if (lnd eq 1) then device,/landscape
      if (encp eq 1) then device,encapsulated = 1
      if (noencp eq 1) then device,encapsulated = 0
      device,/inches,xoffset = 1.0,xsize=6.0,yoffset=1.0,ysize=ysz
      if (clr eq 1) then device,color=1
      if (nclr eq 1) then device,color=0
      get_graphics
      device,/close
      set_plot,'x'
      !p.font=-1
      !p.title=''
      !x.title='x'
      !y.title='y'
      !z.title='z'

      widget_control,/destroy, ev.top
      return
      end;GO

endcase

WIDGET_CONTROL, ev.top, get_uvalue=out_text
WIDGET_CONTROL, out_text , set_value=name + ': ' + string(val),/append


return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO get_ps_options
;-------------------------------------------------------------------
@common

wps = widget_base(title='PostScript',/row,/frame)
lc = widget_base(wps,/column,/frame,space=10)
rc = widget_base(wps,/column,space=10)
w1 = widget_label(lc,value='file: ')
w1 = widget_text(lc, /editable,uvalue='FILE')
w1 = widget_label(lc,value='plot title: ')
w1 = widget_text(lc, /editable,uvalue='PTIT')
w1 = widget_label(lc,value='x title: ')
w1 = widget_text(lc, /editable,uvalue='XTIT')
w1 = widget_label(lc,value='y title: ')
w1 = widget_text(lc, /editable,uvalue='YTIT')
w1 = widget_label(lc,value='z title: ')
w1 = widget_text(lc, /editable,uvalue='ZTIT')
t1 = WIDGET_TEXT(rc, xsize=30, ysize=5, /SCROLL, $
      value=[ 'Postscript editor'])
     WIDGET_CONTROL, wps, set_uvalue=t1

xmenu,['encapsulated','non-encasulated'],rc,/frame,/column,/exclusive,$
      space=20,uvalue=['ENCAP','NONENCAP']
xmenu,['portrait','landscape'],rc,/frame,/column,/exclusive,$
      space=20,uvalue=['PORT','LAND']
wysz = widget_slider(rc,title='ysize', uvalue='YSIZE', $
                       maximum = 9, minimum=2)
xmenu,['color','no color'],rc,/frame,/column,/exclusive,$
      space=20,uvalue=['COLOR','NOCOLOR']
w8 = widget_button(rc,value = 'GO',uvalue='GO')


widget_control,/realize,wps,/hourglass
xmanager,'post',wps


return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO view_3d,nfil,nfr
;-------------------------------------------------------------------
@common

device,retain=2,decomposed=0

nfiles=nfil
nfrm = nfiles*nfr

numclr = !d.n_colors-1
loadct,39
vclr=0
vcomp=0
pln='xy'
slcxy=0
slcxz=0
slcyz=0
frm=1
nlev=6
graph_type = 'surf'
graph_type_prev = 'surf'
ntype = 1
file = 'b1all_'
psfile = 'idl.ps'
prt = 1 & lnd = 0
encp = 0 & noencp = 1
clr = 0 & nclr = 1
;numclr=!d.n_colors
!x.title='x'
!y.title='y'
!z.title='z'
ysz = 5
vecscl = 1   ;1 for vector, 0 for scalar

read_file,0
minx = 0 & maxx = nx-1
miny = 0 & maxy = ny-1
minz = 0 & maxz = nz-1

read_np

get_view_widget

end
;-------------------------------------------------------------------







