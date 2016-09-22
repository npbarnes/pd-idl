;-------------------------------------------------------------------
pro read_file
;-------------------------------------------------------------------
@common

close,1
openr,1,file(0),/f77_unformatted

frame=0

nt=0
nout=0
nx=0
ny=0
nz=0

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

temparr = fltarr(nx,ny,nz,3,/nozero)

for i = 0,frm-1 do begin
   readu,1,frame
   print,'image #.....',frame
   readu,1,temparr
endfor

arr = temparr

close,1

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_surface
;-------------------------------------------------------------------
@common

case pln of

   'xy': surface,reform(arr(*,*,slcxy-1,vcomp)),charsize=1.5
    
   'xz': surface,reform(arr(*,slcxz-1,*,vcomp)),charsize=1.5

   'yz': surface,reform(arr(slcyz-1,*,*,vcomp)),charsize=1.5

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_contour
;-------------------------------------------------------------------
@common

case pln of

   'xy': contour,reform(arr(*,*,slcxy-1,vcomp)),/c_annotation, $
                 charsize=1.5,nlevels=6
    
   'xz': contour,reform(arr(*,slcxz-1,*,vcomp)),/c_annotation, $
                 charsize=1.5,nlevels=6

   'yz': contour,reform(arr(slcyz-1,*,*,vcomp)),/c_annotation, $
                 charsize=1.5,nlevels=6

endcase

return
end
;-------------------------------------------------------------------

;-------------------------------------------------------------------
PRO VELOVEC,U,V,X,Y, Missing = Missing, Length = length, Dots = dots,  $
        Color=color, _EXTRA = extra
;-------------------------------------------------------------------

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

      vxy = reform(arr(*,*,slcxy-1,0))
      vyx = reform(arr(*,*,slcxy-1,1))

      if not((abs(max(vxy)) eq 0) and (abs(max(vyx)) eq 0)) then begin
         velovec,vxy,vyx,xtitle='x',ytitle = 'y',length=1.0, $
                 title='max ='+string(max(sqrt(vxy^2 + vyx^2))), $
                 charsize=1.5      
      endif
      end;xy

   'xz': begin
 
      vxz = reform(arr(*,slcxz-1,*,0))
      vzx = reform(arr(*,slcxz-1,*,2))

      if not((abs(max(vxz)) eq 0) and (abs(max(vzx)) eq 0)) then begin
         velovec,vxz,vzx,xtitle='x',ytitle = 'z',length=1.0, $
                 title='max ='+string(max(sqrt(vxz^2 + vzx^2))), $
                 charsize=1.5
      endif
      end;'xz'

   'yz': begin

      vyz = reform(arr(slcyz-1,*,*,1))
      vzy = reform(arr(slcyz-1,*,*,2))

      if not((abs(max(vyz)) eq 0) and (abs(max(vzy)) eq 0)) then begin
         velovec,vyz,vzy,xtitle='y',ytitle='z',length=1.0, $
                 title='max ='+string(max(sqrt(vyz^2 + vzy^2))), $
                 charsize=1.5
      endif
      end;yz

endcase

return
end
;------------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_movie_frames
;-------------------------------------------------------------------
@common

close,1
openr,1,file(0),/f77_unformatted

nt=0
nout=0
nx=0
ny=0
nz=0

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

case pln of
   'xy': anmt_arr = fltarr(nx,ny,nfrm)
   'xz': anmt_arr = fltarr(nx,nz,nfrm)
   'yz': anmt_arr = fltarr(ny,nz,nfrm)
endcase

temparr = fltarr(nx,ny,nz,3,/nozero)

i=0
while (not(eof(1))) do begin
   readu,1,frame
   readu,1,temparr
   case pln of
      'xy': anmt_arr(*,*,i) = reform(temparr(*,*,slcxy-1,vcomp))
      'xz': anmt_arr(*,*,i) = reform(temparr(*,slcxz-1,*,vcomp))
      'yz': anmt_arr(*,*,i) = reform(temparr(slcyz-1,*,*,vcomp))
   endcase
   i = i + 1
endwhile

close,1

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
      arr2=congrid(arr1,100,100,nfrm)
      xinteranimate,set=[100,100,nfrm]
      for i=0,nfrm-1 do begin
         xinteranimate,image=arr2(*,*,i),frame=i
      endfor
      xinteranimate
      end;xy
    
   'xz': begin
      get_movie_frames
      arr1=bytscl(anmt_arr)
      arr2=congrid(arr1,100,200,nfrm)
      xinteranimate,set=[100,200,nfrm]
      for i=0,nfrm-1 do begin
         xinteranimate,image=arr2(*,*,i),frame=i
      endfor
      xinteranimate
      end;xz
        
   'yz': begin
      get_movie_frames
      arr1=bytscl(anmt_arr)
      arr2=congrid(arr1,50,200,nfrm)
      xinteranimate,set=[50,200,nfrm]
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
      read_file
      make_surface
      end;surf

   'cont': begin
      read_file
      make_contour
      end;cont

   'anim': begin
      make_animate
;      xinteranimate,/close
      end;anim

   'vect': begin
      read_file
      make_vector
      end;vect

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
      pln=plane
      widget_control,wxy,sensitive=1
      widget_control,wxz,sensitive=0
      widget_control,wyz,sensitive=0
      if not(ev.select eq 0) then get_graphics
      end;XY

   'XZ': begin
      widget_control, ev.id, get_value = plane
      pln=plane
      widget_control,wxy,sensitive=0
      widget_control,wxz,sensitive=1
      widget_control,wyz,sensitive=0
      if not(ev.select eq 0) then get_graphics
      end;XZ

   'YZ': begin
      widget_control, ev.id, get_value = plane
      pln=plane
      widget_control,wxy,sensitive=0
      widget_control,wxz,sensitive=0
      widget_control,wyz,sensitive=1
      if not(ev.select eq 0) then get_graphics
      end;YZ

   'SLICEXY': begin
      widget_control, ev.id, get_value = slice
      slcxy=slice
      get_graphics
      end;SLICE

   'SLICEXZ': begin
      widget_control, ev.id, get_value = slice
      slcxz=slice
      get_graphics
      end;SLICE

   'SLICEYZ': begin
      widget_control, ev.id, get_value = slice
      slcyz=slice
      get_graphics
      end;SLICE

   'FRAME': begin
      widget_control, ev.id, get_value = frame
      frm=frame
      get_graphics
      end;FRAME

   'SURFACE': begin
      graph_type = 'surf'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'surf'
      end;SURFACE

   'CONTOUR': begin
      graph_type = 'cont'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'cont'
      end;CONTOUR

   'VECTOR': begin
      graph_type = 'vect'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'vect'
      end;ANIMATE

   'ANIMATE': begin 
      graph_type = 'anim'
      get_graphics
      graph_type = graph_type_prev
      end;ANIMATE
  
   'VIEWXY': partmovie_xy,'npall.dat',nfrm
   'VIEWXZ': partmovie_xz,'npall.dat',nfrm
   'VIEWYZ': partmovie_yz,'npall.dat',nfrm
  
 
   'B1': begin
      file = 'b1all.dat'
      if not(ev.select eq 0) then get_graphics
      end;B1

   'UF': begin 
      file = 'ufall.dat'
      if not(ev.select eq 0) then get_graphics
      end;UF

   'AJ': begin 
      file = 'ugraduall.dat'
      if not(ev.select eq 0) then get_graphics
      end;AJ
      
   'E' : begin
      file = 'Eall.dat'
      if not(ev.select eq 0) then get_graphics
      end;E      

   'DONE': widget_control,/destroy, ev.top

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_view_widget
;-------------------------------------------------------------------
@common

base = widget_base(title='view_3d_vec',/row)
lc = widget_base(base,/column,space=10)
mc = widget_base(base,/frame,/column)
rc = widget_base(base,/column,space=50)
rct = widget_base(rc,/column,space=5)
rcm = widget_base(rc,/frame,/column,space=5)
rcb = widget_base(rc,/frame,/column,space=5)


xmenu,['x','y','z'],lc,/column,/exclusive,$
      title='vector component',space=50,uvalue=['X','Y','Z']
xmenu,['xy','xz','yz'],lc,/column,/exclusive,$
      title='plane',space=50,uvalue=['XY','XZ','YZ']
xmenu,['surface','contour','vector'],rct,/column, $
      /exclusive,title='graphics options',space=50,$
      uvalue=['SURFACE','CONTOUR','VECTOR']

w1 = widget_label(rcm,value = 'animate fluid')
w1 = widget_button(rcm,value = 'go',uvalue = 'ANIMATE')

w1 = widget_label(rcb,value = 'animate cloud')
w1 = widget_button(rcb,value = 'xy',uvalue='VIEWXY')
w1 = widget_button(rcb,value = 'xz',uvalue='VIEWXZ')
w1 = widget_button(rcb,value = 'yz',uvalue='VIEWYZ')



wxy = widget_slider(lc,title='z slice',uvalue='SLICEXY', $
                         maximum=nz, minimum = 1)
widget_control,wxy,sensitive=0
wxz = widget_slider(lc,title='y slice',uvalue='SLICEXZ', $
                         maximum=ny, minimum=1)
widget_control,wxz,sensitive=0
wyz = widget_slider(lc,title='x slice',uvalue='SLICEYZ', $
                         maximum=nx, minimum=1)
widget_control,wyz,sensitive=0
wfrm = widget_slider(lc,title='frame #',uvalue='FRAME', $
                       maximum=nfrm, minimum = 1)

xmenu,['b1','uf','aj','E'],mc,/row,/exclusive, $
      space=20,uvalue=['B1','UF','AJ','E']

draw = widget_draw(mc,xsize=700,ysize=600)

w1 = widget_button(mc,value = 'done',uvalue='DONE')

widget_control,/realize,base,/hourglass
widget_control,get_value = window,draw
wset,window

xmanager,'view',base

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro file_event,ev
;-------------------------------------------------------------------
@common

widget_control,ev.id,get_value = file
print,file

close,1
openr,1,file(0),/f77_unformatted

nt=0
nout=0
nx=0
ny=0
nz=0

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

close,1
print,nt,nout,nx,ny,nz

widget_control,/destroy,ev.top
close,1

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_file
;-------------------------------------------------------------------
@common
 
;w1 = widget_base(title='get_file',/row,/frame)
;w2 = widget_label(w1,value='Enter file: ')
;w2 = widget_text(w1, /editable,uvalue='FILE')

;widget_control,/realize,w1
;xmanager,'file',w1

close,1
openr,1,file(0),/f77_unformatted

frame=0

nt=0
nout=0
nx=0
ny=0
nz=0

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,ny
readu,1,nz

temparr = fltarr(nx,ny,nz,3,/nozero)

nfrm=0
while not(eof(1)) do begin
   readu,1,frame
;   print,'image #.....',frame
   readu,1,temparr
   nfrm=nfrm+1
endwhile

close,1

return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
; PROGRAM view_3d_vec
;-------------------------------------------------------------------
@common

vcomp=0
pln='xy'
slcxy=1
slcxz=1
slcyz=1
frm=1

graph_type = 'surf'
graph_type_prev = 'surf'
file = 'b1all.dat'

get_file
get_view_widget

end
;-------------------------------------------------------------------






