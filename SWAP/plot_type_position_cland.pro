pro plot_type_position_cland, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type, FILE,color_code
;________________________________________________
;_______Program name: plot_type_position.pro
;_______Author: Heather Alison Elliott and Peggy Sloan
;_______Last Revised:09/22/98
;_______Purpose: Set colors, position plot, and set plot type
;________________________________________________
;set_plot, 'X'
;device, true_color=24
;device, decompose=0
bgcolor=255
;______________set up plot features
pi=4*atan(1)
;________setting color tables
;red=[0,1,1,0,0,1]
;green=[0,1,0,1,0,1]
;blue=[0,1,0,0,1,0]
;TVLCT, 255*red,255*green,255*blue
;_______yellow-5
;_______green-3
;_______white-1
;_______red-2
;_______blue-4
;_______black-6  
;output_type=' '
;print, '   For gif choose W and then use make_gif.'
;read, output_type, prompt='Enter plot type (P or W):'
if (output_type eq 'P') then begin
   print, ' '
   print, 'Generating postscript file'
endif else if (output_type eq 'G' ) then begin
   print, ' '
   print, 'Generating gif file'
endif
;______________initialize plot size 
;yin=11.0
;xin=8.5
ywin=fix(yin*75.0)
xwin=xin*75.0
ximg=xdim/xin
yimg=ydim/yin
mx=ximg*xwin
my=yimg*ywin
xspace=1.0/xin
;print, ' initialize plot size'
;______________________
;---------------------------------------
;  set up user-defined symbol - a circle
;---------------------------------------
if(output_type eq 'G' or output_type eq 'g') then begin
  rad = 1.00 
  thk = 1
endif else if(output_type eq 'W' or output_type eq 'w') then begin
  rad = 1.30 
  thk = 2
endif else begin
  rad = 1.20
  thk = 4
endelse
circle = fltarr(2,11)
deg = 2.0 * pi / 10.0
for ii = 0,10 do begin
   circle(0,ii) = rad * cos(deg*ii)
   circle(1,ii) = rad * sin(deg*ii)
endfor
      ;---------------------------
      ; initialize plot parameters
      ;---------------------------
;_________I think xmn and ymn are offsets
      xmn=x_off/xin
      xmx = xmn + ximg
      ;ymn = 2.5/yin
       ymn=y_off/yin
      ymx = ymn + yimg
      yctr = (ymn + ymx) /2.0
      xd = 0.1 / xin
      yd = 0.1 / yin
      lg_sz = 1.1
      mid_sz = 0.9
      sm_sz = 0.7
      if (output_type eq 'G' or output_type eq 'g') then begin
         lg_sz = lg_sz*0.8
         mid_sz = mid_sz*0.8
         sm_sz = sm_sz*0.8
      endif
      ;-----------------------
      ; set up plotting device
      ;-----------------------
;       print, 'set up plot device'
;       print, output_type, 'out'
      if (output_type eq 'W' or output_type eq 'w') then begin
         ;print, 'w type'
         if (!version.os eq 'MacOS') then begin
            set_plot, 'Mac'
            ul_x = 25
            ul_y = 25
         endif else begin
            set_plot, 'X'
            device, get_screen_size = scrn_size
            ul_x = 25
            ul_y = (scrn_size(1) - ywin) - 25
         endelse
         loadct,39,/silent
         !p.background=bgcolor
         !p.color=color_code
         window,/free, title = '',  $
                 xpos = ul_x, ypos = ul_y, $
                 xsize = xwin, ysize = ywin, $
                 retain = 2,color=bgcolor
         erase
      endif else if (output_type eq 'G' or output_type eq 'g') then begin
         if (!version.os eq 'MacOS') then begin
            set_plot, 'Mac'
            ul_x = 25
            ul_y = 25
         endif else begin
            set_plot, 'X'
            device, get_screen_size = scrn_size
            ul_x = 25
            ul_y = (scrn_size(1) - ywin) - 25
         endelse
         ;loadct,39,/silent
         !p.background=bgcolor
         !p.color=color_code
         window, /free, title = 'TIDE Moments',  $
                 xpos = ul_x, ypos = ul_y, $
                 xsize = xwin, ysize = ywin, $
                 retain = 2, color=bgcolor
         erase
      endif else if (output_type eq 'P' or output_type eq 'p') then begin
         ;loadct,39,/silent
         !p.background=bgcolor
         !p.color=color_code
         set_plot,'ps'
          ;FILE=''
          ;PRINT,'Enter PostScript Filename:'
          ;READ,FILE
        print,' Starting PostScript File : '+FILE
;________easy way to move plot on page just change xoff and yoff
        DEVICE,FILE=FILE,/INC,XSI=xin ,YSI=yin, XOFF=0.0, YOFF=11.0,$
          /HELVETICA, /color, /landscape, bits=24
        !P.FONT=0
      endif
return
end
