pro plot_eff,  time, eff,  plot_path,scape, output_type, sub_path, yrange,xrange,plot_label,poso, ilab, it_str, eff_fact,ytitle

if scape eq 'portrait' then begin 
  xin=8.5
  yin=11.0
  xdim=5.8
  ydim=3.5
  x_off=1.53;1.13
  y_off=5;8.2
  color_code=0
  file=plot_path+sub_path+'swap_detect_eff.ps'
  PLOT_TYPE_POSITION_C, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type, FILE,color_code
  !x.ticklen=-.03
  !y.ticklen=-.02
endif else begin 
  yin=8.5
  xin=11.0
  xdim=8.3
  ydim=3.7
  x_off=1.36;1.5
  y_off=3
  color_code=0
  file=plot_path+sub_path+'swap_detect_eff.ps'
  PLOT_TYPE_POSITION_CLAND, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type, FILE,color_code
  !x.ticklen=-.03
  !y.ticklen=-.02
endelse
norm_wid=xmx-xmn
inch_wid=xdim
poso=[xmn, ymn, xmx, ymx]
thick=4
symsize=0.4
pos=poso
ht = pos(3)-pos(1)
gap=.012
  
blank=strarr(60)
blank(*)=' '
 

MAKEDOTS
plot, time, eff, xthick=thick, ythick=thick, pos=pos, /normal,psym=-8, symsize=symsize, ystyle=1,yrange=yrange, xstyle=1, xrange=xrange,/noerase, ytitle=ytitle, ycharsize=1.2,thick=6,xcharsize=1.2,/nodata,/ylog
for i=0, n_elements(time)-2 do begin
  t1=time(i)
  t2=time(i+1)
  if t1 lt !x.crange(0) then t1= !x.crange(0) 
  if t2 gt !x.crange(1) then t2= !x.crange(1) 
  y1=eff(i)
  y2=eff(i+1)
  plots, [t1,t2], [eff(i), eff(i+1)], color=252, thick=6,psym=-8
  ;print, t1, t2, eff(i), eff(i)
endfor
xyouts, pos(0), pos(3)+.05, plot_label,color=0,/normal
oplot, [it_str.tyrr, it_str.tyrr], 10^[!y.crange]
oplot, [!x.crange], [eff_fact,eff_fact]
xyouts, pos(0)+.02, pos(3)-.03, 'Detect Eff. at '+ilab+': '+strcompress(string(eff_fact, format='(F30.5)'),/remove_all),/normal
if output_type eq 'P' or output_type eq 'p' then device, /close
end
