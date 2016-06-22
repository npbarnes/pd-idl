pro plot_vr_sc_b,  time_traj, vmag_all,vr_all, plot_path,scape, output_type, sub_path, yrange,xrange,plot_label,poso, ilab, it_str, it_vr,ytitle

if scape eq 'portrait' then begin 
  xin=8.5
  yin=11.0
  xdim=5.9
  ydim=3.5
  x_off=1.53;1.13
  y_off=5;8.2
  color_code=0
  file=plot_path+sub_path+'vr_sc.ps'
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
  file=plot_path+sub_path+'vr_sc.ps'
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

ytitle='Speed [km/s]'
MAKEDOTS
plot, time_traj.tyrr, vmag_all, xthick=thick, ythick=thick, pos=pos, /normal,psym=-8, symsize=symsize, ystyle=1,yrange=yrange, xstyle=1, xrange=xrange,/noerase, ytitle=ytitle,yminor=2, ycharsize=1.2,thick=6,xcharsize=1.2
oplot, time_traj.tyrr, vr_all, color=252, thick=6,psym=-8 ,symsize=symsize
xyouts, pos(2)-.2, pos(3)-.03, 'NH Speed', color=0,/normal, charsize=0.9
xyouts, pos(2)-.2, pos(3)-.06, 'NH Speed Along Radial', color=252,/normal, charsize=0.9
xyouts, pos(0), pos(3)+.05, plot_label,color=0,/normal
oplot, [it_str.tyrr, it_str.tyrr], [!y.crange], color=90
oplot, [!x.crange], [it_vr,it_vr], color=90
xyouts, pos(0)+.02, pos(3)-.03, 'Vr. at '+ilab+': '+strcompress(string(it_vr, format='(F30.5)'),/remove_all)+' [km/s]',/normal, color=90
if output_type eq 'P' or output_type eq 'p' then device, /close
end
