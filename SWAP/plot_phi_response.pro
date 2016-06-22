pro plot_phi_response, plot_path, output_type, scape, wphi
;__________________________________
;____Plot Phi Response
;__________________________________
file=plot_path+'phi_response.ps'
  if scape eq 'portrait' then begin 
  !x.ticklen=.025
  xdim=5
  ydim=2.5
  x_off=1.5
  y_off=3.1
  xin=8.5;6.692913
  yin=11;8.976378
  file_name=file
  color_code=0
  PLOT_TYPE_POSITION_C, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type,file, color_code 
endif else begin 
  !x.ticklen=.05;.05
  xdim=6
  ydim=2.5
  x_off=1.5
  y_off=1.5
  xin=11.0
  yin=8.0
  file_name=file
  color_code=0
  plot_type_position_cland, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type, FILE,color_code
endelse
plot, wphi.phi, wphi.w, xthick=3, ythick=3, thick=3, /nodata, pos=[xmn, ymn, xmx, ymx], xtitle='Phi Angle [deg]', ytitle='Phi Response', charsize=1.3
MAKEDOTS
oplot, wphi.phi, wphi.w,psym=-8, thick=3, color=252

end
