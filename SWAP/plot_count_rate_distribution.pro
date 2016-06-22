pro plot_count_rate_distribution, plot_path,output_type, scape, sim_count_spin, sim_count, einst
if scape eq 'portrait' then begin 
  xin=8.5
  yin=11.0
  xdim=5.8
  ydim=3.5
  x_off=1.5;1.13
  y_off=5;8.2
  color_code=0
  file=plot_path+'count_rate_energy.ps'
  PLOT_TYPE_POSITION_C, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type, FILE,color_code
  !x.ticklen=-.03
  !y.ticklen=-.02
endif else begin 
  yin=8.5
  xin=11.0
  xdim=8.3
  ydim=3.7
  x_off=1.5;1.5
  y_off=3
  color_code=0
  file=plot_path+'count_rate_energy.ps'
  PLOT_TYPE_POSITION_CLAND, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type, FILE,color_code
  !x.ticklen=-.03
  !y.ticklen=-.02
endelse
norm_wid=xmx-xmn
inch_wid=xdim
poso=[xmn, ymn, xmx, ymx]
loadct, 39,/silent
MAKEDOTS
iz=where(sim_count_spin gt 0)
plot, einst(iz), sim_count_spin(iz),/noerase, /xlog,xrange=[100,2000], yrange=[1,1E6],/ylog, yminor=0,yticklen=-0.025, ystyle=1, xstyle=1, color=0, psym=8,symsize=.7,pos=poso,ytitle='Simulated!CCoincidence Rate [Hz]', xtitle='Energy/Charge [eV/q]',ycharsize=1.2, xcharsize=1.2,xthick=3, ythick=3
iz=where(sim_count gt 0)
oplot, einst(iz), sim_count(iz), color=252,psym=8,symsize=.7
xyouts, poso(0)+.02, poso(3) -.03 , 'Spinning Example',/normal
xyouts, poso(0)+.02, poso(3) -.06 , 'Same Example, but with Theta and Phi fixed at 3.6 deg', color=252,/normal
if strlowcase(output_type) eq 'p' then device, /close

end
