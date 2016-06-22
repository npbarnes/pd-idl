 pro plot_color_4_eb_ea_final, ecen, i, x, xmn, xmx, ymn, ymx, arr_ebea, new_ebea, ct_label, mmin,mmax
posa=[xmn, ymn, xmx, ymx]
yshift=0.08
xshift=0.15
posb=posa
posb(0)=posb(0)+xmx-xmn+xshift
posb(2)=posb(2)+xmx-xmn+xshift

posa1=posa
posa1(1)=posa1(1)-(ymx-ymn)-yshift
posa1(3)=posa1(3)-(ymx-ymn)-yshift
posb1=posa1
posb1(0)=posb1(0)+xmx-xmn+xshift
posb1(2)=posb1(2)+xmx-xmn+xshift
!y.tickinterval=0
!x.tickinterval=5
!x.thick=3
!y.thick=3
!x.ticklen=-0.03
!y.ticklen=-0.03
!y.minor=0
!y.charsize=1
  
del=(new_ebea(i,1)-new_ebea(i,0))
yrange=[min(new_ebea(i,*))-del/2.0,max(new_ebea(i,*))+del/2.0]
  
xrange=[min(x)-1.0/2.0, max(x)+1.0/2.0]
plot, xrange,yrange,  $
ytitle = 'E!DB!N/(kV!DESA!N)', xtitle = 'Azimuth [deg]', $
  position = posa, xrange=xrange, yrange=[yrange(0),yrange(1)], xstyle=1, ystyle=1
img=transpose(reverse(arr_ebea(*,*,i),2))
loadct, 39,/silent
tvimage,(smooth(bytscl(alog10(img),/nan,min = mmin, max = mmax),1)), $
    /noin,/overplot
loadct, 39,/silent
  contour, transpose(arr_ebea(*,*,i)), x, (new_ebea(i,*)),ytitle='',xtitle='',ystyle=1,xrange=xrange, yrange=yrange, xstyle=1,/noerase, pos=posa,/normal, xticklen=-.03, yticklen=-.03, color=52
xyouts, posa(0)+.12, posa(3)+.02, 'SIMION Model Results       Beam Energy = ' + strcompress(string(ecen(i), format='(F20.3)'),/remove_all)+' [eV]',/normal
!y.ticklen=-.15
!y.minor=2
!y.charsize=.8

!y.tickinterval=0
!x.tickinterval=5
!x.thick=3
!y.thick=3
!x.ticklen=-0.03
!y.ticklen=-0.03
!y.minor=0
!y.charsize=1
  
del=(new_ebea(i,1)-new_ebea(i,0))
yrange=[min(new_ebea(i,*))-del/2.0,max(new_ebea(i,*))+del/2.0]
  
xrange=-[max(x)+1.0/2.0,min(x)-1.0/2.0]
plot, xrange,yrange,  $
  ytitle = 'E!DB!N/(kV!DESA!N)', xtitle = '!9q!X [deg]', $
  position = posb, xrange=xrange, yrange=[yrange(0),yrange(1)], xstyle=1, ystyle=1,/noerase
img=transpose((arr_ebea(*,*,i)))
 
loadct, 39,/silent
tvimage,(smooth(bytscl(alog10(img),/nan,min = mmin, max = mmax),1)), $
    /noin,/overplot
loadct, 39,/silent
  contour, transpose(arr_ebea(*,*,i)), -x, (new_ebea(i,*)),ytitle='',xtitle='',ystyle=1,xrange=xrange, yrange=yrange, xstyle=1,/noerase, pos=posb,/normal, xticklen=-.03, yticklen=-.03, color=52
!y.ticklen=-.15
!y.minor=2
!y.charsize=.8
!y.tickinterval=0
!x.tickinterval=0
color_bar, [mmin, mmax], ct_label,[posb(2)+.04,posb(3)-.13,posb(2)+.06,posb(3)-.04]

!x.tickinterval=0
!y.tickinterval=5.
!x.thick=3
!y.thick=3
!x.ticklen=-0.03
!y.ticklen=-0.03
!y.minor=0
!y.charsize=1
del=(new_ebea(i,1)-new_ebea(i,0))
yrange=[min(new_ebea(i,*))-del/2.0,max(new_ebea(i,*))+del/2.0]

xrange=[min(x)-1.0/2.0, max(x)+1.0/2.0]
plot, yrange,xrange,  $
  xtitle = 'E!DB!N/(kV!DESA!N)', ytitle = 'Azimuth [deg]',/noerase, $
  position = posa1, yrange=xrange, xrange=[yrange(0),yrange(1)], xstyle=1, ystyle=1
img=(reverse(arr_ebea(*,*,i),2))
 
loadct, 39,/silent
tvimage,(smooth(bytscl(alog10(img),/nan,min = mmin, max = mmax),1)), $
    /noin,/overplot
loadct, 39,/silent
  contour, arr_ebea(*,*,i), (new_ebea(i,*)),x,ytitle='',xtitle='',ystyle=1,yrange=xrange, xrange=yrange, xstyle=1,/noerase, pos=posa1,/normal, xticklen=-.03, yticklen=-.03, color=52
;stop
!x.tickinterval=0
!y.tickinterval=5.
!x.thick=3
!y.thick=3
!x.ticklen=-0.03
!y.ticklen=-0.03
!y.minor=0
!y.charsize=1

del=(new_ebea(i,1)-new_ebea(i,0))
yrange=[min(new_ebea(i,*))-del/2.0,max(new_ebea(i,*))+del/2.0]
xrange=-[max(x)+1.0/2.0,min(x)-1.0/2.0]
plot, yrange,xrange,  $
  xtitle = 'E!DB!N/(kV!DESA!N)', ytitle = '!9q!X [deg]', $
  position = posb1, yrange=xrange, xrange=[yrange(0),yrange(1)], xstyle=1, ystyle=1,/noerase

img=((arr_ebea(*,*,i)))
loadct, 39,/silent
tvimage,(smooth(bytscl(alog10(img),/nan,min = mmin, max = mmax),1)), $
    /noin,/overplot
loadct, 39,/silent
  contour, arr_ebea(*,*,i), (new_ebea(i,*)),-x,ytitle='',xtitle='',ystyle=1,yrange=xrange, xrange=yrange, xstyle=1,/noerase, pos=posb1,/normal, xticklen=-.03, yticklen=-.03, color=52
!y.ticklen=-.15
!y.minor=2
!y.charsize=.8
!y.tickinterval=0
!x.tickinterval=0
color_bar, [mmin, mmax], ct_label,[posb1(2)+.04,posb1(3)-.13,posb1(2)+.06,posb1(3)-.04]

end
