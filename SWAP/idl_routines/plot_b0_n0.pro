; read 3d coordinate data and plot ambient plasma densities.

f_read_coord,'coord.dat',x,y,z,nx,ny,nz

f_read_3d,'enf.dat',enf
f_read_3d,'b0.dat',b0

;enf=fltarr(nx,ny,nz,/nozero)
;openr,2,efile,/f77_unformatted
;readu,2,enf
;close,2

;shade_surf,reform(alog(enf(*,0,*))),reform(x(*,0,*)),reform(z(*,0,*))

tx=reform(x(*,0,*))
tz=reform(z(*,0,*))
tenf=reform(alog(enf(*,0,*)))
tb0=reform(b0(*,0,*))
triangulate,tx,tz,tr
img1=trigrid(tx,tz,tenf,tr)
img2=trigrid(tx,tz,tb0,tr)
;image_plot,img1
;oplot,x(*,0,*),z(*,0,*)*nz/max(z),psym=3,color=200
;image_plot,img2
;oplot,x(*,0,*),z(*,0,*)*nz/max(z),psym=3,color=200
plot,reform(z(0,0,*)),reform(b0(0,0,*)),psym=1,xtitle='Dist along B (km)', $
     ytitle='B (T)'
plot_io,z(0,0,*),enf(0,0,*),psym=1,xtitle='Dist along B (km)', $
        ytitle='no (m^-3)'

mu=!pi*4*1e-7
mO=16*1.67e-27

valf = b0(0,0,*)/sqrt(mu*mO*enf(0,0,*))
plot,z(0,0,*),valf/1e3,psym=1,xtitle='Dist along B (km)', $
     ytitle='Va (km/s)',title='Alfven velocity along B'

end





