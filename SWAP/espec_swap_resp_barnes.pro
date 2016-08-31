;---------------------------------------------------------------------------;
pro img_cont_ylog, a, x, y,dunit,xpos1,xpos2,ypos1,ypos2,evst_r,WINDOW_SCALE = window_scale, $
                ASPECT = aspect, INTERP = interp,POSTSCRIPT=postscript
;+
; NAME:
;	IMAGE_CONT
; PURPOSE:
;	Overlay an image and a contour plot.
; CATEGORY:
;	General graphics.
; CALLING SEQUENCE:
;	IMAGE_CONT, A
; INPUTS:
;	A = 2 dimensional array to display.
; KEYWORD PARAMETERS:
;	/WINDOW_SCALE = set to scale the window size to the image size,
;		otherwise the image size is scaled to the window size.
;		Ignored when outputting to devices with scalable pixels.
;	/ASPECT = set to retain image's aspect ratio.  Assumes square
;		pixels.  If /WINDOW_SCALE is set, the aspect ratio is
;		retained.
;	/INTERP = set to bi-linear interpolate if image is resampled.
; OUTPUTS:
;	No explicit outputs.
; COMMON BLOCKS:
;	none.
; SIDE EFFECTS:
;	The currently selected display is affected.
; RESTRICTIONS:
;	None that are obvious.
; PROCEDURE:
;	If the device has scalable pixels then the image is written over
;	the plot window.
; MODIFICATION HISTORY:
;	DMS, May, 1988.
;-

nclr = !d.n_colors

if keyword_set(postscript) then begin
set_plot,'ps
device,filename='img.eps
device,/encapsulated
!p.font=0
;device,/hevletica
device,bits=8
device,/color
;loadct,33
device,/inches,xsize=6.5,ysize=6.5,xoffset =1.0, yoffset=1.0
!p.thick=2.0
!x.thick=2.0
!y.thick=2.0
endif

maxa=max(a)
mina=min(a)

sz=size(a)

!p.charsize=1.2
!y.margin = [10,10]
!x.margin = [10,10]

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour

contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,$
            ytitle='Energy/q (eV/q)',$
            xtitle='x (Rp)',$
            xrange=[max(x),min(x)],yrange=[24,max(y)],/ylog,$
        position=[xpos1,ypos1,xpos2,ypos2]
;        xticks=1,$
;        xtickname = [' ',' '],$
;        xtickv=[min(x),max(x)],xminor=1
;            xrange=[min(x),max(x)],yrange=[10,5000],/ylog
;            title='t = '+strmid(strtrim(string(0.5*nfrm*100.),2),0,6) + ' (s)',$
;            /isotropic
;            xrange=[min(x),800],yrange=[min(y),max(y)+1]
p1 = !P & x1 = !X & y1= !Y

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
  ;print,max(a)
  tv,a,px(0),py(0),xsize = swx, ysize = swy, /device


endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tv,a,px(0),py(0)	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
        ;print,max(a)
	tv,poly_2d((a),$	;Have to resample image
                   [[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
                   keyword_set(interp),swx,swy), $
           px(0),py(0)

	endelse			;window_scale
  endelse			;scalable pixels

     axis,xaxis=2,xstyle=1
     axis,yaxis=2,ystyle=1
;     axis,yaxis=1,ystyle=1
;     axis,xaxis=1,xstyle=1
;;     yticks=2,$
;;     ytickv=[min(y),(min(y)+max(y))/2,max(y)], $
;;     ytickname=[string(min(a)),string((min(a)+max(a))/2),string(max(a))], $
;;     ytitle='E!dz!n (mV/m)',$
;;     yrange = [mina*(2.3e-25/1.6e-19)*1e6,maxa*(2.3e-25/1.6e-19*1e6)]
;     ytitle='(u!de!n)!dz!n (km/s)',$
;     yrange = [mina,maxa]


mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
if !d.name eq 'PS' then colors = mx - colors ;invert line colors for pstscrp
m = max(a)
;contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1
;contour,a,/noerase,/nodata,/yst,$	;Do the contour
;	  pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
;          c_color =  colors, charsize=1.5, $
;          xtitle='x (km)', ytitle='y (km)', $
;          levels=[0.02*m,0.04*m,0.06*m,0.2*m,0.4*m,0.5*m,0.7*m,0.9*m]
;;	  xrange=[min(ax),max(ax)], xstyle=1
;dx = !x.crange(1) - !x.crange(0)
;dy = !y.crange(1) - !y.crange(0)
;xyouts,0.05*dx+!x.crange(0),0.85*dy+!y.crange(0),tit,/data
;!p.multi=[0,1,2]

!P = p1 & !X = x1 & !Y = y1

;contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,/noerase

contour,[[0,0],[1,1]],/nodata, xstyle=1, ystyle = 1,$
            ytitle='Energy/q (eV/q)',$
            xtitle='x (Rp)',$
            xrange=[max(x),min(x)],yrange=[24,max(y)],/ylog,$
        position=[xpos1,ypos1,xpos2,ypos2],/noerase

cgColorbar, position=[!y.window(0),!x.window(1)+0.05,!y.window(1),$
                      !x.window(1)+0.08],$
           ncolors=255,$
           maxrange=[exp(max(evst_r))],minrange=[1],$
           title='Counts/s',$
           charsize=1.4,/xlog,xtickformat='log_label'

if keyword_set(postscript) then begin
device,/close
set_plot,'x'
!p.font=-1
endif

return
end
;---------------------------------------------------------------------------;


;----------------------------------------------------------------
PRO read_part,file,nfrm,Ni_max,xp
;----------------------------------------------------------------

; read the nfrm'th frame of the vector particle file 'file' with Ni_max particles.
; store result in xp

;Ni_max=long(0)
;nt=0l
;ntout=0l
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,nt,ntout,Ni_max




xp=fltarr(Ni_max,3,/nozero)


readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end

;----------------------------------------------------------------


;----------------------------------------------------------------
PRO read_part_scalar,file,nfrm,Ni_max,xp
;----------------------------------------------------------------
; read the nfrm'th frame of the scalar particle file 'file' with Ni_max particles.
; store result in xp

;Ni_max=long(0)
;nt=0
;ntout=0
frm=0l

file = file+'.dat'
print,' reading...',file
openr,1,file,/f77_unformatted
;readu,1,nt
;readu,1,ntout
;readu,1,Ni_max
;print,Ni_max

xp=fltarr(Ni_max,/nozero)

readu,1,frm
print,'  image #.....',frm
readu,1,xp
frmcnt = 1

while (frmcnt lt nfrm) do begin

   readu,1,frm
   print,'  image #.....',frm
   readu,1,xp
   frmcnt = frmcnt + 1

endwhile

close,1

return
end
;----------------------------------------------------------------





;------------------------------------------------------------
function get_swap_resp,vc,theta,phi,w,s4,eff
;------------------------------------------------------------
  aeff=0.033
  aeff=aeff/0.0882d
  dect_eff=eff
  aeff=dect_eff*aeff ;cm^2
  aeff = aeff/1e10   ;km^2
  mp=1.67262158D-27
  kb=1.380658D-23
  e=1.6E-19 
  ee=.5*mp/e*(vc*1000.)^2 ;vc in km/s...ee in eV

  wp=interpol(w.w, w.phi,phi)
  junk=min(abs(ee-s4.ecen),iee)

  junk=min(abs(theta-(-s4.x)),ix)

  er1=s4.y(iee,*)*ee

  tt=s4.arr(*,ix,iee)

  ttnew=interpol(tt,er1, ee)
  if (ee lt min(er1) and ee gt max(er1)) then ttnew = 0.0d

  res = wp*ttnew*aeff

return,res
end
;------------------------------------------------------------



;------------------------------------------------------------
function get_dect_eff, eff
;------------------------------------------------------------
common NH_traj_info,traj_data,time_traj,it_str,file_path,traj_met


;________________________________________________
;____Loading in the Dectector Efficiency vs Time
;________________________________________________
; __________________________
; ____Set Types
; _________________________
; _____byte 1
; _____int 2
; _____long 3
; _____float 4
; _____double 5
; _____complex 6
; _____string 7 
fieldtypes=dblarr(3)
csvfile=file_path+'er_table_lookup_3_2015c.csv'
fieldtypes(*)=7
fieldtypes(0)=5
fieldtypes(2)=5
fin=READ_CSV_FIELDTYPES(strcompress(csvfile,/remove_all), fieldtypes)
convert_date_time_to_year_arr_spice, fin.utc, junktime,ef
met_mid=double(interpol(traj_met, time_traj.tyrr, it_str.tyrr))

;________________________________________________________
;____Interpolate to find the correct detector efficiency
;____   a specific time.  A plot is give on slide 6.
;____   Note that you have to calibraiton data taken 
;____   with the current operational voltage to 
;____   calculate the current efficiency. The table
;____   has been set up to handle this issue. 
;_______________________________________________________
met_junk=min(abs(met_mid-fin.met),imet)
if fin.met(imet) gt met_mid then iarr=[imet-1, imet] else iarr=[imet, imet+1]
edetector=fin.eff(iarr(0))+(fin.eff(iarr(1))-fin.eff(iarr(0)))/(fin.met(iarr(1))-fin.met(iarr(0)))*(double(met_mid)-fin.met(iarr(0)))
if fin.met(imet) eq met_mid then edetector = fin.eff(imet)
;_________________________________________________________________


return,edetector
end
;------------------------------------------------------------


;------------------------------------------------------------
function get_NH_vr
;------------------------------------------------------------
common NH_traj_info,traj_data,time_traj,it_str,file_path,traj_met

;_________________________________________________
;___Calculate the speed along the sun-spaceraft
;___   vector (radial speed) for interpolated values
;___   See slide 3. 
;__________________________________________________
is=where( traj_data.nh_hgi_d.r gt 0)
vr_all=dblarr(n_elements(traj_data.NH_HGI_D))
vmag_all=dblarr(n_elements(traj_data.NH_HGI_D))
vr_all(*)=-1000000.0
vmag_all(*)=-1000000.0
for i=0L, n_elements(traj_data.NH_HGI_D)-1L do begin
  if traj_data.nh_hgi_d(i).x gt -1000000.0 then begin
    v_temp=[traj_data.nh_hgi_d(i).vx,traj_data.nh_hgi_d(i).vy,traj_data.nh_hgi_d(i).vz]
    vmag_all(i)=sqrt(v_temp(0)^2+v_temp(1)^2+v_temp(2)^2)
    x_temp=[traj_data.nh_hgi_d(i).x,traj_data.nh_hgi_d(i).y,traj_data.nh_hgi_d(i).z]
    xmag=sqrt(x_temp(0)^2+x_temp(1)^2+x_temp(2)^2)
    vr_all(i)=transpose(v_temp)#x_temp/xmag
  endif
endfor 
;_________________________________________________________________
;___ Interpolate to find the speed for the time
;___     of interest 11/01/2008 18:17:38.734331
;_________________________________________________________________
it_vr=interpol( vr_all(is),time_traj(is).tyrr, it_str.tyrr) 

return,it_vr
end
;------------------------------------------------------------


;----------------------------------------------------------------------
pro get_instrument_look,vpp1,vpp2,resp,s4,wphi,eff
;----------------------------------------------------------------------
  
  ; factor in instrument response
  lphi = atan(-vpp1(1),-vpp1(0))*!radeg
  ltheta = atan(sqrt(vpp1(0)^2 + vpp1(1)^2),vpp1(2))*!radeg

  ltheta = ltheta ;- 90.0

  dphi = 0.0 ;bump instrument look direction.
  lphi = lphi+dphi

  dtheta = 90.0
  ltheta = ltheta-dtheta


  resp = get_swap_resp(vpp2,ltheta,lphi,wphi,s4,eff)


return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
pro get_e_spec,xcur,ycur,zcur,x,y,z,xp,vp,mrat,beta_p,ndx,lxyz, $
               lth,upx,clr,beta,eff,lxE,levst,tags
;----------------------------------------------------------------------
; ARGUMENTS:
; Input:
;   xcur,ycur,zcur: The current location of New Horizons
;   x,y,z: Defines the grid
;   xp,vp,mrat,beta_p,tags: Hybrid code output
; Output:
;   lxE: An array containing the centerpoints of each bin of the spectrogram (SWAP binning)
;   levst: The energy histogram with bins defined by lxE
; Who knows:
;   ndx,lxyz,lth,upx,clr,beta,eff
common fit_info,f_lxyz,f_lth,fit_arr,f_ani,s4,wphi

vr = get_NH_vr()

m1 = 1 ;hybrid simulation mass scaling
 
dx = x(1)-x(0)

dE = 1.0
hmin = 1.0
hmax = 100000.

nh = fix((hmax-hmin)/dE)+1

xsz=1000
ysz=1000

j = ycur
i = xcur
k = zcur

dV = (x(i)+dx - (x(i)-dx))*(y(j)+dx - (y(j)-dx))*(z(k)+dx - (z(k)-dx))

wh = where((xp(*,0) ge x(i)-dx) and (xp(*,0) le x(i)+dx) and $
           (xp(*,1) ge y(j)-dx) and (xp(*,1) le y(j)+dx) and $ 
           (xp(*,2) ge z(k)-dx) and (xp(*,2) le z(k)+dx) and $
           (mrat(*) le 1.0) and (tags(*) ge 1.0))


if (wh(0) gt -1) then begin
   e_arr = 0
   cnt_arr = 0
   for l = 0ll,n_elements(wh)-1 do begin

      
      vpp = reform(vp(wh(l),*))
      vpp(0) = vpp(0)+vr
      vpp2 = sqrt(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)
      vpp1 = vpp/vpp2
      vdotl = transpose(vpp1(*))#lxyz
      
      get_instrument_look,vpp1,vpp2,resp,s4,wphi,eff

      nv = vpp2/(dV*beta)

      ; energy of each macro particle within volume
      e_arr = [e_arr,(vpp(0)^2 + vpp(1)^2 + vpp(2)^2)/mrat(wh(l))]
      ; number of micro particles
      cnt_arr = [cnt_arr,nv*resp/beta_p(wh(l))]
      
   endfor
endif

e_arr = 0.5*m1*1.67e-27*e_arr*1e6/1.6e-19 ;convert to eV

; h is the energy histogram in terms of number of macro particle count
; this will be converted to micro particles in the next for loop.
h = histogram(e_arr,binsize=dE,min = hmin, max = hmax,reverse_indices=ri)
h1 = fltarr(n_elements(h))
for i = 0,n_elements(h)-2 do begin
   if (ri[i] ne ri[i+1]) then begin
      ; h1(i) is the total number of micro particles contributing to the i'th bin of the
      ; energy histogram h.
      h1(i) = total(cnt_arr(ri(ri(i):ri(i+1)-1)))
   endif
endfor

h = h1

; xE is the energy values of each bin of h in order
xE = dE*(findgen(n_elements(h)) + hmin)

;rebin to SWAP energy bins (log scale)
s = {energy_bin, e_mid: 0.0, e_min: 0.0, e_max: 0.0}
close,3
openr,3,'swap_e_bin.dat'

; These two values will end up being the last elements of the rebinned
; histogram, but they will be removed after the loop.
levst = 0
lxE = 10
while not(eof(3)) do begin
   readf,3,s
   emin = s.e_min
   emax = s.e_max      
   if ((emin ge 1.0) and (emax le 100000)) then begin
      wh = where((xE gt emin) and (xE le emax))
      if (wh(0) ge 0) then levst = [total(h(wh)),levst]
      if (wh(0) eq -1) then levst = [0,levst]
      lxE = [(emax+emin)/2,lxE]
   endif
endwhile

emin = 1
emax = 1

; Drop the last element since it doesn't come from data
lxE = lxE(0:n_elements(lxE)-2)
levst = levst(0:n_elements(levst)-2)

print,'max levst...',max(levst)


return
end
;----------------------------------------------------------------------


;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------

common fit_info,f_lxyz,f_lth,fit_arr,f_ani,s4,wphi
common NH_traj_info,traj_data,time_traj,it_str,file_path,traj_met


restore,'fin_arr_ebea_ang_eb_bg_corr.sav'
restore,'w_phi.sav'

file_path='./'

;_________________________________________________________________
;____Use the Date and Time of Interest to find the detector efficiency
;____   eff and the spacecraft speed. 11/01/2008 18:17:38.734331d
;____   Detector Efficienc described on slide 6.  
;_________________________________________________________________
idate='11/01/2008'
itime='18:17:38.734331'
CONVERT_DATE_TIME_TO_YEAR_ARR,idate,itime,it_str
ilab=idate + ' ' +itime
;_________________________________________________
;____Find Spacecraft Speed Along Radial Direction
;____   See slide 3.
;_________________________________________________

;_________________________________________________
;____Load the trajectory information to find the speed
;_________________________________________________
traj_file='20060120_20151231_600.fit' ;10 minute resolution
traj_file=file_path+traj_file
ext_names=['NH_HGI_']
READ_TRAJ, ext_names, traj_file,old_ext_names, traj_data,time_traj
traj_met=double(strmid(traj_data.nh_hgi_d.met,2,10)) +2.d*double(strmid(traj_data.nh_hgi_d.met,13,5))/1D5
i_et=double(interpol(traj_met, time_traj.tyrr, it_str.tyrr))


eff = get_dect_eff()

!p.multi=[0,1,1]

mp=1.67e-27
nfile = '1'
nfrm = 38
procnum=12
ndx = 2.0
lth = 20.0
f_lth = lth
ani = 1.
f_ani=ani

rio = 1800./40.
rpl = 1100.

xsz=1100
ysz=1000

file = 'vdist.mp4'
width = xsz
height = ysz
frames = 180
fps = 30
speed = 2
samplerate = 22050L

        ; Create object and initialize video/audio streams
oVid = IDLffVideoWrite(file)
vidStream = oVid.AddVideoStream(width, height, fps)

dirl = COMMAND_LINE_ARGS()
dir = dirl[0]

read_para,dir

restore,filename=dir+'para.sav'


read_coords,dir,x,y,z

xx = (x - x(nx/2+25))/rpl
yy = (y - y(ny/2))/rpl

device,decomposed=0
loadct,39

file = 'c.np_3d_'+strtrim(string(procnum),2)
c_read_3d_m_32,dir,file,nfrm,np

contour,np(*,*,3)/1e15<0.5,xx,yy,nlev=100,/fill,/xsty,/ysty,/isotropic

im = contour(np(*,*,3)/1e15<0.5,xx,yy,rgb_table=33,$
             xtitle='x (Rp)',$
             xstyle=1,ytitle='y (Rp)',ystyle=1,$
             xtickdir=1,ytickdir=1,dimensions=[1000,1000],axis_style=2,$
             font_size=12,/fill,n_levels=10,margin=0.15,name='espec1',aspect_ratio=1,$
            yrange=[-40,40])
im = contour(np(*,*,3)<0.5e15,xx,yy,rgb_table=33,$
             xtitle='x (Rp)',$
             xstyle=1,ytitle='y (Rp)',ystyle=1,$
             xtickdir=1,ytickdir=1,dimensions=[1000,1000],axis_style=2,$
             font_size=12,/fill,n_levels=250,margin=0.15,name='espec',aspect_ratio=1,/overplot,$
            yrange=[-40,40])

cb = colorbar(target=espec1,title='Density (cm!u-3!n)')


yy0 = -12.0
xx0 = 0.0

yy1=0.0
xx1=-44.0

slp =  (yy1-yy0)/(xx1-xx0)

xtr = xx
ytr = slp*xtr - 12.

itr = 0
wh = where(abs(ytr(0)-yy) eq min(abs(ytr(0)-yy)))
jtr = wh(0)
for i = 1,n_elements(xtr)-1 do begin
   wh = where(abs(ytr(i)-yy) eq min(abs(ytr(i)-yy)))
   itr = [itr,i]
   jtr = [jtr,wh(0)]
   print,wh(0)

endfor

plots,xx(itr),yy(jtr),/data

x0 = 145
x2 = 4
y0 =37+20
y2 = 70+20


slp = (y2-y0)/(x2-x0)


nfrm=nfrm/2

nfil = 0
xfile = dir+'c.xp_'+strtrim(string(procnum),2)
vfile = dir+'c.vp_'+strtrim(string(procnum),2)
mratfile = dir+'c.mrat_'+strtrim(string(procnum),2)
beta_p_file = dir+'c.beta_p_'+strtrim(string(procnum),2)
tags_file = dir+'c.tags_'+strtrim(string(procnum),2)


read_part,xfile,nfrm,Ni_max,xp
read_part,vfile,nfrm,Ni_max,vp
read_part_scalar,mratfile,nfrm,Ni_max,mrat
read_part_scalar,beta_p_file,nfrm,Ni_max,beta_p
read_part_scalar,tags_file,nfrm,Ni_max,tags


xcur = itr(nx-5)
ycur = jtr(nx-5)


phi = 0*!dtor  ;phi = 0 is x direction
theta = 90*!dtor
lxyz = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
f_lxyz = lxyz
upx = -403.0
get_e_spec,xcur,ycur,nz/2.,x,y,z,xp,vp,mrat,beta_p,ndx,lxyz,lth,upx,'blue',beta,eff,lxE,levst,tags

help,levst

levst_arr = fltarr(150,n_elements(levst))
levst_arr(0,*) = levst

xpl = (x - x(nx/2+25))/rpl
ypl = (y - y(ny/2))/rpl
x_arr = xpl(xcur)
y_arr = ypl(ycur)

cnt = 0

for i = itr(nx-5),itr(5),-2 do begin
   cnt = cnt+1
   print,i
   
   xcur=itr(i)
   ycur=jtr(i)
   print,'xcur,ycur...',xcur,ycur
   
   
   phi = 0*!dtor              ;phi = 0 is x direction
   theta = 90*!dtor
   lxyz = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
   f_lxyz = lxyz
   !p.multi=[0,1,1]
   upx = -403.0
   get_e_spec,xcur,ycur,nz/2.,x,y,z,xp,vp,mrat,beta_p,ndx,lxyz,lth,upx,'blue',beta,eff,lxE,levst,tags
   
   levst_arr(cnt,*) = levst
   x_arr = [x_arr,xpl(i)]
   
   contour,alog(levst_arr(0:cnt,*)>1),x_arr(0:cnt),lxE,/ylog,$
           xtitle='x (Rp)',yrange=[24,100000],$
           xstyle=1,ytitle='Energy/q [eV/q]',/fill,nlev=50

endfor

im = contour(alog(levst_arr(0:cnt,*)>1),x_arr(0:cnt),lxE,/ylog,$
             xtitle='x (Rp)',yrange=[24,100000],$;xrange=[max(x_arr),min(x_arr)],$
             xstyle=1,ytitle='Energy/q [eV/q]',ystyle=1,$
             xtickdir=1,ytickdir=1,dimensions=[500,500],axis_style=2,$
             font_size=12,/fill,rgb_table=33,n_levels=250,margin=0.15,name='espec')


cb = colorbar(target=espec,title='Counts')
cb.Save,"espec.pdf",border=30

xpos1=0.15
xpos2= 0.85
ypos1 = 0.2
ypos2 = 0.8

evst_r = alog(levst_arr(1:cnt,*)>1)
evst_r(0,0) = max(evst_r);10.8
evst = bytscl(evst_r)
wh = where(evst eq 0b)
evst(wh) = 255b

img_cont_ylog,evst(*,*),x_arr(*),lxE(*),dunit,xpos1,xpos2,ypos1,ypos2,evst_r,/postscript

cnt_arr = levst_arr(0:cnt,*)
x = x_arr(0:cnt)
y = lxE

save,filename='espec.sav',cnt_arr,x,y


end
;----------------------------------------------------------------------
