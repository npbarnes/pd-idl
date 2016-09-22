;-------------------------------------------------------------------
pro np_prof_th_g1,nfiles,nfrm,maxnp,POSTSCRIPT=postscript
;-------------------------------------------------------------------

nfrm = nfiles*nfrm
;xinteranimate,set=[500,400,nfrm]
WINDOW, 0, XSIZE=500, YSIZE=700, TITLE='Ion Distro'
f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz
ri=4
x = x-(ri*(x(1)-x(0)))
print,x
nprof = dblarr(nx)

if keyword_set(postscript) then begin
set_plot,'ps
!p.font=0
device,filename='iondistro_g1.eps'
device,/palatino
device,/portrait
device,/encapsulated
device,/inches,xsize=6.0,xoffset=1.0,ysize=8.0,yoffset=1.0
end

close,1
file = 'npall_'

openr,1,file+'1.dat',/f77_unformatted

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

close,1

icld = fltarr(nx,ny,nz,/nozero)
ncld = fltarr(nx,ny,nz,/nozero) 

frm = 1
frm_cnt = 0
dy = 0.8/float(nfrm)
print,'dy...',dy
for j = 0,nfiles-1 do begin

   mfile = ((frm-1)/(nfrm/nfiles)) + 1
   print,mfile,frm
   frmcount=1
   frmcount=1+((frm-1) mod (nfrm/nfiles))
   openr,1,file+strtrim(string(mfile),1)+'.dat',/f77_unformatted
   readu,1,nt
   readu,1,nout
   readu,1,nx
   readu,1,ny
   readu,1,nz

   while not(eof(1)) do begin

      readu,1,frame
      print,file+strtrim(string(mfile),1)+' image #.....',frame
      readu,1,icld
      frm_cnt = frm_cnt + 1
      frm = frm + 1
      for i = 0,nx-1 do begin
         nprof(i) = total(icld(i,*,*))/1.0e22
      endfor

      !p.position=[0.1,0.1+dy*(nfrm - frm_cnt),0.9,0.1+dy*(nfrm - frm_cnt) + dy]
      print,!p.position
      !x.ticks = 1
;      !x.tickv = [min(x),max(x)]
      !x.tickv = [0,15]
      !x.tickname = [' ',' ']
      !x.title = ' '
      !x.range = [0,15]
      !x.style = 1
      !y.style = 1
      !y.range = [0,maxnp]
      !y.title = ' ' 
      !y.ticks = 3
      !y.tickv = [0,10,20,maxnp]
      !y.tickname = [' ','10','20',' ']
      !p.title = ' '	 
      !p.subtitle = ' ' 
      if (frm_cnt eq 1) then !p.title = 'Ion distribution along satellite path'
      if (frm_cnt eq nfrm) then begin
         !x.title = 'x (km)'
	 !x.ticks = 3
	 !x.tickv = [0,5,10,15]
         !x.tickname = ['0','5','10','15']
      endif 
      if (frm_cnt eq nfrm/2) then !y.title = '10!e22!n ions per 0.4 km'
      plot,x,nprof, charsize = 1.2,/noerase
      tm = strmid(strtrim(string(0.2+(frm_cnt)*25*0.004),2),0,3)
      xyouts,0.8,0.1+dy*(nfrm - frm_cnt)+0.6*dy,tm+' s',charsize=1.2,/normal
;      xinteranimate,frame=frm_cnt-1,window=0

   endwhile
   close,1

endfor

;xinteranimate

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
   !p.font=-1
endif

return
end
;-------------------------------------------------------------------
