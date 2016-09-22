;-------------------------------------------------------------------
pro np_prof,nfiles,nfrm,maxnp
;-------------------------------------------------------------------

nfrm = nfiles*nfrm
xinteranimate,set=[500,400,nfrm]
WINDOW, 0, XSIZE=500, YSIZE=400, TITLE='Ion Distro'
f_read_coord,'coord.dat',x,y,z,dzg,dzc,nx,ny,nz
nprof = dblarr(nx)

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
         nprof(i) = total(icld(i,*,*))/1.0e23
      endfor

      ri=4
      plot,x-(ri*(x(1)-x(0))),nprof,xtitle='x (km)',$
           ytitle = '10!e23!n ions per 0.4 km',$
           charsize = 1.5,title='Ion distribution along satellite path',$
           yrange = [0,maxnp],ystyle=1
      xinteranimate,frame=frm_cnt-1,window=0

   endwhile
   close,1

endfor

xinteranimate

return
end
;-------------------------------------------------------------------
