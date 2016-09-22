; read 3d coordinate data 
PRO c_read_2dxy_vec_m_32,dir,file,nfrm,x_xy,x_xz
;file='coord.dat'

  close,1
  restore,filename=dir+'para.sav'

  file = dir+file+'.dat'
  print,' reading...',file
  openr,1,file,/f77_unformatted

  x_xy=fltarr(nx,ny,3,/nozero)
  x_xz=fltarr(nx,nz,3,/nozero)

  frm = 0l
  readu,1,frm
  print,'  image #.....',frm
  readu,1,x_xz,x_xy
  frmcnt = 1

  while (frmcnt lt nfrm) do begin
     
     frm=0l
     readu,1,frm
     print,'  image #.....',frm
     readu,1,x_xz,x_xy
     frmcnt = frmcnt + 1
     
  endwhile
  
  close,1
  
end
