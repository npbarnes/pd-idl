; read 2d arrays
PRO c_read_2dxy_m_32,dir,file,nfrm,x_xy,x_xz

  close,1
  restore,filename=dir+'para.sav'

  file = dir+file+'.dat'
  openr,1,file,/f77_unformatted
  print,'reading....',file

  frm = 0l
  x_xy=fltarr(fix(nx),fix(ny),/nozero)
  x_xz=fltarr(fix(nx),fix(nz),/nozero)
  readu,1,frm
  print,'  image #.....',frm
  readu,1,x_xz,x_xy
  frmcnt = 1
  
  while (frmcnt lt nfrm) do begin
     frm = 0l
     readu,1,frm
     print,'  image #.....',frm
     readu,1,x_xz,x_xy
     frmcnt = frmcnt + 1
  endwhile
  
  close,1

  return
end
