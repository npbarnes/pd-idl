pro read_coords,dir,qx,qy,qz

  close,1
  openr,1,dir+'c.coord.dat',/f77_unformatted
  nnx = 0
  readu,1,nnx

  nny = 0
  readu,1,nny

  nnz = 0
  readu,1,nnz

  qx = fltarr(nnx)
  readu,1,qx

  qy = fltarr(nny)
  readu,1,qy

  qz = fltarr(nnz)
  readu,1,qz

  close,1
return
end
