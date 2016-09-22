; plots unformatted FORTRAN 
; file:  name of file
; n,m:   dimensions of variable


PRO fcontour,file,n,m

openr,1,file,/f77_unformatted
arr=fltarr(n,m,/nozero)
readu,1,arr
contour,arr
close,1

end