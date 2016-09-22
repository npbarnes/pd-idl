; plots unformatted FORTRAN 
; file:  name of file
; n:     dimensions of variable


PRO fplotxy,file,n

openr,1,file,/f77_unformatted
x=fltarr(n)
y=fltarr(n)
readu,1,x
readu,1,y
plot,x,y
close,1

end