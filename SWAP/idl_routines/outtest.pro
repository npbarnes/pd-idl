close,1

a=fltarr(10,/nozero)

openr,1,'outtest.dat',/f77_unformatted

readu,1,a

close,1

print,a

end