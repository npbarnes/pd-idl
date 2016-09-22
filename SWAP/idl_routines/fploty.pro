PRO fploty,file,y

n=0

openr,1,file,/f77_unformatted
readu,1,n
print,n
y=fltarr(n,/nozero)
readu,1,y
plot,y,psym=1
close,1

return
end