pro satnp_m,nfiles,POSTSCRIPT=postscript

!p.font=-1
if keyword_set(postscript) then begin
  !p.font=0
  set_plot,'ps
  device,filename='satnp.ps',/palatino
  device,/landscape
end


n=0
dt=0.004
satnp_arr=0.0
t=0.2
time=0.2

close,1
for m=1,nfiles do begin

   file='satnp_'+strtrim(string(m),1)+'.dat'

   openr,1,file,/f77_unformatted
   print,'opening...',file
;   readu,1,n
;   print,n

   while (not(eof(1))) do begin
      readu,1,y 
      satnp_arr=[satnp_arr,y]
      t=t+dt
      print,'np...',t,y   
      time=[time,t]
   endwhile

close,1
endfor

time=time(1:*)
;satnp_arr=median(satnp_arr(1:*),3)
satnp_arr=smooth(satnp_arr,10)

plot,time,satnp_arr/1e20,ytitle = 'n!IBa!n [10!e20!n km!e-3!n]', xtitle = 'release time (sec)', $
     title='G1:  Barium density at satellite position',charsize=1.5

if keyword_set(postscript) then device,/close
!p.font=-1
end

