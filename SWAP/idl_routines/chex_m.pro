pro chex_m,nfiles,POSTSCRIPT=postscript

if keyword_set(postscript) then begin
  set_plot,'ps
  !p.font=0
  device,/palatino
  device,/encapsulated
  device,filename='chex.eps'
;  device,/landscape
end


n=0
dt=0.004
chex_rate=0.0
t=0.2
time=0.2

close,1
for m=1,nfiles do begin

   file='chex_'+strtrim(string(m),1)+'.dat'

   openr,1,file,/f77_unformatted
;   readu,1,n
;   print,n

   while (not(eof(1))) do begin
      readu,1,y 
      chex_rate=[chex_rate,y]
      t=t+dt
      time=[time,t]
   endwhile

close,1
endfor

time=time(1:*)
;chex_rate=median(chex_rate(1:*),3)

;!p.multi=[0,1,2]
plot,time,chex_rate,ytitle = 'dN!dce!n/dt', xtitle = 't (s)', $
     title='Charge exchange rate',charsize=1.4
;plot_io,time,chex_rate,ytitle = 'dNi/dt', xtitle = 't (s)', $
;     title='Charge exchange rates',charsize=1.5,ystyle=1,$
;     yrange=[1e22,1e25]

;xyouts,0.15,0.9,'% released material = ' $
;       +strtrim(string(100.*total(chex_rate*dt)/2.6e24),1),/normal

print,100.*total(chex_rate*dt)/2.6e24

;No=2.4e24
;tau=28.0

;dnidt = (No/tau)*exp(-time/tau)
;oplot,time,dnidt,linestyle=1

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
   !p.font=-1
endif
end