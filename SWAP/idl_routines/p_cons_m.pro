PRO p_cons_m,nfiles,dt,tot_surf,POSTSCRIPT=postscript

;mom_m.pro
;plots momentum as a function of time for hybrid code

if keyword_set(postscript) then begin
  set_plot,'ps
  !p.font=7
  device,/times
  device,filename='p_cons.ps'
  device,xsize=6.5,ysize=9.0,xoffset=1.0,yoffset=1.0,/inches
end

;file=''
close,1
openr,1,'p_conserve_1.dat',/f77_unformatted

nt=0ll
nout=0ll
ts=0ll
sts=0ll
;dt = 0.00005
;dt = 0.2
surf_tot=dblarr(3,/nozero)
graduu_tot=dblarr(3,/nozero)
ugradu_tot=dblarr(3,/nozero)

readu,1,nt
readu,1,nout

print,nt,nout,nt/nout
readu,1,ts
readu,1,sts
tstep = ts
;print,'time step #.....',ts,sts
readu,1,surf_tot,graduu_tot,ugradu_tot
print,surf_tot,graduu_tot,ugradu_tot

surf_totx=surf_tot(0)
surf_toty=surf_tot(1)
surf_totz=surf_tot(2)

graduu_totx=graduu_tot(0)
graduu_toty=graduu_tot(1)
graduu_totz=graduu_tot(2)

ugradu_totx = ugradu_tot(0)
ugradu_toty = ugradu_tot(1)
ugradu_totz = ugradu_tot(2)

while (not(eof(1))) do begin
   readu,1,ts
   readu,1,sts
   tstep = [tstep,(ts-1)+sts]
   print,'time step #.....',(ts-1)*10+sts

   readu,1,surf_tot,graduu_tot,ugradu_tot

   surf_totx = [surf_totx,surf_tot(0)]
   surf_toty = [surf_toty,surf_tot(1)]
   surf_totz = [surf_totz,surf_tot(2)]

   graduu_totx = [graduu_totx,graduu_tot(0)]
   graduu_toty = [graduu_toty,graduu_tot(1)]
   graduu_totz = [graduu_totz,graduu_tot(2)]

   ugradu_totx = [ugradu_totx, ugradu_tot(0)]
   ugradu_toty = [ugradu_toty, ugradu_tot(1)]
   ugradu_totz = [ugradu_totz, ugradu_tot(2)]

endwhile

close,1

nfile = 2
while (nfile le nfiles) do begin
file = 'p_conserve_'+strtrim(string(nfile),1)+'.dat'
close,1
openr,1,file,/f77_unformatted

readu,1,nt
readu,1,nout

print,nt,nout,nt/nout

while (not(eof(1))) do begin
   readu,1,ts
   readu,1,sts
   tstep = [tstep,(ts-1)+sts]
   print,'time step #.....',(ts-1)+sts

   readu,1,surf_tot,graduu_tot,ugradu_tot

   surf_totx = [surf_totx,surf_tot(0)]
   surf_toty = [surf_toty,surf_tot(1)]
   surf_totz = [surf_totz,surf_tot(2)]

   graduu_totx = [graduu_totx,graduu_tot(0)]
   graduu_toty = [graduu_toty,graduu_tot(1)]
   graduu_totz = [graduu_totz,graduu_tot(2)]

   ugradu_totx = [ugradu_totx, ugradu_tot(0)]
   ugradu_toty = [ugradu_toty, ugradu_tot(1)]
   ugradu_totz = [ugradu_totz, ugradu_tot(2)]

endwhile
close,1
nfile = nfile+1
endwhile

nx = n_elements(surf_totx)

stotx = dblarr(nx)
stoty = dblarr(nx)
stotz = dblarr(nx)

gtotx = dblarr(nx)
gtoty = dblarr(nx)
gtotz = dblarr(nx)

utotx = dblarr(nx)
utoty = dblarr(nx)
utotz = dblarr(nx)

stotx(0) = surf_totx(0)*dt
stoty(0) = surf_toty(0)*dt
stotz(0) = surf_totz(0)*dt

gtotx(0) = graduu_totx(0)*dt
gtoty(0) = graduu_toty(0)*dt
gtotz(0) = graduu_totz(0)*dt

utotx(0) = ugradu_totx(0)*dt
utoty(0) = ugradu_toty(0)*dt
utotz(0) = ugradu_totz(0)*dt

for i=1,nx-1 do begin
   stotx(i) = stotx(i-1) + surf_totx(i)*dt
   stoty(i) = stoty(i-1) + surf_toty(i)*dt
   stotz(i) = stotz(i-1) + surf_totz(i)*dt

   gtotx(i) = gtotx(i-1) + graduu_totx(i)*dt
   gtoty(i) = gtoty(i-1) + graduu_toty(i)*dt
   gtotz(i) = gtotz(i-1) + graduu_totz(i)*dt

   utotx(i) = utotx(i-1) + ugradu_totx(i)*dt
   utoty(i) = utoty(i-1) + ugradu_toty(i)*dt
   utotz(i) = utotz(i-1) + ugradu_totz(i)*dt

endfor




!p.multi=[0,3,4]
plot,tstep,stotx,title='surf_totx',charsize=1.5
plot,tstep,stoty,title='surf_toty',charsize=1.5
plot,tstep,stotz,title='surf_totz',charsize=1.5
plot,tstep,gtotx,title='graduu_totx',charsize=1.5
plot,tstep,gtoty,title='graduu_toty',charsize=1.5
plot,tstep,gtotz,title='graduu_totz',charsize=1.5
plot,tstep,utotx,title='ugradu_totx ',charsize=1.5
plot,tstep,utoty,title='ugradu_toty',charsize=1.5
plot,tstep,utotz,title='ugradu_totz',charsize=1.5

;totpx = (stotx+gtotx)
totpx = (stotx)
plot,tstep,totpx,title='total x',charsize=1.5

;totpy = (stoty+gtoty) 
totpy = (stoty) 
plot,tstep,totpy,title='total y',charsize=1.5

;totpz = (stotz+gtotz)
totpz = (stotz)
plot,tstep,totpz,title='total z',charsize=1.5

!p.multi=[0,1,1]

if keyword_set(postscript) then begin
   device,/close
   set_plot,'x
endif

tot_surf=dblarr(n_elements(tstep),3)
tot_surf(*,0) = totpx
tot_surf(*,1) = totpy
tot_surf(*,2) = totpz

return
end






