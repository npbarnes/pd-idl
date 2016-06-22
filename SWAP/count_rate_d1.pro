function count_rate_d1,vc, p


;_____________________________
;____phi respose
;____________________________
restore, filename='w_phi.sav'
;_________________________
;__weenie curves
;________________________
restore,filename='fin_arr_ebea_ang_eb_bg_corr.sav'
get_lun, unit1
openw, unit1, 'inputs.txt'
degrad=!dpi/180.D

no=p(0)
vo=p(1)
to=p(2)
dth=p(3)
dph=p(4)
eff=p(10)

mp=1.67262158D-27
kb=1.380658D-23
vtho=sqrt((2.0D)*kb*to/mp)/1000.d
ango=asin(vtho/vo)*180./!dpi
vo=vo*1000.d *100.d
printf, unit1, format='(3F30.8)', no, vo, to;, thetao, phio
close, unit1
free_lun, unit1


get_lun, unit2
openr , unit2, 'angles.txt'
degrad=!dpi/180.D
thetac=dblarr(n_elements(vc))
phic=dblarr(n_elements(vc))
et=dblarr(n_elements(vc))
c=0.d
a=0.0d
b=0.0d
for i=0, n_elements(vc)-1  do begin 
   readf, unit2, format='(3D30.8)',a, b, c
   et(i)=c
   thetac(i)=(a+dth)*degrad
   phic(i)=(b+dph)*degrad
endfor

close, unit2
free_lun, unit2
res=dblarr(n_elements(vc))
for i=0, n_elements(vc)-1 do begin
 temp=CALCULATE_FLUX_D(vc(i),thetac(i),phic(i), wphi, s4,eff)
 res(i)=temp
endfor

return, res
end
