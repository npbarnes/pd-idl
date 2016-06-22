function count_rate_c2,vc, p

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
et=dblarr(n_elements(vc))
vsc=dblarr(n_elements(vc))
c=0.d
d=0.d
for i=0, n_elements(vc)-1  do begin 
   readf, unit2, format='(3D30.8)',a, b, c
   et(i)=c
endfor

cl0=p(7)
w=360.d/double(p(6))
cl=cl0+w*(double(double(et)-double(p(8)) ))
clr=cl*!dpi/180.D
thetac=p(3) +p(5)*cos(clr+p(9))
phic=p(4) +p(5)*cos(clr+p(9)+!dpi/2.d)

thetac=thetac*degrad
phic=phic*degrad
close, unit2
free_lun, unit2
res=dblarr(n_elements(vc))
for i=0, n_elements(vc)-1 do begin
 temp=CALCULATE_FLUX_D(vc(i),thetac(i),phic(i), wphi, s4,eff)
 res(i)=temp
endfor

return, res
end
