;-----------------------------------------------------------------
function get_swap_res_hybrid,vc,theta,phi,w,s4,eff
;-----------------------------------------------------------------

file_path='/Users/delamere/projects/SWAP/swap_calibration_model/'
;_________________________________________________________________
;____Reference Material is provided in the file called
;____    swap_instrument_response_b.pptx.
;____  
;____To load SWAP instrument response and spacraft speed information
;____   Note that the detector efficiency and the spacecraft speed
;____       are a function of time. 
;_________________________________________________________________

;_________________________________________________________________
;____Use the Date and Time of Interest to find the detector efficiency
;____   eff and the spacecraft speed. 11/01/2008 18:17:38.734331d
;____   Detector Efficienc described on slide 6.  
;_________________________________________________________________
idate='11/01/2008'
itime='18:17:38.734331'
CONVERT_DATE_TIME_TO_YEAR_ARR,idate,itime,it_str
ilab=idate + ' ' +itime
;_________________________________________________
;____Find Spacecraft Speed Along Radial Direction
;____   See slide 3.
;_________________________________________________

;_________________________________________________
;____Load the trajectory information to find the speed
;_________________________________________________
traj_file='20060120_20151231_600.fit' ;10 minute resolution
traj_file=file_path+traj_file
ext_names=['NH_HGI_']
READ_TRAJ, ext_names, traj_file,old_ext_names, traj_data,time_traj
traj_met=double(strmid(traj_data.nh_hgi_d.met,2,10)) +2.d*double(strmid(traj_data.nh_hgi_d.met,13,5))/1D5
i_et=double(interpol(traj_met, time_traj.tyrr, it_str.tyrr))

;_________________________________________________
;___Calculate the speed along the sun-spaceraft
;___   vector (radial speed) for interpolated values
;___   See slide 3. 
;__________________________________________________
is=where( traj_data.nh_hgi_d.r gt 0)
vr_all=dblarr(n_elements(traj_data.NH_HGI_D))
vmag_all=dblarr(n_elements(traj_data.NH_HGI_D))
vr_all(*)=-1000000.0
vmag_all(*)=-1000000.0
for i=0L, n_elements(traj_data.NH_HGI_D)-1L do begin
  if traj_data.nh_hgi_d(i).x gt -1000000.0 then begin
    v_temp=[traj_data.nh_hgi_d(i).vx,traj_data.nh_hgi_d(i).vy,traj_data.nh_hgi_d(i).vz]
    vmag_all(i)=sqrt(v_temp(0)^2+v_temp(1)^2+v_temp(2)^2)
    x_temp=[traj_data.nh_hgi_d(i).x,traj_data.nh_hgi_d(i).y,traj_data.nh_hgi_d(i).z]
    xmag=sqrt(x_temp(0)^2+x_temp(1)^2+x_temp(2)^2)
    vr_all(i)=transpose(v_temp)#x_temp/xmag
  endif
endfor 
;_________________________________________________________________
;___ Interpolate to find the speed for the time
;___     of interest 11/01/2008 18:17:38.734331
;_________________________________________________________________
it_vr=interpol( vr_all(is),time_traj(is).tyrr, it_str.tyrr) 


;_______________________________________________________________
;____aeff is the effective area and has cm^2 units. See slide 5.
;_______________________________________________________________

;_______________________________________________________________
;____eff (now labeled edetector to match slides)  is the
;____   detector efficiency and is dimensional less, it varies
;____   with time as the gain for the detectors and detector
;____   operational voltage is adjusted (slides 6).
;_______________________________________________________________

;_______________________________________________________________
;_____I fold edetector and aeff together into one because 
;_____   in the lab they typically measure the combination. 
;_______________________________________________________________

;________________________________________________
;____Loading in the Dectector Efficiency vs Time
;________________________________________________
; __________________________
; ____Set Types
; _________________________
; _____byte 1
; _____int 2
; _____long 3
; _____float 4
; _____double 5
; _____complex 6
; _____string 7 
fieldtypes=dblarr(3)
csvfile='er_table_lookup_3.csv'
fieldtypes(*)=7
fieldtypes(0)=5
fieldtypes(2)=5
fin=READ_CSV_FIELDTYPES(strcompress(csvfile,/remove_all), fieldtypes)
convert_date_time_to_year_arr_spice, fin.utc, junktime,ef
met_mid=double(interpol(traj_met, time_traj.tyrr, it_str.tyrr))

;________________________________________________________
;____Interpolate to find the correct detector efficiency
;____   a specific time.  A plot is give on slide 6.
;____   Note that you have to calibraiton data taken 
;____   with the current operational voltage to 
;____   calculate the current efficiency. The table
;____   has been set up to handle this issue. 
;_______________________________________________________
met_junk=min(abs(met_mid-fin.met),imet)
if fin.met(imet) gt met_mid then iarr=[imet-1, imet] else iarr=[imet, imet+1]
edetector=fin.eff(iarr(0))+(fin.eff(iarr(1))-fin.eff(iarr(0)))/(fin.met(iarr(1))-fin.met(iarr(0)))*(double(met_mid)-fin.met(iarr(0)))
if fin.met(imet) eq met_mid then edetector = fin.eff(imet)
;_________________________________________________________________

;_________________________________________________________________
;______Load Example of SWAP energy steps and elapsed times at 11/01/2008 18:17:38.734331
;______   Note that the energy SWAP detects is not the ambient energy.
;______   NH is moving away from the Sun SWAP detects a lower solar
;______   wind energy and speed than the ambient. Spacecraft speed
;______   is described on slide 3.
;__________________________________________________________________
LOAD_EXAMPLE, einst, st_et, sp_et

;_______________________________________________
;_____Converting energies to speed
;_______________________________________________
e=1.6D-19 
mp=1.67262158D-27
vinst=sqrt(2.d*einst*e/mp)/1000.d
sz_mod=n_elements(vinst)
et=(double(st_et)+double(sp_et))/2.D  ; find middle time

;_________________________________________________________________
;_____Load the w(E/Es,theta, Es) data. An example of this array at
;_____   ~1000 eV is given on slide 7. 
;_____   We created a table of w(E/Es,theta) for the coarse scan
;_____   energy steps (Es).  This data is stored in a 3-D array
;_____   called s4.array. This array is 
;_________________________________________________________________
;restore,'fin_arr_ebea_ang_eb_bg_corr.sav'

;_________________________________________________________________
;_____Load the p(phi) data; top plot on slide 8
;_________________________________________________________________
;restore,'w_phi.sav'
;_________________________________________________________________

;______________________________________
;___Energy Ratio Range limits of integration
;______________________________________
min_rat=min(s4.y)
max_rat=max(s4.y)
;______________________________________
;___Field of View limits of integration in equation on slide 4
;______________________________________
min_theta=-6.0
max_theta=6.0
min_phi=-276.0/2.0
max_phi=276.0/2.0
npts_phi=138.0
npts_theta=30.0
npts_speed=100.0
theta=min_theta+findgen(npts_theta)*(max_theta-min_theta)/(npts_theta-1.d)
phi=min_phi+findgen(npts_phi)*(max_phi-min_phi)/(npts_phi-1.d)
;____________________________________________________________________
;____Below I construct an array called "response" combines 
;____  Aeff*w*p information  into 1 array as funciton of 
;____  phi, theta, and speed for a specitived time. 
;____  response(phi, theta, speed)
;____
;____For each energy step the fuction "response" holds Aeff*w*p in the
;____   equation on slide 4.
;____Aeff has a time dependence because the detector efficiency varies
;____   with time.
;____
;____To perform the integral on slide 4 you need to shift your
;____   distribution function into the instrument frame using 
;____   the spacecraft velocity information (slide 3) at the
;____   time you are trying to simulate. You have to have the
;____   speed (energy) at which SWAP would see it to be.
;____
;____Here I'm using the energy into the SWAP instrument; the
;____   energy that SWAP sees (detects).
;____
;____Also note that here the theta and phi are referring to the
;____   variables that are integrated. When the spacecraft is spinning
;____   the solar wind beam has a different theta_o and phi_o at which
;____   the beam peaks at every energy step. The solar wind beam
;____   will essentially make a circle in theta-phi space when 
;____   the spacecraft is spinning. I actually
;____   fit the theta and phi sun angles (derived from the attitude)
;____   and on a theta and phi plot it is extremely circular.
;____   See slides 11 and 12. 

;___________________________________________________________________
for i = 0,n_elements(vinst)-1 do begin  ;instrument energy step loop 
   ;____________________________________________
   ;___Finding the speed (energy) limits of integration
   ;___for the passband of a given energy step
   ;___in the equation on slide 4.
   ;______________________________________________
   energy_min=min_rat*einst(i)
   energy_max=max_rat*einst(i)
   speed_min=sqrt(2.d*energy_min*e/mp)/1000.d; km/s
   speed_max=sqrt(2.d*energy_max*e/mp)/1000.d; km/s
   response=dblarr(npts_phi, npts_theta, npts_speed)
   for k=0, npts_theta-1 do begin ; theta loop spanning full FOV
     for l=0, npts_phi-1 do begin ; phi loop spanning full FOV
       r = 0.0
       energy=0.0
       speed=0.0
       for j=0, npts_speed-1 do begin ;energy passband loop spanning full energy range
         vnew=speed_min+j*(speed_max-speed_min)/(npts_speed-1.d)
         res=GET_SWAP_RESP(vinst(i),theta(k),phi(l),wphi,s4,edetector, vnew)
         response(l, k, j)=res
         r = [r,res]
         speed=[speed, vnew]
         enew=.5*mp/e*(vnew*1000.0)^2
         energy=[energy, enew]
         print,'SWAP response...',i,j,l, k, vinst(i), einst(i), vnew, enew,  theta(k), phi(l), res
       endfor ;energy pass band loop
       set_plot, 'X'
       window, 3, retain=2
       plot, energy,r,ytitle='$w(E,\theta)p(\phi)A_{eff}$',xtitle='Energy [eV/q]', title='Energy Step: '+strcompress(string(einst(i),format='(F12.3)'),/remove_all)+' [eV/q]'+' Theta: '+strcompress(string(theta(k),format='(F12.3)'),/remove_all)+' [deg]'+' Phi: '+strcompress(string(phi(l),format='(F12.3)'),/remove_all)+' [deg]', /xlog, xrange=[20,10000], pos=[.1-.03,.15,.45-.03,.8],/normal,charsize=0.85, xstyle=1
       plot, speed,r,ytitle='$w(E,\theta)p(\phi)A_{eff}$',xtitle='Speed [km/s]', title='Speed Step: '+strcompress(string(vinst(i),format='(F12.3)'),/remove_all)+' [km/s]'+' Theta: '+strcompress(string(theta(k),format='(F12.3)'),/remove_all)+' [deg]'+' Phi: '+strcompress(string(phi(l),format='(F12.3)'),/remove_all)+' [deg]', xrange=[50,1500], pos=[.6-.03,.15,.95-.03,.8],/normal,/noerase,charsiz=0.85, xstyle=1
     endfor ;phi loop
  endfor; theta loop
endfor;energy step loop

  
end
