;__________________________________________________________________
;____Purpose: To load and plot the SWAP instrument response
;____        functions and use that to plot simulated
;____         count rates for a given sweep.This code uses 
;____         a Maxwellian distribution.  See documentation in 
;____         swap_instrument_response_b.pptx. 
;__________________________________________________________________
plot_path='plots/'
scape='portrait'
set_plot, 'X'
device, true_color=24
device, decompose=0
loadct, 39, /silent
output_type='P'
traj_file='20060120_20151231_600.fit' ;10 minute resolution
;traj_file='20060120_20151231_3600.fit' ;1 hour resolution
;traj_file='20060120_20151231_86400.fit' ;1 day resolution
file_path=''
;_________________________________________________________________
;____Date and Time of Interest
;_________________________________________________________________
idate='11/01/2008'
itime='18:17:38.734331'
CONVERT_DATE_TIME_TO_YEAR_ARR,idate,itime,it_str
ilab=idate + ' ' +itime
;_________________________________________________________________
;____Date range for time series speed and detector efficiency plot
;_________________________________________________________________
start_date='01/20/2006'
start_time='00:00:00'
stop_date='01/01/2016'
stop_time='00:00:00'
;___________________________________________________________________
;____Convert start and stop time for time sereies plots to year and doy
;___________________________________________________________________
CONVERT_DATE_TIME_TO_YEAR_ARR, start_date,start_time,st
CONVERT_DATE_TIME_TO_YEAR_ARR, stop_date, stop_time, sp
plot_label='START:'+start_date + '('+string(st.daynumber,format='(I03)') +') '+start_time + '!CSTOP:  '+ stop_date + '('+string(sp.daynumber,format='(I03)') +') '+stop_time +'!C'
;_________________________________________________
;____Find Spacecraft Speed Along Radial Direction
;_________________________________________________
;____Load the trajectory information
;_________________________________________________
traj_file=file_path+traj_file
ext_names=['NH_HGI_']
READ_TRAJ, ext_names, traj_file,old_ext_names, traj_data,time_traj
traj_met=double(strmid(traj_data.nh_hgi_d.met,2,10)) +2.d*double(strmid(traj_data.nh_hgi_d.met,13,5))/1D5
st_met=double(interpol(traj_met, time_traj.tyrr, st.tyrr))
sp_met=double(interpol(traj_met, time_traj.tyrr, sp.tyrr))
i_et=double(interpol(traj_met, time_traj.tyrr, it_str.tyrr))
plot_label=plot_label +'MET: '+strcompress(string(st_met,format='(F17.2)'),/remove_all)+ ' to ' + strcompress(string(sp_met,format='(F17.2)'),/remove_all)
;_________________________________________________
;___Calculate the speed along the sun-spaceraft
;___   vector (radial speed) for interplated values
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
;______________________________________________________
;____Plot a time series of the spacecraft radial speed
;____       and the total spacecraft speed.
;______________________________________________________
 ytitle='Speed [km/s]'
it_vr=interpol( vr_all(is),time_traj(is).tyrr, it_str.tyrr)
PLOT_VR_SC_B,time_traj(is), vmag_all(is),vr_all(is), plot_path,scape, output_type, '', [0,45], [st.tyrr,sp.tyrr],plot_label,poso, ilab, it_str,it_vr,ytitle 
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
;____   a specific time
;_______________________________________________________
met_junk=min(abs(met_mid-fin.met),imet)
if fin.met(imet) gt met_mid then iarr=[imet-1, imet] else iarr=[imet, imet+1]
eff_fact=fin.eff(iarr(0))+(fin.eff(iarr(1))-fin.eff(iarr(0)))/(fin.met(iarr(1))-fin.met(iarr(0)))*(double(met_mid)-fin.met(iarr(0)))
if fin.met(imet) eq met_mid then eff_fact = fin.eff(imet)
;_______________________________________________________
;___Plot a time series of the detector efficiency
;___     Added cross hairs for the specified time of
;___     11/01/2008 18:17:38.734331.
;________________________________________________________
yrange=[.01,.2]
xrange=[st.tyrr,sp.tyrr]
ytitle='<Nc>!U2!N/(<Ns><Np>)'
PLOT_EFF,  ef.tyrr, fin.eff,  plot_path,scape, output_type, '', yrange,xrange,plot_label,pos, ilab, it_str, eff_fact, ytitle
;______________________________________________________
;____Calculate Overall Effective Area 
;____(Geometric Effective Area + Detector Efficiency)
;____For specified time 11/01/2008 18:17:38.734331
;______________________________________________________
aeff=0.033 ;McComas et al. 2008
aeff=aeff/0.0882d ; Geometric Part of Effective Area
dect_eff=eff_fact; detector efficiency
aeff=dect_eff*aeff; combined detector efficiency and geometric effective area
;__________________________________
;____Load Phi Response 
;____loading structure called wphi
;___________________________________
restore, filename='w_phi.sav'
PLOT_PHI_RESPONSE, plot_path, output_type, scape, wphi
;___________________________________
;____Loading Energy-Theta Resonse
;____Loading Structure called s4
;___________________________________
restore,filename='fin_arr_ebea_ang_eb_bg_corr.sav'
do_norm=1
if do_norm eq 1 then begin
  ct_label='Log(Norm Coin Cnt Rate)'
  tag='norm'
endif
;_____________________________________________
;____Plot interpolated array as a function of Eb/Ea
;_____________________________________________
file=plot_path+'energy_theta_response.ps'
  if scape eq 'portrait' then begin 
  !x.ticklen=.025
  xdim=2.2
  ydim=2.2
  x_off=1.1
  y_off=8.1
  xin=8.5;6.692913
  yin=11;8.976378
  file_name=file
  color_code=0
  PLOT_TYPE_POSITION_C, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type,file, color_code 
endif else begin 
  !x.ticklen=.05;.05
  xdim=2.2
  ydim=2.2
  x_off=1.5
  y_off=4.5
  xin=11.0
  yin=8.0
  file_name=file
  color_code=0
  plot_type_position_cland, yin, xin, xdim, ydim, x_off, y_off, ymn, xmn, ymx, xmx, lg_sz, mid_sz, sm_sz, thk, rad, output_type, FILE,color_code
endelse

itt=where(finite(alog10(s4.arr)) eq 1)

;______________________________________________________________________
;____The energy-theta response array is s4.arr(energy sub element,theta, energy step)
;____
;____The first element is Ebeam/Eanalyzer where Ebeam refers to the
;____calibration beam energy or the energy of the ambient beam entering
;____then instrument.
;____
;____The second element refers to the azimuth angle. The theta angle
;____    is the negative of the azimuth angle.  theta=-azimuth=-s4.x
;____
;____The 3rd element is the analyzer energy step (Eanalyzer) eanalyzer=s4.ecenter
;____   The values for the 3rd element in s4.arr are s4.ecen
;____  
;____
;___________________________________________________________________ 
azimuth=s4.x
theta=-s4.x
;___________________________________________________________________ 
;___Find the energy and theta response array needed for a given operation voltage
;___junk=min(abs(s4.ecen - energy_step_desired),iz)
;___The 2 resonse_energy_theta=s4.arr(energy,theta, iz)
;___     where theta=-s4.x
;___           energy=s4.ey*energy_step_desired
;___________________________________________________________________ 
;____________________________________
;_____Create Color Plots of the Array
;_____Only plotting every 100th one 
;____________________________________
;____Color bar scale
mmin=-3
mmax=0.5
increment=100  ; set to 1 to see plots for all energies.
for i=0,1000-1,100 do begin 
  PLOT_COLOR_4_EB_EA_FINAL,s4.ecen, i, azimuth,xmn, xmx, ymn, ymx, s4.arr, s4.y, ct_label,mmin,mmax
endfor

;___________________________________________________________________
;____Example instrument coincidence count rate energy distribution
;____   for 11/01/2008 18:17:38.734331.
;____   Here just grabbed and pasted the energy and time information
;____   to do a quick example.
;___________________________________________________________________
;_____Define input array for the function used to model the count rate 
;_____     distribution
;_____For spinning intervals the solar wind makes essentially a circle
;_____      in SWAP FOV which can essentially be represented as a
;_____      circle in SWAP theta and phi angles where theta and
;_____      phi correspond to x and y in the circle parametric 
;_____      equations.
;_____To model spinning intervals use count_rate_c2.pro
;_____ The ET time information is read in through the angles.txt file.
;______p(0) -density
;______p(1) -speed as measured ambient speed - radial spacecraft speed
;______p(2) -temperature
;______p(3) -origin of solar wind circle in swap theta [deg]
;______p(4) -origin of solar wiind circle in swap phi [deg]
;______p(5) -radius of solar wind circle in [deg]
;______p(6) -spin period [sec] typically 5rpm ~12 sec
;______p(7) -reference angle in [deg]
;______p(8) - reference ephemeris time in ET [sec]
;______p(9) - spin phase [deg]
;______p(10) -detector efficiency at the time of interest 

;_____To model 3-axis intervals use count_rate_d1.pro and the theta
;_____and phi angles for each energy step are stored along with the ET time
;_____in the angles.txt file
;______p(0) -density
;______p(1) -speed as measured ambient speed - radial spacecraft speed
;______p(2) -temperature
;______p(3) -delta theta, solar wind offset from Sun theta angle
;______p(4) -delta phi, solar wind offset from Sun phi angle
;______p(5) -N/A
;______p(6) -N/A
;______p(7) -N/A
;______p(8) -N/A
;______p(9) -N/A
;______p(10) -detector efficiency at the time of interest 
;____________________________________________________________
;_____Example at 11/01/2008 18:17:38.734331 is spinning.
;_____     Using count_rate_c2.pro
;_____For p using the detector efficiency at  11/01/2008 18:17:38.734331
p=[0.0513505764d, 452.882382297230d,  8934.1155226985d,  -0.1075242000d, $
0.1669704000d, 3.6409709942d, 12.0424680985d, 9.6645848020d, $
278835532.1118609309d,   0.0857176207d, eff_fact]  
;_______________________________________________________________
;____Reduce the ambient wind speed by the spacecraft speed at 
;____      11/01/2008 18:17:38.734331 to
;____      simulate the measured speed. 
;_______________________________________________________________
;***************************
p(1)=p(1)-it_vr;************
;***************************

;_________________________________________________________________
;______Load Example at 11/01/2008 18:17:38.734331
;__________________________________________________________________
LOAD_EXAMPLE, einst, st_et, sp_et

;_______________________________________________
;_____Converting energies to speed
;_______________________________________________
e=1.6D-19 
mp=1.67262158D-27
vinst=sqrt(2.d*einst*e/mp)/1000.d
sz_mod=n_elements(vinst)
et=(double(st_et)+double(sp_et))/2.D
close, 111

;_________________________________________________________________________
;_______Writing out time information that the count_rate_c2.pro routine
;_______reads in to create a model rate based on a Maxwellian.
;__________________________________________________________________________ 
openw,111,'angles.txt'
a=0
b=0
for i=0,n_elements(vinst)-1 do begin 
  printf,111, format='(3D30.8)', a, a, et((i))
endfor
close, 111
;_________________________________________________________________
;_______Calling count_rate_c2.pro to calcuate simulated count rates
;_______   using a Maxwellian and the various SWAP Response Functions
;_______   This is for spinning cases.
;_______
;_______count_rate_c2.pro Summary:
;_______   1) loads the w and response curves on 
;_______      slides 7 and 8.
;_______   2) the detector efficiency at the time of interest
;_______      is input via 10th parameter in p(10) and is held fixed
;_______      during the fitting process so it is a constant and not 
;_______      a varied fit parameter.
;_______   3) The parameters for the Maxwellian and sun circle
;_______      created via spinning are given in p the order is
;_______      specified above. When fitting the data these other
;_______      parameters are variables.
;_______   4) This subroutine routine defines the circle using
;_______      parameteric equations for a circle.
;_______   5) A subroutine called calculate_flux_d.pro is called to
;_______      calculate the integrand for the equation on slide 4.
;_______   6) The subroutine calculate_flux_d.pro combines the
;_______      geometric effective area with the detector efficiency
;_______      to get the overall effective area, calculates the
;_______      Maxwellian (via flux.pro), and preforms the integration.
;_______   7) The subroutine calculate_flux_d.pro calls a subroutine
;_______      called flux.pro that calculates the Maxwellian and 
;_______      the cos(theta) and V^3 part of the integrad in 
;_______      the equation on slide 4. 
;_________________________________________________________________
sim_count_spin= COUNT_RATE_C2(vinst, p)
if i eq 0 then begin
   ystyle=1
   xstyle=1
endif else begin 
   ystyle=5
   xstyle=5
 endelse 
;_________________________________________________________________________
;_______Writing out time, theta, and phi information that the count_rate_d1.pro routine
;_______   reads in to create a model rate based on a Maxwellian.
;_______Just setting the sun theta to 3.6 deg and phi to 3.6 deg 
;_______   for every energy step in the sweep. That is the sun is in
;_______   the same location for every swap measurement in the sweep
;_______   pair. This information is transfered via a file called angles.txt.
;_______I also assume that the solar wind has a fixed offset
;_______   frome the sun direction and the offset theta and phi angles
;_______   are p(3) and p(4) respectively.
;__________________________________________________________________________ 
openw,111,'angles.txt'
for i=0,n_elements(vinst)-1 do begin 
  printf,111, format='(3D30.8)', 3.6, 3.6, et((i))
endfor
close, 111
;_________________________________________________________________
;_______Calling count_rate_d1.pro to calcuate simulated count rates
;_______   using a Maxwellian and the various SWAP Response Functions
;_______   This is 3-axis stablized cases. For this example I assume 
;_______   the pointing is fixed for the entire sweep pair. 
;_______
;_______count_rate_d1.pro Summary:
;_______   1) loads the w and response curves on 
;_______      slides 7 and 8.
;_______   2) the detector efficiency at the time of interest
;_______      is input via 10th parameter in p(10) and is held fixed
;_______      during the fitting process so it is a constant and not 
;_______      a varied fit parameter.
;_______   3) The parameters for the Maxwellian in p where the order
;_______      is specified above. When fitting the data, these
;_______      parameters are variables.
;_______   4) The Sun theta and phi location for each energy
;_______      is  read in via the angle.txt file.
;_______   5) The solar wind offsets in theta and phi are stored in
;_______      p(3) and p(4) respectively. When fitting the data,  
;_______      both of these parameters are variables.
;_______   5) A subroutine called calculate_flux_d.pro is called to
;_______      calculate the integrand for the equation on slide 4.
;_______   6) The subroutine calculate_flux_d.pro combines the
;_______      geometric effective area with the detector efficiency
;_______      to get the overall effective area, calculates the
;_______      Maxwellian (via flux.pro), and preforms the integration.
;_______   7) The subroutine calculate_flux_d.pro calls a subroutine
;_______      called flux.pro that calculates the Maxwellian and 
;_______      the cos(theta) and V^3 part of the integrad in 
;_______      the equation on slide 4. 
;___________________________________________________________________
sim_count= COUNT_RATE_D1(vinst, p)
PLOT_COUNT_RATE_DISTRIBUTION, plot_path,output_type, scape, sim_count_spin, sim_count, einst

end
