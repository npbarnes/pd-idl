;----------------------------------------------------------------
PRO read_part,ffile,nfrm,Ni_max,xp
;----------------------------------------------------------------
; read the nfrm'th frame of the vector particle file 'file' with Ni_max particles.
; store result in xp

    frm=0l

    file = ffile+'.dat'
    print,' reading...',file
    openr,1,file,/f77_unformatted
    xp=fltarr(Ni_max,3,/nozero)


    readu,1,frm
    print,'  image #.....',frm
    readu,1,xp
    frmcnt = 1

    while (frmcnt lt nfrm) do begin

       readu,1,frm
       print,'  image #.....',frm
       readu,1,xp
       frmcnt = frmcnt + 1

    endwhile

    close,1

    return
end

;----------------------------------------------------------------


;----------------------------------------------------------------
PRO read_part_scalar,ffile,nfrm,Ni_max,xp
;----------------------------------------------------------------
; read the nfrm'th frame of the scalar particle file 'file' with Ni_max particles.
; store result in xp

    frm=0l

    file = ffile+'.dat'
    print,' reading...',file
    openr,1,file,/f77_unformatted
    xp=fltarr(Ni_max,/nozero)

    readu,1,frm
    print,'  image #.....',frm
    readu,1,xp
    frmcnt = 1

    while (frmcnt lt nfrm) do begin

       readu,1,frm
       print,'  image #.....',frm
       readu,1,xp
       frmcnt = frmcnt + 1

    endwhile

    close,1

    return
end
;----------------------------------------------------------------

;------------------------------------------------------------
function get_swap_resp,ee,theta,phi,w,s4,eff
;------------------------------------------------------------
    ;Compute effective area
    aeff=0.033
    aeff=aeff/0.0882d
    dect_eff=eff
    aeff=dect_eff*aeff ;cm^2
    aeff = aeff/1e10   ;km^2

    ;Tansmission coef at angle phi
    wp=interpol(w.w, w.phi,phi)

    ;get x bin from theta 
    junk=min(abs(theta-(-s4.x)),ix)

    ;get energy bin of ee
    junk=min(abs(ee-s4.ecen),iee)

    ;Energies represented within the bin
    er1=s4.y(iee,*)*ee

    ;Transmission for this theta within the energy bin
    tt=s4.arr(*,ix,iee)

    ;Tansmission for this theta with that energy
    ttnew=interpol(tt,er1, ee)
    if (ee lt min(er1) or ee gt max(er1)) then ttnew = 0.0d

    ;SWAP response (transmission as a function of phi)*(transmission as a function of E and theta)*(effective area)
    res = wp*ttnew*aeff

    return,res
end
;------------------------------------------------------------



;------------------------------------------------------------
function get_dect_eff, eff
;------------------------------------------------------------
common NH_traj_info,traj_data,time_traj,it_str,file_path,traj_met


    ;____Loading in the Dectector Efficiency vs Time
    ; ____Set Types
    ; _____byte 1
    ; _____int 2
    ; _____long 3
    ; _____float 4
    ; _____double 5
    ; _____complex 6
    ; _____string 7 
    fieldtypes=dblarr(3)
    csvfile=file_path+'er_table_lookup_3_2015c.csv'
    fieldtypes(*)=7
    fieldtypes(0)=5
    fieldtypes(2)=5
    fin=READ_CSV_FIELDTYPES(strcompress(csvfile,/remove_all), fieldtypes)
    convert_date_time_to_year_arr_spice, fin.utc, junktime,ef
    met_mid=double(interpol(traj_met, time_traj.tyrr, it_str.tyrr))

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


    return,edetector
end
;------------------------------------------------------------


;------------------------------------------------------------
function get_NH_vr
;------------------------------------------------------------
    common NH_traj_info,traj_data,time_traj,it_str,file_path,traj_met

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

    return,it_vr
end
;------------------------------------------------------------

; These are the right-handed passive rotation matricies
; (i.e. rotate the coodinates with a fixed vector)
; Take their transpose (or equivelently let alpha -> -alpha) 
; for the right-handed active roataion.
; (i.e. rotate the vector with fixed coodinates)
function rotXmat, alpha
  return, [[1.0,0.0,0.0],[0.0,cos(alpha),sin(alpha)],[0.0,-sin(alpha),cos(alpha)]]
end

function rotYmat, alpha
  return, [[cos(alpha),0.0,-sin(alpha)],[0.0,1.0,0.0],[sin(alpha),0.0,cos(alpha)]]
end

function rotZmat, alpha
  return, [[cos(alpha),sin(alpha),0.0],[-sin(alpha),cos(alpha),0.0],[0.0,0.0,1.0]]
end

;----------------------------------------------------------------------
pro get_instrument_look, orientation, v1, lphi, ltheta
;----------------------------------------------------------------------
; ARGUMENTS:
; Input:
;   vpp1: Normalized velocity vector of the macro paricle in space coordinates xyz
;           i.e. NOT the velocity with respect to the instrument coordinates that
;           are also called xyz
;   vpp2: Velocity squared of the macro particle
;   eff: Detector efficiency
;   wphi: A structure containing phi values(bins) and the transmission coef at each
;           of those phi values
; Output:
;   resp: The SWAP response to the macroparticle
; Who knows:
;   s4
    ; spacecraft orientation (i.e. sunward direction)

    stheta = orientation[0]
    sphi = orientation[1]
    spin = orientation[2]

    ; convert direction of particle motion in pluto coords, vpp1, to direction of particle
    ; motion in SWAP coords
    ; first rotate +y to +x (this puts it in swap coords for orientation (theta=0,phi=0,spin=0)
    v1_swap = rotZmat(-90.0*!DtoR)##transpose(v1) 
    ; then apply rotations of the spacecraft in turn
    v1_swap = rotZmat(sphi*!DtoR)##rotXmat(stheta*!DtoR)##rotYmat(spin*!DtoR)##v1_swap

    ; get the look direction
    v1_swap_look = -v1_swap

    ; get angles of the look direction
    lphi = atan(v1_swap_look(0),v1_swap_look(1))*!radeg
    ltheta = -atan(v1_swap_look(2),sqrt(v1_swap_look(0)^2 + v1_swap_look(1)^2))*!radeg



    return
end
;----------------------------------------------------------------------

pro get_NH_local_particles, xp, x, y, z, radius, tags, mrat, species, weights
    hydrogen = where((sqrt( (xp(*,0)-x)^2 + (xp(*,1)-y)^2 + (xp(*,2)-z)^2 ) le radius) and $
               (mrat(*) eq 1.0) and $
               (tags(*) ne 0.0), /null)
    helium = where((sqrt( (xp(*,0)-x)^2 + (xp(*,1)-y)^2 + (xp(*,2)-z)^2 ) le radius) and $
               (mrat(*) eq 0.5) and $
               (tags(*) ne 0.0), /null)
    methane = where((sqrt( (xp(*,0)-x)^2 + (xp(*,1)-y)^2 + (xp(*,2)-z)^2 ) le radius) and $
               (mrat(*) eq 1.0/16.0) and $
               (tags(*) ne 0.0), /null)

    species = List(hydrogen, helium, methane)
    weights = [1.0, 2.0, 1.0]
end

;----------------------------------------------------------------------
pro make_e_spectrum,xcur,ycur,zcur,xp,vp,mrat,beta_p, $
               beta,eff,bins,tags,particles,spectrum, orientation
;----------------------------------------------------------------------
; ARGUMENTS:
; Input:
;   xcur,ycur,zcur: The current location of New Horizons
;   xp,vp,mrat,beta_p,beta,tags: Hybrid code output
;   bins: An array containing the left endpoints of each bin of the spectrogram (SWAP binning)
;   eff: Detector efficiency
;   heavy: Set keyword to restrict output to heavy particles
; Output:
;   levst: The energy histogram of particle counts
    common fit_info,s4,wphi

    mp = 1.6726215E-27 ; kg
    e =  1.6021766E-19 ; C 

    vr = get_NH_vr()

    m1 = 1 ;hybrid simulation mass scaling
     
    dE = 1.0
    hmin = 1.0
    hmax = 100000.

    nh = fix((hmax-hmin)/dE)+1

    xsz=1000
    ysz=1000

    radius = 2000.0
    dV = (4.0/3.0)*!DPI*radius^3

    count = 0l


    e_arr = 0
    cnt_arr = 0
    for l = 0ll,n_elements(particles)-1 do begin

      
      ;particle velocity
      v = reform(vp(particles(l),*))
      ;add spacecraft velocity
      v(0) = v(0)+vr
      ;magnitude
      vmag = sqrt(v(0)^2 + v(1)^2 + v(2)^2)
      ;velocity direction unit vector
      v1 = v/vmag
      ;energy/charge in electron volts per fundamental charge
      ;equivilent to joules per coulomb
      ee=.5*((mp/e)/mrat(particles(l)))*(vmag*1000.)^2
      
      get_instrument_look, orientation, v1, lphi, ltheta
      resp = get_swap_resp(ee,ltheta,lphi,wphi,s4,eff)

      nv = vmag/(dV*beta*beta_p(particles(l)))

      ; energy of each macro particle within volume
      e_arr = [e_arr,ee]
      ; number of micro particles for each macro particle
      cnt_arr = [cnt_arr,nv*resp]
      
    endfor

    ; build histogram of microparticle counts
    spectrum = fltarr(n_elements(bins))
    foreach bin, bins, i do begin
       part_in_bin = where(e_arr ge bin.e_min and e_arr lt bin.e_max, count)
       if (count ne 0) then begin
           spectrum(i) = total(cnt_arr(part_in_bin))
       endif
    endforeach

    return
end
;----------------------------------------------------------------------
pro get_pluto_position, para, pluto_position
    if (tag_exist(para,"pluto_offset")) then begin
        pluto_position = para.qx(n_elements(para.qx)/2 + para.pluto_offset)
    endif else begin
        print, "pluto offset assumed 30"
        pluto_position = para.qx(n_elements(para.qx)/2 + 30)
    endelse
end

function y_from_x, x, para, p0, p1

    rpl = para.RIo

    maxx = para.qx(-1)
    maxy = para.qy(-1)
    maxz = para.qz(-1)

    slp =  (p1.y-p0.y)/(p1.x-p0.x)

    y = -slp*x + p0.y*rpl + maxy/2

    return, y
end

pro make_flyby_e_spectrograms, dir, p, t1, t2, xp, vp, mrat, beta_p, eff, bins, tags,$
    spectrograms_by_species, times, positions, orientations
;----------------------------------------------------------------------
; ARGUMENTS:
; Input:
;   p: parameter structure.
;   x: Structure in the form {x0:*, x1:*, y0:*, y1:*}. Two points that NH passes through on its trajectory.
;       In units of Rp with pluto at the origin. x-y pluto coordinates.
; Output:
;   levst_arr: Synthetic energy spectrogram

    common python_interpreter, NH_tools
    ; radius of pluto
    rpl = 1187.

    delta_t = 60.0D
    num_spectra = floor((t2-t1)/delta_t) + 1

    maxx = p.qx(-1)
    maxy = p.qy(-1)
    maxz = p.qz(-1)


    ; Build a histogram of micro particle counts using the SWAP energy bins (log scale)
    ; First read what the SWAP bins are.
    readbins, bins
    ; We now have the bin values
    lxE = bins.e_mid

    spectrograms_by_species = fltarr(3,num_spectra,n_elements(bins))

    cnt = 0

    zcur = maxz/2.0
    radius = 2000.0
    get_pluto_position, p, pluto_position


    i = 0
    times = dblarr(num_spectra)
    positions = dblarr(num_spectra,3)
    orientations = dblarr(num_spectra, 3)
    for t = t1, t2, delta_t do begin

       
       xcur = NH_tools.pos_at_time(t)
       ycur=y_from_x(xcur, p, {point, x:0.,y:12.}, {point, x:158., y:-30.})
       print, "Progress: ", (1 - (t2-t)/(t2-t1))*100, "%", xcur, ycur
       o = NH_tools.orientation_at_time(t)
       
       times[i] = t
       positions[i,0] = xcur
       positions[i,1] = ycur
       positions[i,2] = zcur
       orientations[i,0] = o[0]
       orientations[i,1] = o[1]
       orientations[i,2] = o[2]

       get_NH_local_particles, xp, xcur+pluto_position, ycur, zcur, radius, tags, mrat, species, weights

       for j=0, 2 do begin
           make_e_spectrum,xcur+pluto_position,ycur,zcur,xp,vp,mrat,beta_p,p.beta,eff,bins,tags,species(j),spec,o
           spectrograms_by_species(j,i,*) = spec*weights(j)
       endfor

       i = i+1
    endfor
    print, "Progress: 100%"
end

pro get_ave_v, beta_p, vp, vave
    vave = [0.,0.,0.]
    
    ; dV and volume will cancel out in the next calculation anyway
    ; so I just left them out.
    num_parts = 1./(beta_p(*))

    for i=0,2 do begin
        vave(i) = total(num_parts(*)*vp(*,i))/total(num_parts(*))
    endfor
    print, vave
end

;----------------------------------------------------------------------
;main program
;----------------------------------------------------------------------

common fit_info,s4,wphi
common NH_traj_info,traj_data,time_traj,it_str,file_path,traj_met
common python_interpreter, NH_tools

NH_tools = Python.Import('NH_tools')
t1 = NH_tools.close_start
t2 = NH_tools.close_end


restore,'fin_arr_ebea_ang_eb_bg_corr.sav'
restore,'w_phi.sav'

file_path='./'

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


eff = get_dect_eff()


args = COMMAND_LINE_ARGS()
dir = args[0]

read_para,file_dirname(dir),p

procnum=ceil(float(p.num_proc)/2.0)

file = 'vdist.mp4'
width = 1100l
height = 1000l
fps = 30

; Create object and initialize video/audio streams
oVid = IDLffVideoWrite(file)
vidStream = oVid.AddVideoStream(width, height, fps)
device,decomposed=0
loadct,39

if(tag_exist(p,"part_nout")) then begin
    nfrm=p.nt/p.part_nout
endif else begin
    print, "No part_nout: assuming 30 frames"
    nfrm=30
endelse
nfrm=20

xfile = dir+'c.xp_'+strtrim(string(procnum),2)
vfile = dir+'c.vp_'+strtrim(string(procnum),2)
mratfile = dir+'c.mrat_'+strtrim(string(procnum),2)
beta_p_file = dir+'c.beta_p_'+strtrim(string(procnum),2)
tags_file = dir+'c.tags_'+strtrim(string(procnum),2)

read_part,xfile,nfrm,p.Ni_max,xp
read_part,vfile,nfrm,p.Ni_max,vp
read_part_scalar,mratfile,nfrm,p.Ni_max,mrat
read_part_scalar,beta_p_file,nfrm,p.Ni_max,beta_p
read_part_scalar,tags_file,nfrm,p.Ni_max,tags

make_flyby_e_spectrograms, dir, p, t1, t2, xp, vp, mrat, beta_p, eff, bins, tags,$
                           spectrograms_by_species, times, positions, orientations

;count = 1
;for i=nfrm-1, nfrm-2, -1 do begin
;    read_part,xfile,i,p.Ni_max,xp
;    read_part,vfile,i,p.Ni_max,vp
;    read_part_scalar,mratfile,i,p.Ni_max,mrat
;    read_part_scalar,beta_p_file,i,p.Ni_max,beta_p
;    read_part_scalar,tags_file,i,p.Ni_max,tags
;
;    make_flyby_e_spectrograms, dir, p, traj, xp, vp, mrat, beta_p, eff, bins, tags, new_sp_by_sp
;    spectrograms_by_species += new_sp_by_sp
;    count += 1
;endfor
;spectrograms_by_species = spectrograms_by_species/count


readbins, ebins
ebins = ebins.e_mid

save, description=dir, filename='espec-'+args[1]+'.sav',spectrograms_by_species,times,positions,orientations,ebins

end
;----------------------------------------------------------------------
