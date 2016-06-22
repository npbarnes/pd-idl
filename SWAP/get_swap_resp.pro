;------------------------------------------------------------
function get_swap_resp,vc,theta,phi,p,s4,eff, vnew
;------------------------------------------------------------

  aeff=0.033; cm^2 from McComas et al., [2008] 
  ;________________________________________________
  ;____Here aeff is the effective area 
  ;____    including the geometric and detector 
  ;____    efficiency on the launch date.
  ;________________________________________________
  aeff=aeff/0.0882d  ; denominator from Phil Valek 
  ;________________________________________________
  ;____Now the aeff above is only the geometrical 
  ;____    portion of the effective area.
  ;________________________________________________
  dect_eff=eff
  ;________________________________________________
  ;___eff is the detector efficiency
  ;___    it is a function of time see slide 5 and 6
  ;________________________________________________
  aeff=dect_eff*aeff

  ;________________________________________________
  ;____Here aeff is the effective area 
  ;____    including the geometric and detector 
  ;____    efficiency on the date you specified 
  ;____    in the main routine
  ;________________________________________________
  mp=1.67262158D-27
  kb=1.380658D-23
  e=1.6E-19 
  ee=.5*mp/e*(vc*1000.)^2 ; convert step speed to step energy
  enew=.5*mp/e*(vnew*1000.0)^2
  ;________________________________________________
  ;___ Interpolte p to find the value at phi specified
  ;___   in the main routine. See slide 8.
  ;___________________________________________________
  p_phi=interpol(p.w, p.phi,phi)

  ;___________________________________________________
  ;_____Find the energy step index so you can look up 
  ;_____ the theta phi 
  ;___________________________________________________
  junk=min(abs(ee-s4.ecen),iee)

  ;___________________________________________________
  ;____Extract angle information for w response
  ;___________________________________________________
  azimuth_lab=s4.x
  theta_lab=-azimuth_lab

  ;___________________________________________________
  ;_____w(E/Es, theta, Es) =s4.array
  ;___________________________________________________
  w=s4.arr
  ;____________________________________________________
  ;_____In the lab they use an angle called azimuth
  ;_____  theta as defined on slide 1 has the opposite
  ;_____  sign of azimuth.
  ;_____Find index (ix) for w(E/Es, theta, Es)  that corresponds to
  ;_____ the theta of interest specifed in the main program
  ;____________________________________________________
  junk=min(abs(theta-(theta_lab)),ix)

  ;____________________________________________________
  ;_____Find the ambient energy for the w response
  ;____________________________________________________
  er1=s4.y(iee,*)*ee

  ;________________________________________________________
  ;____Interpolating the energy to find the response at enew(vnew)
  ;__________________________________________________________
  ttnew=interpol(w(*,ix,iee),er1, enew)
  if enew lt min(er1) and enew gt max(er1) then ttnew=0.0d
  res = p_phi*ttnew*aeff
  return, res
end
;------------------------------------------------------------
