pro read_traj, ext_names, traj_file,old_ext_names, traj_data,time_traj
;_________________________
;_____Read Trajectory Information
;_________________________

old_ext_names=ext_names
master_names=['NH_HGI_','NH_HAE_DATE_','NH_HAE_J2000_','NH_HG_','NH_HEE_','NH_HEEQ_','NH_J2000_','NH_PLUTO_J2000_', 'NH_JUPITER_J2000_', 'NH_PLUTO_IAU_', 'NH_JUPITER_IAU_']

master_ext_num=findgen(n_elements(master_names))+1
ext_num=dblarr(n_elements(ext_names))
for i=0, n_elements(ext_names)-1 do begin 
    iw=where(ext_names(i) eq master_names)
    if iw(0) eq -1 then print, ext_names(i) +' is not a valid coordinate system.'
    if iw(0) gt -1 then ext_num(i)=master_ext_num(iw) else ext_num(i)=-99999
endfor
iv=where(ext_num gt -99999)
if iv(0) gt -1 then begin 
  ext_names=ext_names(iv)
  ext_num=ext_num(iv)
  str_traj='traj_data={'
  loop_lim= n_elements(ext_num)-1 
  for i=0, loop_lim do begin 
    ext=ext_num(i)
    str_data=ext_names(i)+'d=mrdfits(traj_file, ext_num(i), header,/unsigned,/silent)'
    junk=execute(str_data)
    if i le  n_elements(ext_num)-2 then begin 
      str_traj=str_traj+ ext_names(i)+'d:' +ext_names(i)+'d,'
    endif else begin 
     str_traj=str_traj+ ext_names(i)+'d:' +ext_names(i)+'d}'
    endelse
  endfor 
 junk=execute(str_traj)
 str_temp=' CONVERT_DATE_TIME_TO_YEAR_ARR_SPICE, traj_data.'+ext_names(0)+'d.utc, junk, time_traj'
 print, str_temp
 junk=execute(str_temp)
endif
  


end
