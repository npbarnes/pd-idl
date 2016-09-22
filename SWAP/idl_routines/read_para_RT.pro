pro read_para_RT,dir

  close,1
  openr,1,dir+'para.dat',/f77_unformatted

  nx = 0l
  ny = 0l
  nz = 0l
  readu,1,nx,ny,nz,dx,dy,delz
  print,nx,ny,nz,dx,dy,delz


  nt = 0l
  ntsub = 0l
  nout = 0l
  readu,1,nt,dtsub_init,ntsub,dt,nout
; print,nt,dtsub_init,ntsub,dt,nout


  out_dir = '                               '

  readu,1,out_dir
  print,out_dir

;  model_choice = '             '

;  readu,1,model_choice
;;  print,model_choice

;  readu,1,nf_init,b0_init
;  print,nf_init,b0_init
  
;  readu,1,nu_init,lww2,lww1
;  print,nu_init,lww2,lww1

;  readu,1,Mdot,Mdot_part
;  print,Mdot,Mdot_part

  readu,1,vtop,vbottom
  print,vtop,vbottom
stop

  Ni_max = 0l
  readu,1,Ni_max
  
  mproton = 0.0d
  m_fluid=0.0d
  m_pu=0.0d
  readu,1,mproton,m_fluid,m_pu

  tempf0 = 0.0d
  readu,1,tempf0

  RIo = 0.0d
  readu,1,RIo

  alpha = 0.0d
  beta = 0.0
  readu,1,alpha,beta

  comm_sz = 0l
  io_proc = 0l
  readu,1,comm_sz,io_proc

  close,1
  

  save,filename=dir+'para.sav',nx,ny,nz,dx,dy,delz,$
       nt,dtsub_init,ntsub,dt,nout,$
       out_dir,$
       model_choice,$
       nf_init,b0_init,$
       nu_init,lww2,lww1,$
       Mdot,Mdot_part,$
       vsw,$
       Ni_max,$
       mproton,m_fluid,m_pu,$
       tempf0,$
       RIo,$
       alpha,beta,$
       comm_sz,io_proc


return
end
