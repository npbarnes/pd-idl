;-----------------------------------------------------------------------
pro get_stream_lines,data,whx,why,whz,ny,nz,outverts,outconn,outnormals,$
                     vertcolors




return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro Io_STREAM,LINES = lines, TUBES = tubes

temp_val = 150 ;eV
den_val = 0.08 ;c
b_val = 0.4e-9

Rio = 1200.
ctbl = 13
ach = 0.75

;file_dir = '/Volumes/MacD97-2/hybrid/3d_buf/run_test/'
file_dir = '/Volumes/Scratch/hybrid/Pluto/run_test4/'
;file_dir = '/Volumes/MacD97-2/hybrid/3d_buf/pleiades/'
;file_dir = '/Volumes/MacD97-2/hybrid/3d_buf/janus/run_1/'
;file_dir = '/Volumes/MacD97-2/hybrid/Io/local/'
;restore,file_dir+'Xianzhe_100km_nt400.sav'
;restore,file_dir+'XIANZHE_3D_local.sav'
;restore,file_dir+'LARGE_3D_local.sav'
restore,'/Volumes/Scratch/Dols_output/trajectory_0_24_27_31.sav'
restore,file_dir+'3d_stitch_Pluto.sav'


;ufarr(*,*,*,0) = smooth(uiarr(*,*,*,0),2)
;ufarr(*,*,*,1) = smooth(uiarr(*,*,*,1),2)
;ufarr(*,*,*,2) = smooth(uiarr(*,*,*,2),2)

b1arr(*,*,*,0) = smooth(b1arr(*,*,*,0),2)
b1arr(*,*,*,1) = smooth(b1arr(*,*,*,1),2)
b1arr(*,*,*,2) = smooth(b1arr(*,*,*,2),2)

nparr = smooth(nparr,2)



x = x/Rio
y = y/Rio
z = z/Rio


;minx = min(x)
;maxx = max(x)

minx = -100.0
maxx = 200.0
 
miny = -150.0
maxy = 150.0

;miny = min(y)
;maxy = max(y)

;minz = min(z)
;maxz = max(z)


minz = -200
maxz = 200



;define subvolume

whx = where((x gt minx) and (x lt maxx))
why = where((y gt miny) and (y lt maxy))
whz = where((z gt minz) and (z lt maxz))

x = x(whx)
y = y(why)
z = z(whz)


; Read the global map file data.
;IMAGE = read_png('./io_image.png', r, g, b)
IMAGE = read_png('./pluto.png', r, g, b)

MESH_OBJ, 4, vertices, polygons, REPLICATE(1.0, 101, 101)  

oModel = OBJ_NEW('IDLgrModel')  

xax = OBJ_NEW('IDLgrAxis',direction=0,range=[minx,maxx],location=[0,0,0],$
              /notext,/exact,thick=3,ticklen=1)
yax = OBJ_NEW('IDLgrAxis',direction=1,range=[miny,maxy],location=[0,0,0],$
              /notext,/exact,thick=3,ticklen=1)
zax = OBJ_NEW('IDLgrAxis',direction=2,range=[minz,maxz],location=[0,0,0],$
              /notext,/exact,thick=3,ticklen=1)

xtit = OBJ_NEW('IDLgrText','x',baseline = [1.0,0,0], $
               locations=[-35,5,0],$
               alignment=0.5,/enable_formatting,color=!color.black,$
              char_dimensions = [10,10])

ytit = OBJ_NEW('IDLgrText','y',baseline = [0,1.0,0], $
               updir = [-1.0,0,0], locations=[-5,-75,0],$
               char_dimensions=[10,10],$
               alignment=0.5,/enable_formatting,color=!color.black)
oModel->Add, ytit

ztit = OBJ_NEW('IDLgrText','z',baseline = [0,0,1.0], $
               updir = [-1.0,0, 0],locations=[-5,0,95],$
               char_dimensions=[10,10],$
               alignment=0.5,/enable_formatting,color=!color.black)
oModel->Add, ztit


;oModel->Add,xax
;oModel->Add,yax
;oModel->Add,zax

oPalette = OBJ_NEW('IDLgrPalette')  
oPalette -> LOADCT, ctbl  
oPalette -> SetRGB, 255, 255, 255, 255  
oImage = OBJ_NEW('IDLgrImage', image, PALETTE = oPalette)  

vector = FINDGEN(101)/100.  
;vector = shift(vector,20)
texure_coordinates = FLTARR(2, 101, 101)  
texure_coordinates[0, *, *] = vector # shift(REPLICATE(1., 101),0)  
texure_coordinates[1, *, *] = REPLICATE(1., 101) # vector  

oPolygons = OBJ_NEW('IDLgrPolygon', SHADING = 1, $  
   DATA = vertices, POLYGONS = polygons, $  
   COLOR = [255, 255, 255], $  
   TEXTURE_COORD = texure_coordinates, $  
   TEXTURE_MAP = oImage, /TEXTURE_INTERP)  

;oPolygons->SetProperty, SHADER=oShader

;oModel -> ADD, oPolygons  

;add Jupiter
;Jup_img = read_png('./jupiter_css.png', r, g, b)
;MESH_OBJ, 4, vertices, polygons, REPLICATE(10.0, 101, 101)  
;Jup_Image = OBJ_NEW('IDLgrImage', Jup_img, PALETTE = oPalette)  
;vertices(1,*) = vertices(1,*)+2.5*10.0
;vertices(2,*) = vertices(2,*)-0.5*10.0

;oJupiter = OBJ_NEW('IDLgrPolygon', SHADING = 1, $  
;   DATA = vertices, POLYGONS = polygons, $  
;   COLOR = [255, 255, 255], $  
;   TEXTURE_COORD = texure_coordinates, $  
;   TEXTURE_MAP = Jup_Image, /TEXTURE_INTERP)  
;oModel -> ADD, oJupiter


data = FLTARR(3,n_elements(whx), n_elements(why), n_elements(whz))

data(0,*,*,*) = smooth(uparr(whx,why,whz,0),4)
data(1,*,*,*) = smooth(uparr(whx,why,whz,1),4)
data(2,*,*,*) = smooth(uparr(whx,why,whz,2),4)


ystep = 10
zstep = 10
nseeds = 6*LONG((ny*20*nz)/(ystep*zstep))
seeds = FLTARR(nseeds)
iseed=0L
;whx = where(x gt minx)
for k = 0,n_elements(whz)-1 do begin
   for j = 0,n_elements(why)-1 do begin
      if ((k mod zstep) eq 0) and ((j mod ystep) eq 0) and $
         (iseed lt (nseeds-2)) and (y(j) gt -40) and (y(j) lt 40) $
         and (z(k) gt -50) and (z(k) lt 50) then begin
         seeds(iseed) = FLOAT(0)
         seeds(iseed+1) = FLOAT(j)
         seeds(iseed+2) = FLOAT(k)
         iseed = iseed+3
      endif
   endfor
endfor

maxIterations=4000
stepSize=2.0

PARTICLE_TRACE,data,seeds,outverts,outconn,outnormals, $
               MAX_ITERATIONS=maxIterations, MAX_STEPSIZE=stepSize,  $
               INTEGRATION=0,ANISOTROPY=[1,1,1], SEED_NORMAL=[0, 0, 1]

outnormals(0,*)= median(reform(outnormals(0,*)),100)
outnormals(1,*)= median(reform(outnormals(1,*)),100)
outnormals(2,*)= median(reform(outnormals(2,*)),100)

mon = sqrt(outnormals(0,*)^2 + outnormals(1,*)^2 + outnormals(2,*)^2)
outnormals(0,*) = outnormals(0,*)/mon
outnormals(1,*) = outnormals(1,*)/mon
outnormals(2,*) = outnormals(2,*)/mon


magdata = SQRT(data[0,*, *,*]^2 + data[1,*, *,*]^2 + data[2,*, *,*]^2)
vertX =  REFORM(outverts[0,*],N_ELEMENTS(outverts)/3)
vertY =  REFORM(outverts[1,*],N_ELEMENTS(outverts)/3)
vertZ =  REFORM(outverts[2,*],N_ELEMENTS(outverts)/3)
vertcolors = BYTSCL(INTERPOLATE(magdata,vertX, vertY, vertZ))


if(KEYWORD_SET(tubes)) then begin    ;square profile for stream-tubes.
   width=0.2
   profile = [[-1,-1],[-1,1],[1,1],[1,-1],[-1,-1]]
   profile=fltarr(2,10)
   for i = 0,9 do begin
      profile(0,i) = cos(i*40*!dtor)
      profile(1,i) = sin(i*40*!dtor)
   endfor
   profile(1,*) = stepsize*profile(1,*)

   nverts = N_ELEMENTS(outverts)/3
   STREAMLINE, TEMPORARY(outverts),TEMPORARY(outconn), $
               outnormals*width,outverts,outconn, PROFILE=profile

   magdata = SQRT(data[0,*, *,*]^2 + data[1,*, *,*]^2 + data[2,*, *,*]^2)
   vertX =  REFORM(outverts[0,*],N_ELEMENTS(outverts)/3)
   vertY =  REFORM(outverts[1,*],N_ELEMENTS(outverts)/3)
   vertZ =  REFORM(outverts[2,*],N_ELEMENTS(outverts)/3)
   vertcolors = BYTSCL(INTERPOLATE(magdata,vertX, vertY, vertZ))


   ovx = reform(outverts(0,*))
   ovy = reform(outverts(1,*))
   ovz = reform(outverts(2,*))
   
;truncated float to interger (e.g. 2.56 = 2)
   fovx = fix(ovx)
   fovy = fix(ovy)
   fovz = fix(ovz)
   
;interpolate floating point indices to intermediate grid locations
;switch x and z directions
   px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
   py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
   pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 
   
   outverts(0,*) = px
   outverts(1,*) = py
   outverts(2,*) = pz
   
   
   oStreamlines = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=1.0)



endif


if(KEYWORD_SET(lines)) then begin $   ;square profile for stream-tubes.

   ovx = reform(outverts(0,*))
   ovy = reform(outverts(1,*))
   ovz = reform(outverts(2,*))
   
;truncated float to interger (e.g. 2.56 = 2)
   fovx = fix(ovx)
   fovy = fix(ovy)
   fovz = fix(ovz)
   
;interpolate floating point indices to intermediate grid locations
;switch x and z directions
   px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
   py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
   pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 
   
   outverts(0,*) = px
   outverts(1,*) = py
   outverts(2,*) = pz
   
   outverts=MESH_SMOOTH(outverts,outconn)

   oStreamlines = OBJ_NEW('IDLgrPolyline',outverts, $
                          POLYLINES=outconn) ;, $
;                       LABEL_OBJECTS=oLblSymbols, $
;                       LABEL_POLYLINES=lblPolys, $
;                       LABEL_OFFSETS=lblOffsets, $
;                       /LABEL_USE_VERTEX_COLOR, $
;                       /LABEL_NOGAPS)

endif

oPalette = OBJ_NEW('IDLgrPalette')
oPalette->LOADCT, ctbl
oStreamlines->SetProperty, PALETTE = oPalette, VERT_COLORS = vertcolors

;oModel->Add, oStreamlines 

uf2 = sqrt(uparr(*,*,*,0)^2 + uparr(*,*,*,1)^2 + uparr(*,*,*,2)^2)
oCont = OBJ_NEW('IDLgrContour',uf2(whx,why,10),geomx = x,geomy = y, geomz = z(10),$
                planar=1,fill=1,n_levels=255,palette=oPalette,shading=1,$
               alpha_channel=ach,c_color=bytscl(indgen(255)))
;oModel->Add,oCont

;ocb = OBJ_NEW('IDLgrColorbar',uf2(whx,why,40),Palette=oPalette)
;ocb = OBJ_NEW('IDLgrColorbar',parent=oCont,Palette=oPalette)
;oModel->Add,ocb


ntot = nparr
ntot = smooth(ntot,2)
isosurface,ntot(whx(0):whx(n_elements(whx)-1),$
                 why(0):why(n_elements(why)-1),$
                 whz(0):whz(n_elements(whz)-1))/1e15,den_val,$
           outverts,outconn

outverts=MESH_SMOOTH(outverts,outconn)

ovx = reform(outverts(0,*))
ovy = reform(outverts(1,*))
ovz = reform(outverts(2,*))

;truncated float to interger (e.g. 2.56 = 2)
fovx = fix(ovx)
fovy = fix(ovy)
fovz = fix(ovz)

px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

outverts(0,*) = px
outverts(1,*) = py
outverts(2,*) = pz

oden = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=0.7,palette=opalette,$
			color=!color.grey)

;oModel->Add, oden


ntot= nparr
;ti = (nfarr/ntot)*(1e6*pfarr/nfarr)/1.6e-19 $;
;       + ((nparr/ntot)*temp_p_arr>0)
;ti = 1e6*(pfarr/nfarr)/1.6e-19
;temp
ti = temp_p_arr>0
ti = smooth(ti,2)
isosurface,ti(whx(0):whx(n_elements(whx)-1),$
                 why(0):why(n_elements(why)-1),$
                 whz(0):whz(n_elements(whz)-1)),temp_val,$
           outverts,outconn


outverts=MESH_SMOOTH(outverts,outconn)

ovx = reform(outverts(0,*))
ovy = reform(outverts(1,*))
ovz = reform(outverts(2,*))

;truncated float to interger (e.g. 2.56 = 2)
fovx = fix(ovx)
fovy = fix(ovy)
fovz = fix(ovz)

px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

outverts(0,*) = px
outverts(1,*) = py
outverts(2,*) = pz

otemp = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=ach,palette=opalette,$
	  color=!color.blue)

;oModel->Add, otemp

@get_const
bx = reform(smooth(b1arr(*,*,*,0),2))*mp/q
;bx = abs(bx)
;bx

isosurface,(bx(whx(0):whx(n_elements(whx)-1),$
                 why(0):why(n_elements(why)-1),$
                 whz(0):whz(n_elements(whz)-1))),b_val,$
           outverts,outconn


outverts=MESH_SMOOTH(outverts,outconn)

ovx = reform(outverts(0,*))
ovy = reform(outverts(1,*))
ovz = reform(outverts(2,*))

;truncated float to interger (e.g. 2.56 = 2)
fovx = fix(ovx)
fovy = fix(ovy)
fovz = fix(ovz)

px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

outverts(0,*) = px
outverts(1,*) = py
outverts(2,*) = pz

obx = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=0.7,palette=opalette,$
			color=!color.green)


by = reform(smooth(b1arr(*,*,*,0),2))*mp/q
;bx = abs(bx)
;bx

isosurface,by(whx(0):whx(n_elements(whx)-1),$
                 why(0):why(n_elements(why)-1),$
                 whz(0):whz(n_elements(whz)-1)),-b_val,$
           outverts,outconn


outverts=MESH_SMOOTH(outverts,outconn)

ovx = reform(outverts(0,*))
ovy = reform(outverts(1,*))
ovz = reform(outverts(2,*))

;truncated float to interger (e.g. 2.56 = 2)
fovx = fix(ovx)
fovy = fix(ovy)
fovz = fix(ovz)

px = (ovx - fovx)*(x(fovx+1)-x(fovx)) + x(fovx) 
py = (ovy - fovy)*(y(fovy+1)-y(fovy)) + y(fovy) 
pz = (ovz - fovz)*(z(fovz+1)-z(fovz)) + z(fovz) 

outverts(0,*) = px
outverts(1,*) = py
outverts(2,*) = pz

oby = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                          SHADING = 1,alpha_channel=0.7,palette=opalette,$
			color=!color.red)

;oModel->Add, oby



;add Galileo trajectory

traj = fltarr(3,n_elements(x24))
for i = 0,n_elements(x24)-1 do begin
    traj(0,i) = x24(i)
    traj(1,i) = y24(i)    
    traj(2,i) = z24(i)
endfor

oI24 = OBJ_NEW('IDLgrPolyline',data=traj,color=!color.white,thick=4)


traj = fltarr(3,n_elements(x27))
for i = 0,n_elements(x27)-1 do begin
    traj(0,i) = x27(i)
    traj(1,i) = y27(i)    
    traj(2,i) = z27(i)
endfor

oI27 = OBJ_NEW('IDLgrPolyline',data=traj,color=!color.white,thick=4)


traj = fltarr(3,n_elements(x31))
for i = 0,n_elements(x31)-1 do begin
    traj(0,i) = x31(i)
    traj(1,i) = y31(i)    
    traj(2,i) = z31(i)
endfor

oI31 = OBJ_NEW('IDLgrPolyline',data=traj,color=!color.white,thick=4)

traj = fltarr(3,n_elements(xJ0))
for i = 0,n_elements(xJ0)-1 do begin
    traj(0,i) = xJ0(i)
    traj(1,i) = yJ0(i)    
    traj(2,i) = zJ0(i)
endfor

oJ0 = OBJ_NEW('IDLgrPolyline',data=traj,color=!color.white,thick=4)

;add all objects here---------------------------------------------

oModel->Add, xtit
oModel->Add,xax
oModel->Add,yax
oModel->Add,zax
oModel -> ADD, oPolygons  ;this is Io
oModel->Add, oStreamlines 
;oModel->Add, oI24
;oModel->Add, oI27
;oModel->Add, oI31
;oModel->Add, oJ0

;oModel -> ADD, oJupiter
;oModel->Add,oCont
;oModel->Add,ocb

;oModel->Add, otemp
oModel->Add, obx
oModel->Add, oby
oModel->Add, oden


;-----------------------------------------------------------------


oModel->Rotate, [1,0,0], -90
oModel->Rotate, [0,1,0], 30*9.5
oModel->Rotate, [1,0,0], 30

oLight = OBJ_NEW('IDLgrLight', loc=[100,-100,0], Inten=1.0, type=2)
oModel->Add,oLight

aLight = OBJ_NEW('IDLgrLight', Inten=1.0, type=0)
oModel->Add,aLight


view = obj_new('IDLgrView',color=[0,0,0],dimensions=[0,0],eye=9,$
 	viewplane_rect=10*[-3.5,-3.5,7.0,7.0],zclip = 10*[7,-1000])
view->add, oModel

;win = obj_new('IDLgrWindow',dimensions=[1200,1200], graphics_tree=view)

;win->draw, view
;win->GetProperty, image_data=im

;XINTERANIMATE, SET=[1200,1200,361], /SHOWLOAD


;for i = 0,360/4. do begin
;;for i = 0,n_elements(x31)-2 do begin

;   omodel->rotate,[0,1,0],1*4.

;;   dx = x31(i+1)-x31(i)
;;   dy = y31(i+1)-y31(i)
;;   dz = z31(i+1)-z31(i)
;;   oModel->translate,-dx,-dy,-dz
;;view->add, oModel
;   win->draw, view
;   win->GetProperty, image_data=im
;   xinteranimate, frame = i, image = im
   
;endfor
;xinteranimate,/keep_pixmaps

XOBJVIEW, oModel, scale=1.5,/BLOCK,xsize=1200,ysize=1200

OBJ_DESTROY, [oModel, oPalette]

end















