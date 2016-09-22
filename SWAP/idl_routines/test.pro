 ; Create some data.
          vol = RANDOMU(s, 40, 40, 40)
          FOR i=0, 10 DO vol = SMOOTH(vol, 3)
          vol = BYTSCL(vol(3:37, 3:37, 3:37))
          opaque = RANDOMU(s, 40, 40, 40)
          FOR i=0, 10 DO opaque = SMOOTH(opaque, 3)
          opaque = BYTSCL(opaque(3:37, 3:37, 3:37), TOP=25B)

       ; Set up the view.
          T3D, /RESET
          T3D, TRANSLATE=[-0.5, -0.5, -0.5]
          T3D, SCALE=[0.7, 0.7, 0.7]
          T3D, ROTATE=[30, -30, 60]
          T3D, TRANSLATE=[0.5, 0.5, 0.5]

       ; Generate and display the image.
          img = PROJECT_VOL(vol, 64, 64, 64, DEPTH_Q=0.7, $
                OPAQUE=opaque, TRANS=(!P.T))
          TVSCL, img

end
