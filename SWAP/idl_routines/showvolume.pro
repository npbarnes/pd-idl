Pro SHOWVOLUME, vol, thresh, LOW = low
;Display the contour surface of a volume.

s = SIZE(vol)	;Get the dimensions of the volume.
print,s(0),s(1),s(2),s(3)
IF s(0) NE 3 THEN print,'Error... must be a 3D array.'

SCALE3, XRANGE=[0, S(1)], YRANGE=[0, S(2)], ZRANGE=[0, S(3)]

;	Use SCALE3 to establish the 3D transforma-
;tion and coordinate ranges.



IF N_ELEMENTS(low) EQ 0 THEN low = 0
;Default = view high side of contour surface.

SHADE_VOLUME, vol, thresh, v, p, LOW = low
;Produce vertices and polygons.

TV, POLYSHADE(v,p,/T3D)	;Produce image of surface and display.

END
