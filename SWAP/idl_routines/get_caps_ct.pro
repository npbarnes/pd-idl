;----------------------------------------------------------------
pro get_caps_ct
;----------------------------------------------------------------

openr,1,'/Users/delamere/idl_routines/caps_ct.dat'

ct = {color_tab, r: 0.0, g: 0.0, b:0.0}

readf,1,ct

r = ct.r
g = ct.g
b = ct.b
while (not(eof(1))) do begin
readf,1,ct

r = [r,ct.r]
g = [g,ct.g]
b = [b,ct.b]
endwhile

r = bytscl(r)
g = bytscl(g)
b = bytscl(b)

r = congrid(r,256)
g = congrid(g,256)
b = congrid(b,256)

r(255) = 255
g(255) = 255
b(255) = 255

d = dist(100)
tvlct,r,g,b

;device,decompose=0
;tvscl,d

close,1

return
end
;----------------------------------------------------------------

