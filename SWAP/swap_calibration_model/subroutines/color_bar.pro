pro color_bar, lim, ytit, pos, ct = ct
;I have finaly gotten smart and wrote a stupid color bar procedure.
;This will place a color bar that ranges the full rang of the
;currently loaded colar table.  
;lim = vector of size 2 (lower bound, upper bound)
;title = the title to put along the y axis
;pos = the 4 element vector ofthe location to place the color bar


if(n_elements(ct) EQ 0) then ct = 39
bar = fltarr(2, 256)
bar[0,*] = findgen(256)
bar[1,*] = findgen(256)

blank = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
loadct, ct, /silent
tvimage, bar, position = pos
loadct, 39, /silent
plot, [0,1], lim, /nodata, /xstyle, /ystyle, $
  position = pos, xtickname = blank, $
  xminor = 1, xticks = 1, /noerase, ytickname = blank
axis, yaxis = 1, ytitle = ytit, /ystyle

end
