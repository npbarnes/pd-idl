PRO READ_TGA, im, filename
;+
; NAME: 
;     READ_TGA
;
; PURPOSE: 
;     Read a file made with a Targa, or Targa+ image
;          capture board - 
;
; CATEGORY: 
;     Input/Output - special format
;
; CALLING SEQUENCE: 
;     READ_TGA, IM, [filename] 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
;      filename - name of TGA file (with all paths & extensions)
;                 if not supplied - routine will prompt for name.
;	
; KEYWORD PARAMETERS:
;
;      none
;
; OUTPUTS:
;
;      IM...........BYTE Array with dimensions of image
;
; OPTIONAL OUTPUTS:
;
;      on screen - gives a list of pertinent paramters
;
; COMMON BLOCKS:
;
;      none
;
; SIDE EFFECTS:
;
;      none
;
; RESTRICTIONS:
;
;      Only two of several modes are supported:
;           8 bit black & white
;          16 bit true color - 5 bit red, green & blue averaged
;                              together, no color map
;
; PROCEDURE:
;
;      Read the header - determine if it is a supported format -
;      find the image size - read it - return the image
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;      9/20/93 Added 8 bit, and spruced things up a bit
;      D Hampton - Geophysical Institute, University of Alaska
;-

on_error, 2                       ; return to caller

IF (N_PARAMS() EQ 0) THEN BEGIN
   PRINT,' '
   PRINT,'	PRO READ_TGA, IMAGE, [FILENAME]'
   PRINT,'	   Reads 16 bit color (Red,Green,Blue all 5 bit) or '
   PRINT,'	          8 bit black & white TGA images.'
   PRINT,'	CALLING SEQUENCE:'
   PRINT,'	   READ_TGA, IMAGE, [FILENAME]'
   PRINT,'	WHERE - '
   PRINT,'	   IMAGE.......is the returned image - is a BYTE array'
   PRINT,'	      with dimensions of the stored image in TGA file.'
   PRINT,'	   FILENAME....is the optional filename, if it is'
   PRINT,'	      not passed, the routine will prompt for the name'
   PRINT,' '
   GOTO, FINISH
ENDIF

IF (N_PARAMS() LT 2) THEN BEGIN	; Must prompt for TGAfile...
    TGAfile = ''
    READ,' TGA file to open : ', TGAfile
ENDIF ELSE BEGIN		; Passed the filename...
    TGAfile = filename
ENDELSE

OPENR,11, TGAfile, error = ioerr
IF IOERR NE 0 THEN BEGIN
    PRINT, '%%% Error opening file - ',+TGAfile +' %%%'
    CLOSE,11
    GOTO, FINISH
ENDIF
n_char = 0b
head = bytarr(17)
READU,11,n_char
READU,11,head

col_map_type = byte(head,0)	&  image_type = byte(head,1)
col_map_orig = fix(head,2)  	&  col_map_len = fix(head,4)
col_map_entry = byte(head,6)
x_orig = fix(head,7)  		&  y_orig = fix(head,9)
im_width = fix(head,11)  	&  im_height = fix(head,13)
im_pix_size = byte(head,15)

PRINT,'  '
PRINT,'****************************************'
PRINT,'      Header Data for ' + TGAfile
PRINT,' Color Map Type : ', col_map_type
PRINT,' Image Type : ', image_type
PRINT,' Color Map Origin : ', col_map_orig
PRINT,' Color Map Length : ', col_map_len
PRINT,' Color Map Entry : ', col_map_entry
PRINT,' X Origin : ', x_orig
PRINT,' Y Origin : ', y_orig
PRINT,' Image Width : ', im_width
PRINT,' Image Height : ', im_height
PRINT,' Image Pixel Size : ', im_pix_size
PRINT,' Overlay info :', byte(head,16)
IF n_char GT 0 THEN BEGIN
    im_ident = bytarr(n_char)
    READU,11, im_ident
    PRINT,' n_char : ',n_char
    PRINT,' Image Identifier String --"' + STRING(im_ident) + '"'
    PRINT,'****************************************'
    PRINT,' '
ENDIF ELSE BEGIN
    PRINT,'****************************************'
    PRINT,' '
ENDELSE

IF (col_map_len NE 0) THEN BEGIN	; if color map - read...
    CASE col_map_entry OF
	8b: cmap = BYTARR(col_map_len)
	16b: cmap = INTARR(col_map_len)
	32b: cmap = FLTARR(col_map_len)
    ELSE: BEGIN
	PRINT,' Strange col_map_entry', col_map_entry
	nmaps = (col_map_entry/8)*col_map_len
	PRINT,'nmaps = ', nmaps
	cmap = bytarr(nmaps)
    	  END
    ENDCASE
    READU,11, cmap
ENDIF
;----------------- Now read images...
IF (im_pix_size EQ 16b) THEN BEGIN
				 	; assume no color table
    image = intarr(im_width,im_height)
    readu,11,image

    CLOSE,11
    IF (image_type NE 2) THEN BEGIN  	; and 5 bit rgb
	rd= (image and 31744)/512
	gn= (image and 992)/16
	bl= (image and 31)*2
	im = BYTE((rd+gn+bl+2)/3)
    ENDIF ELSE BEGIN
	im = image
    ENDELSE
    GOTO, FINISH
ENDIF
					; 8 bit black and white
IF (im_pix_size EQ 8b) THEN BEGIN ; AND (image_type EQ 3) THEN BEGIN	
    im = bytarr(im_width,im_height)
    readu,11,im
    CLOSE,11
    GOTO, FINISH
ENDIF

PRINT, 'Pixel size not supported ... yet'

FINISH:
CLOSE,11
END
