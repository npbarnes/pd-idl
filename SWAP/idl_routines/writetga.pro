PRO WRITETGA, im, format, filename
;+
; NAME: 
;     WRITETGA
;
; PURPOSE: 
;     Write a file in TGA format for 8bit B&W
;          capture board - 
;
; CATEGORY: 
;     Input/Output - special format
;
; CALLING SEQUENCE: 
;     WRITETGA, IM, format, [filename] 
;
; INPUTS:
;     format - gives the TGA format for the image...
;	go with what Hans has set up for the TARGA+ 
;	STAR/STEREO system...
;     im - image to save
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
;      none
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
   PRINT,'	PRO WRITETGA
   PRINT,'	   Reads 16 bit color (Red,Green,Blue all 5 bit) or '
   PRINT,'	          8 bit black & white TGA images.'
   PRINT,'	CALLING SEQUENCE:'
   PRINT,'	   READ_TGA, IMAGE, FORMAT, [FILENAME]'
   PRINT,'	WHERE - '
   PRINT,'	   IMAGE......array to be saved'
   PRINT,'	   FORMAT......is the TGA format to save: only '
   PRINT,'	      format = 3 supported (8bit B&W no colormap)'
   PRINT,'	   FILENAME....is the optional filename, if it is'
   PRINT,'	      not passed, the routine will prompt for the name'
   PRINT,' '
   GOTO, FINISH
ENDIF

SIM = size(im)
IF (SIM(3) NE 1) THEN BEGIN
    PRINT, '(writetga) Converting IMAGE to BYTE array'
    im = BYTSCL(im,top=31)
ENDIF

IF (N_PARAMS() LT 3) THEN BEGIN	; Must prompt for TGAfile...
    TGAfile = ''
    READ,' TGA file to open : ', TGAfile
ENDIF ELSE BEGIN		; Passed the filename...
    TGAfile = filename
ENDELSE

OPENW,11, TGAfile, error = ioerr
IF IOERR NE 0 THEN BEGIN
    PRINT, '%%% Error opening file - ',+TGAfile +' %%%'
    CLOSE,11
    GOTO, FINISH
ENDIF
text1 = '  256x243 5 bit B/W  -- Output by WRITETGA.PRO'
im_ident = [BYTE(text1),bytarr(192)]
n_char = BYTE(n_elements(im_ident))

col_map_type = 0b	&  image_type = 3b
col_map_orig = 0  	&  col_map_len = 0
col_map_entry = 0b
x_orig = 0  		&  y_orig = 0
im_width = FIX(SIM(1)) 	&  im_height = FIX(SIM(2))
im_pix_size = 8b	&  overlay = 0b
;------------------------- Write the Header -------------------------
WRITEU,11,n_char
WRITEU,11,col_map_type
WRITEU,11,image_type
WRITEU,11,col_map_orig
WRITEU,11,col_map_len
WRITEU,11,col_map_entry
WRITEU,11,x_orig
WRITEU,11,y_orig
WRITEU,11,im_width
WRITEU,11,im_height
WRITEU,11,im_pix_size
WRITEU,11,overlay
WRITEU,11,im_ident
;-----------------------------------------------------------------
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

IF n_char GT 0 THEN BEGIN
    PRINT,' Image Identifier String --"' + STRING(im_ident) + '"'
    PRINT,'****************************************'
    PRINT,' '
ENDIF ELSE BEGIN
    PRINT,'****************************************'
    PRINT,' '
ENDELSE

WRITEU,11,im

FINISH:
CLOSE,11
END
