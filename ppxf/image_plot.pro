;----------------------------------------------------------------------
PRO image_plot, z, x, y, WINDOW_SCALE=window_scale, ASPECT=aspect, $
        INTERP=interp, _EXTRA=extra 
compile_opt idl2
;+
;
; NAME:
;   IMAGE_PLOT
;
; PURPOSE:
;   Plot an image with axis.
;
; CATEGORY:
;   General graphics.
;
; CALLING SEQUENCE:
;   IMAGE_PLOT, Z, X, Y
;
; INPUTS:
;   Z:  The two-dimensional array to display.
;   X:  x scale (the syntax is the same as in CONTOUR)
;   Y:  y scale
;
; KEYWORD PARAMETERS:
; WINDOW_SCALE: Set this keyword to scale the window size to the image size.
;       Otherwise, the image size is scaled to the window size.
;       This keyword is ignored when outputting to devices with
;       scalable pixels (e.g., PostScript).
;
;   ASPECT: Set this keyword to retain the image's aspect ratio.
;       Square pixels are assumed.  If WINDOW_SCALE is set, the
;       aspect ratio is automatically retained.
;
;   INTERP: If this keyword is set, bilinear interpolation is used if
;       the image is resized.
;
;   _EXTRA: also accept all options accepted by PLOT
;
; OUTPUTS:
;   No explicit outputs.
;
; COMMON BLOCKS:
;   None.
;
; SIDE EFFECTS:
;   The currently selected display is affected.
;
; RESTRICTIONS:
;   None.
;
; PROCEDURE:
;   If the device has scalable pixels, then the image is written over
;   the plot window.
;
; MODIFICATION HISTORY:
;   V1.0: DMS, May, 1988.
;   V2.0: Michele Cappellari, May 1998:
;       (Adapted from IMAGE_CONT to use PLOT instead of CONTOUR)
;   V2.01: MC, May 2000: added _EXTRA feature
;   V2.02: Allow for optional (X,Y) vectors. MC, Leiden, 6 June 2003
;   V2.03: Included _EXTRA keyword also in the first call to PLOT.
;       Thank you to Jesus Falcon-Barrorso for pointing out this bug.
;       MC, Leiden, 20 July 2006
;   V2.04: Added one unit to the plotted X and Y coordinates, when no axses
;       are provided in input. MC, Oxford, 20 January 2007
;   V2.05: Removed _EXTRA keyword introduced in V2.03. MC, Oxford, 6 February 2007
;   V2.06: Use Coyote Graphics. MC, Oxford, 11 October 2011
;-

ON_ERROR,2                  ;Return to caller if an error occurs
sz = SIZE(z)                ;Size of image
IF sz[0] LT 2 THEN MESSAGE, 'Parameter not 2D'
IF N_PARAMS() EQ 1 THEN BEGIN
    x = indgen(sz[1]+1)
    y = indgen(sz[2]+1)
ENDIF

;set window used by plot
cgPLOT, [[0,0],[1,1]], /NODATA, XSTYLE=4, YSTYLE=4;, _EXTRA=extra

px = !X.WINDOW * !D.X_VSIZE ;Get size of window in device units
py = !Y.WINDOW * !D.Y_VSIZE
swx = px[1]-px[0]           ;Size in x in device units
swy = py[1]-py[0]           ;Size in Y
six = FLOAT(sz[1])          ;Image sizes
siy = FLOAT(sz[2])
aspi = six / siy            ;Image aspect ratio
aspw = swx / swy            ;Window aspect ratio
f = aspi / aspw             ;Ratio of aspect ratios

IF (!D.FLAGS AND 1) NE 0 THEN BEGIN                 ;Scalable pixels?
    IF KEYWORD_SET(aspect) THEN BEGIN               ;Retain aspect ratio?
        IF f GE 1.0 THEN swy=swy/f ELSE swx=swx*f   ;Adjust window size
    ENDIF
    TVSCL, z, px[0], py[0], XSIZE=swx, YSIZE=swy, /DEVICE
ENDIF ELSE BEGIN                                    ;Not scalable pixels
    IF KEYWORD_SET(window_scale) THEN BEGIN         ;Scale window to image?
        TVSCL, z, px[0], py[0]                      ;Output image
        swx = six                                   ;Set window size from image
        swy = siy
    ENDIF ELSE BEGIN                                ;Scale window
        IF KEYWORD_SET(aspect) THEN BEGIN
            IF f GE 1.0 THEN swy=swy/f ELSE swx=swx*f
        ENDIF                                       ;aspect
        TV, POLY_2D(BYTSCL(z), $                    ;Have to resample image
            [[0,0],[six/swx,0]],[[0,siy/swy],[0,0]], $
            KEYWORD_SET(interp), swx, swy), px[0], py[0]
    ENDELSE                                         ;window_scale
ENDELSE                                             ;scalable pixels

cgPLOT, [MIN(x),MAX(x)], [MIN(y),MAX(y)], /NODATA, /NOERASE, /DEVICE, $
       POS=[px[0],py[0],px[0]+swx,py[0]+swy], /XSTYLE, /YSTYLE, _EXTRA=extra

RETURN
END
;----------------------------------------------------------------------
