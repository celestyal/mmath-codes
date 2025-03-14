	PRO FXREAD, FILENAME, DATA, HEADER, P1, P2, P3, P4, P5,		$
		NANVALUE=NANVALUE, PROMPT=PROMPT, AVERAGE=AVERAGE,	$
		YSTEP=Y_STEP, NOSCALE=NOSCALE, NOUPDATE=NOUPDATE,	$
		ERRMSG=ERRMSG
;+
; Project     : SOHO - CDS
;
; Name        : 
;	FXREAD
; Purpose     : 
;	Read basic FITS files.
; Explanation : 
;	Read the primary array from a disk FITS file.  Optionally allows the
;	user to read in only a subarray and/or every Nth pixel.
; Use         : 
;	FXREAD, FILENAME, DATA  [, HEADER  [, I1, I2  [, J1, J2 ]]  [, STEP ]]
; Inputs      : 
;	FILENAME = String containing the name of the file to be read.
; Opt. Inputs : 
;	I1,I2	 = Data range to read in the first dimension.  If passed, then
;		   HEADER must also be passed.  If not passed, or set to 0,0
;		   then the entire range is read.
;	J1,J2	 = Data range to read in the second dimension.  If passed, then
;		   HEADER and I1,J2 must also be passed.  If not passed, or set
;		   to 0,0 then the entire range is read.
;	STEP	 = Step size to use in reading the data.  If passed, then
;		   HEADER must also be passed.  Default value is 1.  Ignored if
;		   less than 1.
; Outputs     : 
;	DATA	 = Data array to be read from the file.
; Opt. Outputs: 
;	HEADER	 = String array containing the header for the FITS file.
; Keywords    : 
;	NANVALUE = Value signalling data dropout.  All points corresponding to
;		   IEEE NaN (not-a-number) are set to this value.  Ignored
;		   unless DATA is of type float or double-precision.
;	PROMPT	 = If set, then the optional parameters are prompted for at the
;		   keyboard.
;	AVERAGE	 = If set, then the array size is reduced by averaging pixels
;		   together rather than by subselecting pixels.  Ignored unless
;		   STEP is nontrivial.  Note:  this is much slower.
;	YSTEP	 = If passed, then STEP is the step size in the 1st dimension,
;		   and YSTEP is the step size in the 2nd dimension.  Otherwise,
;		   STEP applies to both directions.
;	NOSCALE	 = If set, then the output data will not be scaled using the
;		   optional BSCALE and BZERO keywords in the FITS header.
;		   Default is to scale, if and only if BSCALE and BZERO are
;		   present and nontrivial.
;	NOUPDATE = If set, then the optional BSCALE and BZERO keywords in the
;		   optional HEADER array will not be changed.  The default is
;		   to reset these keywords to BSCALE=1, BZERO=0.  Ignored if
;		   NOSCALE is set.
;	ERRMSG   = If defined and passed, then any error messages will be
;		   returned to the user in this parameter rather than
;		   depending on the MESSAGE routine in IDL.  If no errors are
;		   encountered, then a null string is returned.  In order to
;		   use this feature, ERRMSG must be defined first, e.g.
;
;			ERRMSG = ''
;			FXREAD, ERRMSG=ERRMSG, ...
;			IF ERRMSG NE '' THEN ...
;
; Calls       : 
;	GET_DATE, IEEE_TO_HOST, FXADDPAR, FXHREAD, FXPAR, WHERENAN
; Common      : 
;	None.
; Restrictions: 
;	Groups are not supported.
;
;	The optional parameters I1, I2, and STEP only work with one or
;	two-dimensional arrays.  J1 and J2 only work with two-dimensional
;	arrays.
;
;	Use of the AVERAGE keyword is not compatible with arrays with missing
;	pixels.
;
; Side effects: 
;	If the keywords BSCALE and BZERO are present in the FITS header, and
;	have non-trivial values, then the returned array DATA is formed by the
;	equation
;
;			DATA = BSCALE*original + BZERO
;
;	However, this behavior can overridden by using the /NOSCALE keyword.
;
;	If the data is scaled, then the optional HEADER array is changed so
;	that BSCALE=1 and BZERO=0.  This is so that these scaling parameters
;	are not applied to the data a second time by another routine.  Also,
;	history records are added storing the original values of these
;	constants.  Note that only the returned array is modified--the header
;	in the FITS file itself is untouched.
;
;	If the /NOUPDATE keyword is set, however, then the BSCALE and BZERO
;	keywords are not changed.  It is then the user's responsibility to
;	ensure that these parameters are not reapplied to the data.  In
;	particular, these keywords should not be present in any header when
;	writing another FITS file, unless the user wants their values to be
;	applied when the file is read back in.  Otherwise, FITS readers will
;	read in the wrong values for the data array.
;
; Category    : 
;	Data Handling, I/O, FITS, Generic.
; Prev. Hist. : 
;	W. Thompson, May 1992, based in part on READFITS by W. Landsman, and
;			       STSUB by M. Greason and K. Venkatakrishna.
;	W. Thompson, Jun 1992, added code to interpret BSCALE and BZERO
;			       records, and added NOSCALE and NOUPDATE
;			       keywords.
;	W. Thompson, Aug 1992, changed to call FXHREAD, and to add history
;			       records for BZERO, BSCALE.
; Written     : 
;	William Thompson, GSFC, May 1992.
; Modified    : 
;	Version 1, William Thompson, GSFC, 12 April 1993.
;		Incorporated into CDS library.
;	Version 2, William Thompson, GSFC, 17 November 1993.
;		Corrected bug with AVERAGE keyword on non-IEEE compatible
;		machines.
;		Corrected bug with subsampling on VAX machines.
;	Version 3, William Thompson, GSFC, 31 May 1994
;		Added ERRMSG keyword.
;       Version 4, William Thompson, GSFC, 23 June 1994
;               Modified so that ERRMSG is not touched if not defined.
; Version     :
;       Version 4, 23 June 1994
;-
;
	ON_ERROR, 2
;
;  This parameter will be used later in conjunction with the average keyword.
;
	ALREADY_CONVERTED = 0
;
;  Parse the input parameters.
;
	CASE N_PARAMS() OF
		2:  BEGIN & I1=0  & I2=0  & J1=0  & J2=0  & STEP=1  & END
		3:  BEGIN & I1=0  & I2=0  & J1=0  & J2=0  & STEP=1  & END
		4:  BEGIN & I1=0  & I2=0  & J1=0  & J2=0  & STEP=P1 & END
		5:  BEGIN & I1=P1 & I2=P2 & J1=0  & J2=0  & STEP=1  & END
		6:  BEGIN & I1=P1 & I2=P2 & J1=0  & J2=0  & STEP=P3 & END
		7:  BEGIN & I1=P1 & I2=P2 & J1=P3 & J2=P4 & STEP=1  & END
		8:  BEGIN & I1=P1 & I2=P2 & J1=P3 & J2=P4 & STEP=P5 & END
		ELSE:  BEGIN
			MESSAGE = 'Syntax:  FXREAD, FILENAME, DATA ' +	$
				'[, HEADER [, I1, I2 [, J1, J2 ] [, STEP ]]'
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
			END
	ENDCASE
;
;  Get the UNIT number, and open the file.
;
	OPENR, UNIT, FILENAME, /BLOCK, /GET_LUN
;
;  Read in the FITS header.
;
	FXHREAD,UNIT,HEADER,STATUS
	IF STATUS NE 0 THEN BEGIN
		FREE_LUN,UNIT
		MESSAGE = 'Unable to read FITS header'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF
;
;  Extract the keywords BITPIX, NAXIS, NAXIS1, ...
;
	BITPIX = FXPAR(HEADER,'BITPIX')
	NAXIS = FXPAR(HEADER,'NAXIS')
	IF NAXIS EQ 0 THEN BEGIN
		FREE_LUN,UNIT
		MESSAGE = 'File does not contain a primary array'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	ENDIF
	DIMS = FXPAR(HEADER,'NAXIS*')
	N1 = DIMS(0)
	IF NAXIS EQ 2 THEN N2 = DIMS(1) ELSE N2 = 1
;
;  Determine the array type from the keyword BITPIX.
;
	CASE BITPIX OF
		  8:	IDLTYPE = 1	; Byte
		 16:	IDLTYPE = 2	; Integer*2
		 32:	IDLTYPE = 3	; Integer*4
		-32:	IDLTYPE = 4	; Real*4
		-64:	IDLTYPE = 5	; Real*8
	ENDCASE
;
;  Set the default values for the optional parameters.
;
	IF (I1 EQ 0) AND (I2 EQ 0) THEN I2 = N1-1
	IF (J1 EQ 0) AND (J2 EQ 0) THEN J2 = N2-1
;
;  If the prompt keyword was set, the prompt for the parameters.
;
	IF KEYWORD_SET(PROMPT) THEN BEGIN
		ANSWER = ''
		READ,'Enter lower limit for X ['+STRTRIM(I1,2)+']: ', ANSWER
		IF ANSWER NE '' THEN I1 = (ANSWER)
;
		ANSWER = ''
		READ,'Enter upper limit for X ['+STRTRIM(I2,2)+']: ', ANSWER
		IF ANSWER NE '' THEN I2 = LONG(ANSWER)
;
		ANSWER = ''
		READ,'Enter lower limit for Y ['+STRTRIM(J1,2)+']: ', ANSWER
		IF ANSWER NE '' THEN J1 = LONG(ANSWER)
;
		ANSWER = ''
		READ,'Enter upper limit for Y ['+STRTRIM(J2,2)+']: ', ANSWER
		IF ANSWER NE '' THEN J2 = LONG(ANSWER)
;
		ANSWER = ''
		READ,'Enter step size ['+STRTRIM(STEP,2)+']: ', ANSWER
		IF ANSWER NE '' THEN STEP = LONG(ANSWER)
	ENDIF
;
;  Differentiate between XSTEP and YSTEP.
;
	XSTEP = STEP > 1
	IF N_ELEMENTS(Y_STEP) EQ 1 THEN YSTEP = Y_STEP ELSE YSTEP = XSTEP
;
;  If any of the optional parameters were passed, then update the dimensions
;  accordingly.  First check I1 and I2.
;
	IF (I1 NE 0) OR (I2 NE N1-1) THEN BEGIN
		IF NAXIS GT 2 THEN BEGIN
			FREE_LUN,UNIT
			MESSAGE = 'Range parameters can only be set for ' + $
				'one or two-dimensional arrays'
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		ENDIF
		IF (MIN([I1,I2]) LT 0) OR (MAX([I1,I2]) GE DIMS(0)) THEN BEGIN
			FREE_LUN,UNIT
			MESSAGE = 'I1,I2 must be in the range 0 to ' +	$
				STRTRIM(DIMS(0)-1,2)
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		END ELSE IF I1 GE I2 THEN BEGIN
			MESSAGE = 'I2 must be greater than I1'
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		ENDIF
		DIMS(0) = I2 - I1 + 1
	ENDIF
;
;  Next, check J1 and J2.
;
	IF (J1 NE 0) OR (J2 NE N2-1) THEN BEGIN
		IF NAXIS NE 2 THEN BEGIN
			FREE_LUN,UNIT
			MESSAGE = 'J1, J2 can only be set for ' +	$
				'two-dimensional arrays'
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		ENDIF
		IF (MIN([J1,J2]) LT 0) OR (MAX([J1,J2]) GE DIMS(1)) THEN BEGIN
			FREE_LUN,UNIT
			MESSAGE = 'J1,J2 must be in the range 0 to ' +	$
				STRTRIM(DIMS(1)-1,2)
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		END ELSE IF J1 GE J2 THEN BEGIN
			MESSAGE = 'J2 must be greater than J1'
			IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
				ERRMSG = MESSAGE
				RETURN
			END ELSE MESSAGE, MESSAGE
		ENDIF
		DIMS(1) = J2 - J1 + 1
	ENDIF
;
;  Next, check XSTEP.  Note that the dimensions of the final result are
;  somewhat differ depending on whether the keyword AVERAGE is set or not.
;
	IF XSTEP GT 1 THEN BEGIN
	    IF NAXIS GT 2 THEN BEGIN
	        FREE_LUN,UNIT
	        MESSAGE = 'STEP can only be set for one or ' +	$
	            'two-dimensional arrays'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	    END ELSE IF XSTEP NE LONG(XSTEP) THEN BEGIN
	        FREE_LUN,UNIT
	        MESSAGE = 'STEP must be an integer value'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	    END ELSE IF KEYWORD_SET(AVERAGE) THEN BEGIN
	        DIMS(0) = DIMS(0) / LONG(XSTEP)
	    END ELSE BEGIN
	        DIMS(0) = LONG(DIMS(0) + XSTEP - 1) / LONG(XSTEP)
	        INDEX = LINDGEN(DIMS(0))*XSTEP
	    ENDELSE
	ENDIF
;
;  Finally, check YSTEP.  This parameter is ignored for anything other than
;  two-dimensional arrays.
;
	IF (NAXIS EQ 2) AND (YSTEP GT 1) THEN BEGIN
	    IF YSTEP NE LONG(YSTEP) THEN BEGIN
	        FREE_LUN,UNIT
	        MESSAGE = 'YSTEP must be an integer value'
		IF N_ELEMENTS(ERRMSG) NE 0 THEN BEGIN
			ERRMSG = MESSAGE
			RETURN
		END ELSE MESSAGE, MESSAGE
	    END ELSE IF KEYWORD_SET(AVERAGE) THEN BEGIN
	        DIMS(1) = DIMS(1) / LONG(YSTEP)
	    END ELSE BEGIN
	        DIMS(1) = LONG(DIMS(1)+YSTEP-1) / LONG(YSTEP)
	    ENDELSE
	END ELSE YSTEP = 1
;
;  Make the array.
;
	DATA = MAKE_ARRAY(DIMENSION=DIMS,TYPE=IDLTYPE,/NOZERO)
;
;  Find the start of the data to be read in.
;
	POINT_LUN,-UNIT,OFFSET		;Current position
	DELTA = N1*ABS(BITPIX)/8
	IF J1 NE 0 THEN BEGIN
		OFFSET = OFFSET + J1*DELTA
		POINT_LUN,UNIT,OFFSET
	ENDIF
;
;  If the I range, XSTEP or YSTEP is non-trivial, then read in the file line by
;  line.  If pixel averaging, then read in YSTEP lines.
;
	IF (DIMS(0) NE N1) OR (XSTEP GT 1) OR (YSTEP GT 1) THEN BEGIN
	    IF NAXIS EQ 1 THEN NJ = 1 ELSE NJ = DIMS(1)
	    FOR J = 0,NJ-1 DO BEGIN
	        IF YSTEP GT 1 THEN POINT_LUN,UNIT,OFFSET+J*YSTEP*DELTA
	        IF (YSTEP GT 1) AND KEYWORD_SET(AVERAGE) AND (NAXIS EQ 2) $
	            THEN LINE = MAKE_ARRAY(N1,YSTEP,TYPE=IDLTYPE,/NOZERO) $
	            ELSE LINE = MAKE_ARRAY(N1,TYPE=IDLTYPE,/NOZERO)
	        READU,UNIT,LINE
;
;  If I1,I2 do not match the array size, then extract the relevant subarray.
;
	        IF (I1 NE 0) OR (I2 NE N1-1) THEN LINE = LINE(I1:I2,*)
;
;  Suppose that the step size is non-trivial.  If AVERAGE was set, then convert
;  to the host format, and use REBIN to average the data.  (Note that missing
;  pixels are not correctly handled in this case.)  Otherwise, select out the
;  relevant portion of the data.
;
	        IF (XSTEP GT 1) OR (YSTEP GT 1) THEN BEGIN
	            IF KEYWORD_SET(AVERAGE) THEN BEGIN
			IEEE_TO_HOST, LINE
			ALREADY_CONVERTED = 1
	                IF NAXIS EQ 1 THEN BEGIN
	                    DATA(0,J) = REBIN(LINE(0:XSTEP*DIMS(0))-1,DIMS(0))
	                END ELSE BEGIN
	                    DATA(0,J) = REBIN(LINE(0:XSTEP*DIMS(0)-1,*),DIMS(0),1)
	                ENDELSE
;
;  It seems that VAXes don't like to even subsample floating point data that
;  isn't yet in local format.  Therefore it is necessary to convert beforehand.
;  To be compatible with the rest of the routine, the data is then converted
;  back again.
;
	            END ELSE IF (!VERSION.ARCH EQ 'vax') AND (IDLTYPE GE 4) $
			    THEN BEGIN
			TEST = LINE
			IEEE_TO_HOST, TEST
			TEST = TEST(INDEX)
			HOST_TO_IEEE, TEST
			DATA(0,J) = TEST
		    END ELSE DATA(0,J) = LINE(INDEX)
;
;  Otherwise, if the step size is trivial, then simply store the line in the
;  data array.
;
	        END ELSE BEGIN
	            DATA(0,J) = LINE
	        ENDELSE
	    ENDFOR
;
;  Otherwise, if the file doesn't have to be read in line by line, then just
;  read the data array.
;
	END ELSE READU,UNIT,DATA
;
;  Convert the data from IEEE to host format, keeping track of any IEEE NaN
;  values.  Don't do this if the conversion has already taken place.
;
	IF NOT ALREADY_CONVERTED THEN BEGIN
		IF (N_ELEMENTS(NANVALUE) EQ 1) AND (IDLTYPE GE 4) AND	$
			(IDLTYPE LE 6) THEN W = WHERENAN(DATA,COUNT) ELSE $
			COUNT = 0
		IEEE_TO_HOST,DATA
	END ELSE COUNT = 0
;
;  If the parameters BZERO and BSCALE are non-trivial, then adjust the array by
;  these values.
;
	IF NOT KEYWORD_SET(NOSCALE) THEN BEGIN
		BZERO  = FXPAR(HEADER,'BZERO')
		BSCALE = FXPAR(HEADER,'BSCALE')
		GET_DATE,DTE
		IF (BSCALE NE 0) AND (BSCALE NE 1) THEN BEGIN
			DATA = BSCALE*DATA
			IF NOT KEYWORD_SET(NOUPDATE) THEN BEGIN
				FXADDPAR,HEADER,'BSCALE',1.
				FXADDPAR,HEADER,'HISTORY',DTE +		$
					' applied BSCALE = '+ STRTRIM(BSCALE,2)
			ENDIF
		ENDIF
		IF BZERO NE 0 THEN BEGIN
			DATA = DATA + BZERO
			IF NOT KEYWORD_SET(NOUPDATE) THEN BEGIN
				FXADDPAR,HEADER,'BZERO',0.
				FXADDPAR,HEADER,'HISTORY',DTE +		$
					' applied BZERO = '+ STRTRIM(BZERO,2)
			ENDIF
		ENDIF
	ENDIF
;
;  Store NANVALUE everywhere where the data corresponded to IEE NaN.
;
	IF COUNT GT 0 THEN DATA(W) = NANVALUE
;
;  Close the file and return.
;
	FREE_LUN, UNIT
	IF N_ELEMENTS(ERRMSG) NE 0 THEN ERRMSG = ''
	RETURN
	END
