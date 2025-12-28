;+
; :Description:
;   Reads metadata about each Carrington rotation from the file supplied with
;   NSO/Kitt Peak synoptic maps.
;
; :Arguments:
;   rot: out, Array
;     Carrington rotation number.
;   filepath: in, String
;     Path to the file containing the Carrington rotation metadata.
;   date: out, optional, Array
;     Year (in decimal).
;   b0: out, optional, Array
;     B-angle
;   rd: out, optional, Array
;     Sun's semi-diameter in arc-minutes
;   mm: out, optional, Array
;     Month (as an integer)
;   dd: out, optional, Array
;     Day of month, given by day + (hour+min/60.+sec/3600.)/24.
;   yy: out, optional, Array
;     Year (integer)
;   jul: out, optional, Array
;     Julian date
;
;-
pro carrdata, rot, filepath, date=date, b0=b0, rd=rd, mm=mm, day=day, yy=yy, jul=jul
  ;+ Prepare variables to read in observation metadata
  a = 1.0 & b = 1 & c = 1.0 & d = 1.0 & e = 1 & f = 1.0 & g = 1 & h = 1.0
  
  ; Initialise arrays
  if arg_present(date) then date = []
  if arg_present(b0) then b0 = []
  if arg_present(rd) then rd = []
  if arg_present(mm) then mm = []
  if arg_present(dd) then dd = []
  if arg_present(yy) then yy = []
  if arg_present(jul) then jul = []

  get_lun, z
  openr, z, filepath
  
  ;+ Read in metadata line-by-line
  while not EOF(z) do begin
    readf, z, a, b, c, d, e, f, g, h
    if (where(b eq rot) gt -1) then begin
      if arg_present(date) then date = [date, a]
      if arg_present(b0) then b0 = [b0, c]
      if arg_present(rd) then rd = [rd, d]
      if arg_present(mm) then mm = [mm, e]
      if arg_present(dd) then dd = [dd, f]
      if arg_present(yy) then yy = [yy, g]
      if arg_present(jul) then jul = [jul, h]
    endif
  endwhile

  close, z
  free_lun, z
end
