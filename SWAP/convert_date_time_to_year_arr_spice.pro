pro convert_date_time_to_year_arr_spice, date, time,d
  sz=n_elements(date)
  d={month:0L, day:0L, year:0L,hour:0L, min:0L, sec:0D, msec:0D, daynumber:0L, fdoy:0D, tyrr:0D, date:' ', time:'', fp_year:0D, fp_hrs:0.D}
  d=replicate(d,sz)
  d.date=strmid(strcompress(date,/remove_all),0,10)
  d.month=long(strmid(d.date, 5,2))
  d.day=long(strmid(d.date,8,2))
  d.year=long(strmid(d.date,0,4))
  ;print, d.year
  
  d.time=strcompress(strmid(date,11,14), /remove_all)
  d.hour=long(strmid(d.time,0,2))
  d.min=long(strmid(d.time,3,2))
  if strlen(d(0).time) gt 4 then begin 
    d.sec=double(strmid(d.time,6,2))
  endif else begin
    d.sec=0
  endelse
  d.msec=double(strmid(d.time,8,6))
  ;print, strmid(d.time,8,6)
  d.fdoy=julday(d.month,d.day,d.year, d.hour, d.min, d.sec)-julday(1,1,d.year, 0, 0, 0)+1.D
  d.fdoy=d.fdoy+d.msec/24./3600.
  d.daynumber=long(d.fdoy)
  ndays=dblarr(sz)
  ndays(*)=365.D
  imod=where (d.year mod 4 eq 0)
  if imod(0) gt -1 then ndays(imod)=366.D
  d.fp_year=double(d.year)+(d.fdoy-1.D )/ndays
  d.tyrr=d.fp_year
  d.fp_hrs=d.hour+d.min/60.D +(d.sec+d.msec)/3600.D
 
end
