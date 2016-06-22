pro convert_date_time_to_year_arr, date, time,d
  sz=n_elements(date)
  d={month:0L, day:0L, year:0L,hour:0L, min:0L, sec:0D, daynumber:0L, fdoy:0D, tyrr:0D, date:' ', time:'', fp_year:0D, fp_hrs:0.D}
  d=replicate(d,sz)
  d.date=strcompress(date,/remove_all)
  d.month=long(strmid(d.date, 0,2))
  d.day=long(strmid(d.date,3,2))
  d.year=long(strmid(d.date,6,4))
  d.time=strcompress(time, /remove_all)
  d.hour=long(strmid(d.time,0,2))
  d.min=long(strmid(d.time,3,2))

  if strlen(time(0)) gt 8 then begin 
    d.sec=double(strmid(d.time,6,10))
  endif 
  if strlen(time(0)) gt 4 and  strlen(time(0)) le  8 then begin 
     d.sec=double(strmid(d.time,6,2))
  endif
  if strlen(time(0)) le 4 then begin 
    d.sec=0
  endif
  d.fdoy=julday(d.month,d.day,d.year, d.hour, d.min, d.sec)-julday(1,1,d.year, 0, 0, 0)+1.D
  d.daynumber=long(d.fdoy)
  ndays=dblarr(sz)
  ndays(*)=365.D
  imod=where (d.year mod 4 eq 0)
  if imod(0) gt -1 then ndays(imod)=366.D
  d.fp_year=double(d.year)+(d.fdoy-1.D )/ndays
  d.tyrr=d.fp_year
  d.fp_hrs=d.hour+d.min/60.D +d.sec/3600.D
 
end
