function csv_template_type, csvfile, fieldtypes
  compile_opt idl2

  ; Read the header and the first line of data from the CSV file
  openr, lun, csvfile, /GET_LUN
  header = '' & readf, lun, header
  record = '' & readf, lun, record
  free_lun, lun

  ; Define regex patterns for acceptable variable names and sizes
  name_regex = '^[a-z][a-z0-9_]*(\[[1-9][0-9]*])?$'
  size_regex = '\[[0-9]*]'

  ; Split the headers based on commas and remove all whitespace
  tokens = strcompress(strsplit(header, ',', /EXTRACT), /REMOVE_ALL)
  names = strarr(n_elements(tokens))
  sizes = intarr(n_elements(tokens))

  ; Determine the size of each token; if not specified the size is 1
  for i = 0, n_elements(tokens) - 1 do begin
    if not stregex(tokens[i], name_regex, /FOLD_CASE, /BOOLEAN) then begin
      message, 'Field ' + tokens[i] + ' is invalid'
    endif

    pos = stregex(tokens[i], size_regex, LENGTH = length)
    if pos eq -1 then begin
      sizes[i] = 1
      names[i] = tokens[i]
    endif else begin
      sizes[i] = fix(strmid(tokens[i], pos + 1, length - 2))
      names[i] = strmid(tokens[i], 0, pos)
    endelse
  endfor

  ; Create variables that will be used to create the template
  ; For now, default the types to always be double precision
  fieldcount = long(total(sizes))
  fieldgroups = lonarr(fieldcount)
  fieldnames = strarr(fieldcount)
  ;fieldtypes = replicate(5, fieldcount)  
  fieldlocations = strsplit(record, ',') - 1 > 0
  ;fieldtypes[0]=7
 
  ; Duplicate group values and name values if there are arrays in the data
  start_index = 0
  end_index = 0
  for i = 0, n_elements(sizes) - 1 do begin
    end_index = start_index + sizes[i] - 1
    fieldgroups[start_index : end_index] = i
    fieldnames[start_index : end_index] = names[i]
    ;;start_index += sizes[i]  
    start_index=start_index+sizes[i]
  endfor
  
  template = create_struct(             $
      'version', 1.0,                   $
      'datastart', 1L,                  $
      'delimiter', 44B,                 $
      'missingvalue', !VALUES.F_NAN,    $
      'commentsymbol', '',              $
      'fieldcount', fieldcount,         $
      'fieldtypes', fieldtypes,         $
      'fieldnames', fieldnames,         $
      'fieldlocations', fieldlocations, $
      'fieldgroups', fieldgroups        $
  )
  return, template
end
