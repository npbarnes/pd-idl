function read_csv_fieldtypes, csvfile, fieldtypes
  template = CSV_TEMPLATE_TYPE(csvfile,fieldtypes)
  return, READ_ASCII(csvfile, TEMPLATE = template)
end
