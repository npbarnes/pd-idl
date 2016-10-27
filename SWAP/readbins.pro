pro readbins, bins
    bin = {energy_bin, e_mid: 0.0, e_min: 0.0, e_max: 0.0}
    openr,3,'swap_e_bin.dat'
    readf,3,bin
    bins = [bin]
    while not(eof(3)) do begin
       readf,3,bin
       bins = [bin,bins]
    endwhile

    close,3

    return
end
