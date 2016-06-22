main_peter.pro: This program loads the response and sets up loops analagous
		 to those on side 4 of swap_instrument_response_b.pptx
		 This program plots on the screen.

main_plot_swap_response_funct_b.pro: The purpose is to load and plot the SWAP 
		instrument response functions and use that to plot simulated
		count rates for a given sweep.This code uses a Maxwellian 
		distribution.  See documentation in 
		swap_instrument_response_b.pptx. This program plots to
		files located in plots/.

Running Programs:
   First update the paths in run_once.csh.

   linux> source run_once.csh       
   idl> .run main_plot_swap_response_funct_b.pro


Data Files:
   20060120_20151231_600.fit        ;Trajectory information at 10 minute resolution
   20060120_20151231_3600.fit       ;Trajectory information at 1 hour resolution
   20060120_20151231_86400.fit      ;Trajectory information at 1 day resolution
   er_table_lookup_3.csv            ;SWAP detector efficiency 
   fin_arr_ebea_ang_eb_bg_corr.sav  ;SWAP energy-theta response
   w_phi.sav                        ;SWAP phi response

Subroutines:
   get_swap_resp.pro
   color_bar.pro
   load_example.pro
   plot_color_4_eb_ea_final.pro
   plot_count_rate_distribution.pro
   plot_eff.pro
   plot_vr_sc_b.pro
   plot_phi_response.pro
   read_traj.pro
   calculate_flux_d.pro
   count_rate_c2.pro
   count_rate_d1.pro
   flux.pro
   convert_date_time_to_year_arr.pro
   convert_date_time_to_year_arr_spice.pro
   makedots.pro
   plot_type_position_cland.pro
   tvimage.pro
   csv_template_type.pro
   fsc_color.pro
   fsc_brewer.tbl
   query_ascii.pro
   read_ascii.pro
   read_csv_fieldtypes.pro
   tvimage_congrid.pro
   tvimage_error.pro

IDL Routines:
   loadct.pro
   filepath.pro
   interpol.pro
   int_tabulated.pro
   julday.pro
   reverse.pro
   strsplit.pro
   uniq.pro
   path_sep.pro

Astron Routines (Reading FITS Files):
   match.pro
   mrd_hread.pro
   mrd_skip.pro
   fxmove.pro
   fxpar.pro
   fxposit.pro
   gettok.pro
   mrdfits.pro
   mrd_struct.pro




