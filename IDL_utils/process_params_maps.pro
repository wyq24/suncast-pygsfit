; Function to read HDF5 files and create arrays for data, error, and fit maps
FUNCTION read_and_create_arrays, file_list, num_maps

  ;temp structure:
  
  tmap = read_hdf5(file_list[0], /shallow)
  temp_struct = tmap.maps.datamaps.(0)
  ; Initialize arrays
  new_datamaps_array = make_array(num_maps, n_elements(file_list),value=temp_struct)
  new_errmaps_array = make_array(num_maps, n_elements(file_list), value=temp_struct)
  new_fitmaps_array = make_array(num_maps, n_elements(file_list), value=temp_struct)

  ; Loop through files and read data
  FOR f = 0, n_elements(file_list) - 1 DO BEGIN
    cmaps = read_hdf5(file_list[f], /shallow)

    FOR i = 0, num_maps - 1 DO BEGIN
      new_datamaps_array[i, f] = cmaps.maps.datamaps.(i)
      new_errmaps_array[i, f] = cmaps.maps.errmaps.(i)
      new_fitmaps_array[i, f] = cmaps.maps.fitmaps.(i)
    ENDFOR
  ENDFOR

  RETURN, {datamaps: new_datamaps_array, errmaps: new_errmaps_array, fitmaps: new_fitmaps_array}
END

; Main program
PRO process_files, file_list

  num_maps = 50
  num_files = n_elements(file_list)
  template_map = read_hdf5(file_list[0], /shallow)
  ; Read data from files and create arrays
  result = read_and_create_arrays(file_list, num_maps)

  ; Initialize the new_maps structure with dynamic arrays
  new_maps = {B: make_array(num_files, value=template_map.maps.B), $
    CHISQR: make_array(num_files, value=template_map.maps.CHISQR), $
    DATAMAPS: result.datamaps, $
    DELTA: make_array(num_files, value=template_map.maps.DELTA), $
    ERRB: make_array(num_files, value=template_map.maps.ERRB), $
    ERRDELTA: make_array(num_files, value=template_map.maps.ERRDELTA), $
    ERRE_MAX: make_array(num_files, value=template_map.maps.ERRE_MAX), $
    ERRMAPS: result.errmaps, $
    ERRN_NTH: make_array(num_files, value=template_map.maps.ERRN_NTH), $
    ERRN_TH: make_array(num_files, value=template_map.maps.ERRN_TH), $
    ERRPEAKFLUX: make_array(num_files, value=template_map.maps.ERRPEAKFLUX), $
    ERRPEAKFREQ: make_array(num_files, value=template_map.maps.ERRPEAKFREQ), $
    ERRTHETA: make_array(num_files, value=template_map.maps.ERRTHETA), $
    ERRT_E: make_array(num_files, value=template_map.maps.ERRT_E), $
    ERRWB: make_array(num_files, value=template_map.maps.ERRWB), $
    ERRWNTH: make_array(num_files, value=template_map.maps.ERRWNTH), $
    E_MAX: make_array(num_files, value=template_map.maps.E_MAX), $
    FITMAPS: result.fitmaps, $
    N_NTH: make_array(num_files, value=template_map.maps.N_NTH), $
    N_TH: make_array(num_files, value=template_map.maps.N_TH), $
    PEAKFLUX: make_array(num_files, value=template_map.maps.PEAKFLUX), $
    PEAKFREQ: make_array(num_files, value=template_map.maps.PEAKFREQ), $
    RESIDUAL: make_array(num_files, value=template_map.maps.RESIDUAL), $
    THETA: make_array(num_files, value=template_map.maps.THETA), $
    T_E: make_array(num_files, value=template_map.maps.T_E), $
    WB: make_array(num_files, value=template_map.maps.WB), $
    WNTH: make_array(num_files, value=template_map.maps.WNTH)}

  ; Fill the structure fields
  FOR f = 0, num_files - 1 DO BEGIN
    cmaps = read_hdf5(file_list[f], /shallow)
    new_maps.B[f] = cmaps.maps.B
    new_maps.CHISQR[f] = cmaps.maps.CHISQR
    new_maps.DELTA[f] = cmaps.maps.DELTA
    new_maps.ERRB[f] = cmaps.maps.ERRB
    new_maps.ERRDELTA[f] = cmaps.maps.ERRDELTA
    new_maps.ERRE_MAX[f] = cmaps.maps.ERRE_MAX
    new_maps.ERRN_NTH[f] = cmaps.maps.ERRN_NTH
    new_maps.ERRN_TH[f] = cmaps.maps.ERRN_TH
    new_maps.ERRPEAKFLUX[f] = cmaps.maps.ERRPEAKFLUX
    new_maps.ERRPEAKFREQ[f] = cmaps.maps.ERRPEAKFREQ
    new_maps.ERRTHETA[f] = cmaps.maps.ERRTHETA
    new_maps.ERRT_E[f] = cmaps.maps.ERRT_E
    new_maps.ERRWB[f] = cmaps.maps.ERRWB
    new_maps.ERRWNTH[f] = cmaps.maps.ERRWNTH
    new_maps.E_MAX[f] = cmaps.maps.E_MAX
    new_maps.N_NTH[f] = cmaps.maps.N_NTH
    new_maps.N_TH[f] = cmaps.maps.N_TH
    new_maps.PEAKFLUX[f] = cmaps.maps.PEAKFLUX
    new_maps.PEAKFREQ[f] = cmaps.maps.PEAKFREQ
    new_maps.RESIDUAL[f] = cmaps.maps.RESIDUAL
    new_maps.THETA[f] = cmaps.maps.THETA
    new_maps.T_E[f] = cmaps.maps.T_E
    new_maps.WB[f] = cmaps.maps.WB
    new_maps.WNTH[f] = cmaps.maps.WNTH
  ENDFOR

  ; Save the structure
  save, filename='/Users/walterwei/Desktop/test_all.sav', new_maps

END

; Example usage
file_list = ['/Users/walterwei/Downloads/20220511/gsfit_test/pix_t/slf_final_XX_t21_allbd_params_map.h5', $
'/Users/walterwei/Downloads/20220511/gsfit_test/pix_t/slf_final_XX_t23_allbd_params_map.h5', $
'/Users/walterwei/Downloads/20220511/gsfit_test/pix_t/slf_final_XX_t24_allbd_params_map.h5',$
'/Users/walterwei/Downloads/20220511/gsfit_test/pix_t/slf_final_XX_t25_allbd_params_map.h5' ]
process_files, file_list
end