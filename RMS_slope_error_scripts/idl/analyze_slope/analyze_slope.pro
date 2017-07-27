PRO analyze_slope

  pix_4D_x = 999  
  pix_4D_y = 1004
  
  fov_4D_x = 87.23  ;mm
  fov_4D_y = 87.75 ; mm
  
  dx = fov_4D_x / pix_4D_x  ; mm / pix
  dy = fov_4D_y / pix_4D_y  ; mm / pix
  
  print, dx, dy
  
  int_len = 8.0 ; mm
  filt_len = int_len / 2  ; mm
  
  lam = 632.8e-6 ; wavelength in mm
  
  

  h5_path = '/media/byrdie/KINGSTON/7.27.17_M2SurfaceTesting/meas_0.h5'
  
  surf = read_4sight_hdf5(h5_path)
  
  ;atv, surf
  
  m_val = surf[0,0]  ; Grab the value of masked pixels

  surf[WHERE((surf EQ m_val) OR (surf EQ 0.0))] = -0.0/0.0
  
  print, wfe2slope(surf, int_len, filt_len, dx=dx, lambda=lam) / 1e-6
  
  

END