function read_4sight_hdf5, path

  ; Open the hdf5 file
  file_id = H5F_OPEN(path)

  ; Open the appropriate dataset within the file 
  data_id = H5D_OPEN(file_id, 'measurement0/genraw/data')
  
  img = H5D_READ(data_id)
  
  return, img
  

end