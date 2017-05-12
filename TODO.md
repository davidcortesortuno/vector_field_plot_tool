# TODO

- Currently, we interpolate the whole 3D data from an array, VTK
  file, etc. into a 2D mesh grid
  (not good for sampling the 3D field)
  We can slightly solve this by passing a filter (the option is
  implemented) and interpolate the data specifying a range of
  z values

- I have now implemented full 3D interpolation but `scipy.griddata`
  only accepts `linear` methods to interpolate the data. Using 2D
  data allows a `cubic` interpolation which seems more accurate. We
  should think what is best or give the option to interpolate
  using both options.
- 
