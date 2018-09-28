[x] implemented lr4_lowpass and lr4_highpass
[x] factored out helper functions into helpers.h
[x] add an addition parameter MMAPFNAME to control inputs which is used in
    the filename of mmap files in /dev/shm (to help identify different instances
    of the plugins wehn used in PaXoverRack)
[x] rename .so name to t5_3band_parameq_with_shelves
[x] get proper ladspa ids
[x] rename label to 3band_parameq_with_shelves
[x] rename struct and functions from BiquadEq to ThreeBandParaqmEqWithSelhves
