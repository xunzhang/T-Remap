icpc -O3 -g -I//opt/intel/netCDF/include *.cxx -c
icpc -o thu-remap *.o -L/opt/intel/netCDF/lib/ -lnetcdf

