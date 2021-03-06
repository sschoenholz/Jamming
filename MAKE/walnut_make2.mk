

INCLUDE = \
-I$(srcDIR) \
-I/opt/pkg/arpack++/icc/include \
-I/opt/pkg/eigen/3.1.3-gcc/include/eigen3 \
-I/usr/include/suitesparse

LIBRARY = \
-L/opt/pkg/arpack/icc \
-L/opt/pkg/netcdf/4.3.0-icc/lib \
-L/opt/pkg/hdf5/1.8.10-patch1-icc/lib

#SuiteSparseLINK = -lamd -lcholmod -lcolamd -lccolamd -lcamd -lumfpack -blas
#SuiteSparseLINK = -lumfpack -lamd -lufconfig -lcholmod -lcolamd -lblas
SuiteSparseLINK = -lumfpack -lamd -lcholmod -lcolamd -lblas
netCDFLINK = -lnetcdf_c++ -lnetcdf
hdf5LINK = -lhdf5_hl -lhdf5 -lz
intelLINK = -lifcore -limf -lm
ArpackLINK = -larpack_LINUX
LINK = $(SuiteSparseLINK)  -lgfortran $(ArpackLINK) $(netCDFLINK) $(hdf5LINK) $(intelLINK) $(SuiteSparseLINK)



