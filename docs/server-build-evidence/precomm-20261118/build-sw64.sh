#!/bin/bash
set -euo pipefail
export LC_ALL=en_US.utf8
export PATH=/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr/bin:/usr/sw/swgcc/swgcc710-tools-SEA-1449/usr_sw/bin:/usr/sw/mpi/mpi_20220608_SEA/bin:$PATH
export iAMAS_TOOL_DIR=/tmp/iamas-host-tools
src=/home/export/online1/mdt00/shisuan/swlicom1/hk/iamas/precomm-20261118/source-sw64
build=/home/export/online1/mdt00/shisuan/swlicom1/hk/iamas/precomm-20261118/build-sw64
wrapper=/home/export/online1/mdt00/shisuan/swlicom1/softwares_sw/v3_swgomp/bin
swgfortran --version
mpif90 -show
cmake -S "$src" -B "$build"  -DCMAKE_C_COMPILER="$wrapper/mpixcc.py"  -DCMAKE_CXX_COMPILER="$wrapper/mpixcxx.py"  -DCMAKE_Fortran_COMPILER="$wrapper/mpixfort.py"  -DCMAKE_BUILD_TYPE=RelWithDebInfo  -DCMAKE_C_FLAGS="-mieee -mftz -msimd -fopenmp"  -DCMAKE_CXX_FLAGS="-mieee -mftz -msimd -fopenmp -fpermissive"  -DCMAKE_Fortran_FLAGS="-mieee -mftz -msimd -fopenmp -finit-local-zero -I$build/iAMAS/module/GNU/7.1.0 -I$build/iAMAS/module/GNU/7.1.0/core_atmosphere"  -DCMAKE_C_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG" -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG" -DCMAKE_Fortran_FLAGS_RELWITHDEBINFO="-O2 -g -DNDEBUG" -DBUILD_SHARED_LIBS=OFF -DiAMAS_CORES=atmosphere -DiAMAS_OPENMP=ON -DUSE_SW=ON -DUSE_SERIALIZATION=OFF  -DPLATFORM_FLAG=sunway@swgomp  -DNetCDF_ROOT=/usr/sw/yyzlib/local-saca  -DNetCDF_C_CONFIG_EXECUTABLE=/tmp/iamas-netcdf-config/nc-config  -DNetCDF_Fortran_CONFIG_EXECUTABLE=/tmp/iamas-netcdf-config/nf-config  -DPnetCDF_ROOT=/usr/sw/yyzlib/local-saca/pnetcdf-1-12-1  -DPnetCDF_CONFIG_EXE=/tmp/iamas-netcdf-config/pnetcdf-config  -DPIO_ROOT=/tmp/iamas-gcc7-deps/pio  -DKokkos_DIR=/tmp/iamas-gcc7-deps/kokkos/lib64/cmake/Kokkos  -Dflcl_DIR=/tmp/iamas-gcc7-deps/flcl/lib64/cmake/flcl  -DOpenMP_C_FLAGS=-fopenmp -DOpenMP_CXX_FLAGS=-fopenmp -DOpenMP_Fortran_FLAGS=-fopenmp  -DOpenMP_C_LIB_NAMES= -DOpenMP_CXX_LIB_NAMES= -DOpenMP_Fortran_LIB_NAMES=
cmake --build "$build" --target atmosphere -j8
