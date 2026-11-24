#!/bin/bash

# default
OPTION=0
PLATFORM="general"
iAMAS_CORE=atmosphere
NJOBS=4
BUILD_PATH="build"
BUILD_TYPE="Release"
COMPILER=intel
DNN="FALSE"
KPP="FALSE"
USE_SW="FALSE"
USE_SERIALIZATION="TRUE"
# functions
parse_arguments() {
    while [[ "$#" -gt 0 ]]; do
      case $1 in
        -O|--option) OPTION="$2"; shift ;;
        -P|--platform) PLATFORM="$2"; shift ;;
        -C|--iamas_core) iAMAS_CORE="$2"; shift ;;
        -J|--njobs) NJOBS="$2"; shift ;;
        -E|--compiler) COMPILER="$2"; shift ;;
        -B|--build_path) BUILD_PATH="$2"; shift ;;
        -T|--build_type) BUILD_TYPE="$2"; shift ;;
        --dnn) DNN="$2"; shift ;;
        --kpp) KPP="$2"; shift ;;
        *) echo "unknown argument: $1"; exit 1 ;;
      esac
      shift
    done
}
set_compiler() {
    local compiler_id=$1
    # setting compiler
    if [ "$compiler_id" == "intel" ]; then
      if [ "$COMPILER" == "openmpi" ]; then
         export C=mpicc
         export CC=mpicc
         export CXX=mpicxx
         export FC=mpif90
         export FFLAGS=" -cpp" #zhxia add "-fpp" for macro definition in fortran code
         export CFLAGS=""
         export FFLAGS_DEBUG=" -cpp -g -fcheck=all -fbacktrace -Wall -Wextra"
         export CFLAGS_DEBUG="-g -O0"
         export CXXFLAGS_DEBUG="-g -O0"
      elif [ "$COMPILER" == "intel" ]; then
         export C=mpiicc
         export CC=mpiicc
         export CXX=mpicxx
         export FC=mpiifort
         export FFLAGS=" -fpp -fp-model precise" #zhxia add "-fp-model precise" for avoiding radical optimization of floating-point numbers
         export CFLAGS=""
         export CXXFLAGS=""
      	 export FFLAGS_DEBUG=" -fpp -CU -CB -check all -fpe0 -traceback"
         export CFLAGS_DEBUG=" -g -Og"
         export CXXFLAGS_DEBUG=" -g -Og"
      fi
    elif [ "$compiler_id" == "sw64" ]; then
      export C=gcc
      export CC=mpicc
      export CXX=swumpicxx
      export FC=mpif90

      export FFLAGS=" -cpp -mieee -mftz -msimd -mfma -fno-common -finit-local-zero -zerouv -fdefault-real-8"
      export CFLAGS=" -mieee -mftz -msimd -mfma -fno-common -zerouv"
      export CXXFLAGS=" -mieee -mftz -msimd -mfma -fno-common -zerouv"
      export FFLAGS_DEBUG=" -cpp -ffree-line-length-none -fconvert='big-endian' -ffree-form -fboudns-check -fbacktrace -ffpe-trap=invalid,zero,overflow"
      export CFLAGS_DEBUG=""
      export CXXFLAGS_DEBUG=""
    else
      echo "Unknown compiler '$compiler_id'"
      exit 1
    fi
    export CMAKE_C_COMPILER=$CC
    export CMAKE_CXX_COMPILER=$CXX
    export CMAKE_Fortran_COMPILER=$FC

    if [ "$BUILD_TYPE" == "Debug" ]; then
        export CMAKE_C_FLAGS=$CFLAGS_DEBUG
        export CMAKE_CXX_FLAGS=$CFLAGS_DEBUG
        export CMAKE_Fortran_FLAGS=$FFLAGS_DEBUG
        DEBUG_FLAG="TRUE"
    else
        export CMAKE_C_FLAGS=$CFLAGS
        export CMAKE_CXX_FLAGS=$CXXFLAGS
        export CMAKE_Fortran_FLAGS=$FFLAGS
    fi
}

set_platform() {
    local platform_id=$1
    if [ "$platform_id" == "general" ]; then
        set_compiler "intel"
        USE_SERIALIZATION="TRUE"
        export BUILD_SHARED_LIBS="OFF"
        # setting external library
        if [ $DNN == "TRUE" ]; then
            export Torch_DIR=/work/home/zhxia/env/libtorch_cpu
        fi
        if [ $KPP == "TRUE" ]; then
            export KPP_FLEX_LIB_DIR=/work/home/sflin/anaconda3/pkgs/flex-2.6.4-h58526e2_1004/lib
            export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${KPP_FLEX_LIB_DIR}:"
        fi
    elif [ "$platform_id" == "core01" ]; then
        set_compiler "intel"
        USE_SERIALIZATION="TRUE"
        export BUILD_SHARED_LIBS="OFF"
        # setting external library
        if [ $DNN == "TRUE" ]; then
            export Torch_DIR=/home/zhxia/env/libtorch_cpu
        fi
        if [ $KPP == "TRUE" ]; then
            export KPP_FLEX_LIB_DIR=/home/zhxia/anaconda3/envs/keyan/lib
            export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${KPP_FLEX_LIB_DIR}:"
        fi
    elif [ "$platform_id" == "h800" ]; then
        set_compiler "intel"
        USE_SERIALIZATION="TRUE"
        export BUILD_SHARED_LIBS="OFF"
        # setting external library
        if [ $DNN == "TRUE" ]; then
            export Torch_DIR=/test1/lacar_zhxia/env/software/mathlib_for_intel_2022/libtorch
        fi
    elif [ "$platform_id" == "sunway@swustc" ]; then
        set_compiler "sw64"
        USE_SW="TRUE"
        USE_SERIALIZATION="FALSE"
        export BUILD_SHARED_LIBS="OFF"
        export PATH=/home/export/online1/mdt00/shisuan/swustc/qcjiang/hefei/tddft/dependencies/swhmpi:$PATH
        module purge
        module load gcc/8.3.1
        module load swgcc/1449
        module load mpi/swuc
        export NETCDF=/home/export/online1/mdt00/shisuan/swustc/ftao/local-saca
        export PnetCDF_ROOT=/home/export/online1/mdt00/shisuan/swustc/ftao/local-saca/pnetcdf
        export PIO=/home/export/online1/mdt00/shisuan/swustc/ftao/local-saca/pio-2.new
        export HDF5=/home/export/online1/mdt00/shisuan/swustc/ftao/local-saca
        export Kokkos_DIR=/home/export/online1/mdt00/shisuan/swustc/xinqu/env/kokkos
        export flcl_DIR=/home/export/online1/mdt00/shisuan/swustc/xinqu/env/flcl
        export ZLIB=/home/export/online1/mdt00/shisuan/swustc/ftao/local-saca
        export iAMAS_TOOL_DIR=/home/export/online1/mdt00/shisuan/swustc/jungu/dongze/iAMAS-tools-gnu
  elif [ "$platform_id" == "sunway@swustciamas" ]; then
        set_compiler "sw64"
        USE_SW="TRUE"
        USE_SERIALIZATION="FALSE"
        export BUILD_SHARED_LIBS="OFF"
        source $HOME/.local/sw/setenv-release-SEA-1449
        source $HOME/.local/sw/setenv-mpi-swuc

        export NETCDF=/home/export/online1/mdt00/shisuan/swustciamas/data/software/local-saca
        export Kokkos_DIR=/home/export/online1/mdt00/shisuan/swustciamas/data/software/kokkos
        export flcl_DIR=/home/export/online1/mdt00/shisuan/swustciamas/data/software/flcl
        export PnetCDF_ROOT=/home/export/online1/mdt00/shisuan/swustciamas/data/software/local-saca/pnetcdf
        export PIO=/home/export/online1/mdt00/shisuan/swustciamas/data/software/local-saca/pio-2.new
        export HDF5=/home/export/online1/mdt00/shisuan/swustciamas/data/software/local-saca
        export ZLIB=/home/export/online1/mdt00/shisuan/swustciamas/data/software/local-saca
        export iAMAS_TOOL_DIR=/home/export/online1/mdt00/shisuan/swustciamas/data/software/iAMAS-tools-gnu
        if [ $KPP == "TRUE" ]; then
            export KPP_FLEX_LIB_DIR=/home/export/online1/mdt00/shisuan/swustciamas/data/software/flex/lib
            export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${KPP_FLEX_LIB_DIR}:"
        fi
  elif [ "$platform_id" == "sunway@swiamas" ]; then
        set_compiler "sw64"
        USE_SW="TRUE"
        USE_SERIALIZATION="FALSE"
        export BUILD_SHARED_LIBS="OFF"
        export NETCDF=/home/export/online2/ifir/swiamas/data/software/local-saca
        export Kokkos_DIR=/home/export/online2/ifir/swiamas/data/software/kokkos
        export flcl_DIR=/home/export/online2/ifir/swiamas/data/software/flcl
        export PnetCDF_ROOT=/home/export/online2/ifir/swiamas/data/software/local-saca/pnetcdf
        export PIO=/home/export/online2/ifir/swiamas/data/software/local-saca/pio-2.new
        export HDF5=/home/export/online2/ifir/swiamas/data/software/local-saca
        export ZLIB=/home/export/online2/ifir/swiamas/data/software/local-saca
        export iAMAS_TOOL_DIR=/home/export/online2/ifir/swiamas/data/software/iAMAS-tools-gnu
        if [ $KPP == "TRUE" ]; then
            export KPP_FLEX_LIB_DIR=/home/export/online2/ifir/swiamas/data/software/flex/lib
            export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${KPP_FLEX_LIB_DIR}:"
        fi

    else
        echo "Unknown platform '$platform_id'"
        exit 1
    fi
    echo "- Compiler Setting"
    echo "  - CMAKE_C_COMPILER="$CMAKE_C_COMPILER
    echo "  - CMAKE_CXX_COMPILER="$CMAKE_CXX_COMPILER
    echo "  - CMAKE_Fortran_COMPILER="$CMAKE_Fortran_COMPILER
    echo "  - CMAKE_C_FLAGS="$CMAKE_C_FLAGS
    echo "  - CMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS
    echo "  - CMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS
    echo "  - BUILD_SHARED_LIBS="$BUILD_SHARED_LIBS
    if [ $DNN == "TRUE" ]; then
        echo "  - TORCH_FOR_DNN="$Torch_DIR
    fi
    if [ $KPP == "TRUE" ]; then
        echo "  - FLEX_FOR_KPP="$KPP_FLEX_LIB_DIR
    fi
}

buildOrMakeOrClean() {
    local option_id=$1
    local PROJ_PATH=$(pwd)
    if [ "$option_id" == 0 ]; then
        echo "     **Build and make.**    "
        echo "- Log file: "bash.log""
		mkdir -p $BUILD_PATH
        # run CMake
        cmake -S $PROJ_PATH -B $BUILD_PATH \
              -DiAMAS_CORES=$iAMAS_CORE \
              -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
              -DCMAKE_C_COMPILER=$CMAKE_C_COMPILER \
              -DCMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER \
              -DCMAKE_Fortran_COMPILER=$CMAKE_Fortran_COMPILER \
              -DCMAKE_C_FLAGS="$CMAKE_C_FLAGS" \
              -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
              -DCMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS" \
              -DBUILD_SHARED_LIBS=$BUILD_SHARED_LIBS \
              -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
              -DPLATFORM_FLAG=$PLATFORM \
              -DUSE_SW=$USE_SW		\
              -DUSE_SERIALIZATION=$USE_SERIALIZATION		\
              -DUSE_DNN=$DNN \
              -DUSE_KPP=$KPP \
              -DUSE_DEBUG=$DEBUG_FLAG 
        # run make
        cd "$BUILD_PATH"
        make VERBOSE=1 -j"$NJOBS" 
        echo "*************************************"
    elif [ "$option_id" == 1 ]; then
        echo "     **Only Build.**    "
        echo "- Log file: "bash.log""
		mkdir -p $BUILD_PATH
        # run CMake
        cmake -S $PROJ_PATH -B $BUILD_PATH \
              -DiAMAS_CORES=$iAMAS_CORE \
              -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
              -DCMAKE_C_COMPILER=$CMAKE_C_COMPILER \
              -DCMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER \
              -DCMAKE_Fortran_COMPILER=$CMAKE_Fortran_COMPILER \
              -DCMAKE_C_FLAGS="$CMAKE_C_FLAGS" \
              -DCMAKE_CXX_FLAGS="$CMAKE_CXX_FLAGS" \
              -DCMAKE_Fortran_FLAGS="$CMAKE_Fortran_FLAGS" \
              -DBUILD_SHARED_LIBS=$BUILD_SHARED_LIBS \
              -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
              -DPLATFORM_FLAG=$PLATFORM \
              -DUSE_SW=$USE_SW		\
              -DUSE_SERIALIZATION=$USE_SERIALIZATION		\
              -DUSE_DNN=$DNN \
              -DUSE_KPP=$KPP \
              -DUSE_DEBUG=$DEBUG_FLAG 
        echo "*************************************"
    elif [ "$option_id" == 2 ]; then
        echo "     **Only Make.**    "
        echo "- Log file: "bash.log""
        cd "$BUILD_PATH"
        # run make
        make VERBOSE=1 -j"$NJOBS" 
        echo "*************************************"
    elif [ "$option_id" == 3 ]; then
        echo "     **Cleaning project...**    "
        if [ -d "$BUILD_PATH" ]; then
            echo "- Removing build directory: $BUILD_PATH"
            rm -rf "$BUILD_PATH"
        else
            echo "- Build directory not found. Nothing to clean for iAMAS."
        fi

        local KPP_DIR="src/atmosphere/chemistry/schemes/KPP/kpp-3.3.0"
        if [ -f "$KPP_DIR/src/Makefile" ]; then
            echo "- Cleaning KPP source directory: $KPP_DIR"
            (cd "$KPP_DIR/src" && make distclean)
        else
            echo "- KPP Makefile not found. Nothing to clean for KPP."
        fi

        if [ -f "$KPP_DIR/bin/kpp" ]; then
           (cd "$KPP_DIR/bin" && rm kpp)
        fi

        local KPP_MECHANISMS_DIR="src/atmosphere/chemistry/schemes/KPP/mechanisms"
        if [ -d "$KPP_MECHANISMS_DIR" ]; then
            echo "- Cleaning KPP-generated mechanism files from: $KPP_MECHANISMS_DIR"
            for mechanism_path in "$KPP_MECHANISMS_DIR"/*; do
                if [ -d "$mechanism_path" ]; then
                    local mechanism_name=$(basename "$mechanism_path")
                    echo "  - Removing generated files for mechanism: $mechanism_name"
                    rm -f "${mechanism_path}/${mechanism_name}"_*.f90 "${mechanism_path}/${mechanism_name}".log "${mechanism_path}/Makefile_${mechanism_name}"
                fi
            done
        else
            echo "- KPP mechanisms directory not found. Nothing to clean there."
        fi

        if [ -f "bash.log" ]; then
           rm bash.log
        fi
 
        echo "*************************************"
    else
      echo "invalid option: '$option_id'"
      exit 1
    fi
}

#----------------------------------
#          MAIN
#----------------------------------
# arguments of scripts
parse_arguments "$@"

if [ "$OPTION" != 3 ]; then
   # build log
   echo "************Bulding log**************"
   echo "- OPTION="$OPTION"; 0/1/2:Build&Make/Build/Make;"
   echo "- Platform = "$PLATFORM
   echo "- iAMAS_CORE = "$iAMAS_CORE
   echo "- Build path = "$BUILD_PATH
   echo "- Build type = "$BUILD_TYPE
   echo "- Number of Process for MAKE = "$NJOBS

   # setting compiler
   set_platform "$PLATFORM"
   echo "*************************************"
fi

buildOrMakeOrClean $OPTION

