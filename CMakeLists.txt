# iAMAS-Model
cmake_minimum_required(VERSION 3.23)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Functions/iAMAS_Functions.cmake)
get_iamas_version(IAMAS_VERSION)
message(STATUS "iAMAS Version: ${IAMAS_VERSION}")
project(
  iAMAS
  LANGUAGES C CXX Fortran
  VERSION ${IAMAS_VERSION}
  DESCRIPTION "iAMAS - Model for Prediction Across Scales")

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})
include(GNUInstallDirs)

# Options
set(iAMAS_ALL_CORES atmosphere init_atmosphere)
set(iAMAS_CORES
    atmosphere
    CACHE STRING "iAMAS cores to build. Options: ${iAMAS_ALL_CORES}")
if(iAMAS_CORES MATCHES " ") # Convert strings separated with spaces to CMake
  # list separated with ';'
  string(REPLACE " " ";" iAMAS_CORES ${iAMAS_CORES})
  set(iAMAS_CORES
      ${iAMAS_CORES}
      CACHE STRING "iAMAS cores to build. Options: ${iAMAS_ALL_CORES}" FORCE)
endif()
option(DO_PHYSICS "Use built-in physics schemes." TRUE)
option(DO_CHEMISTRY "Use built-in chemistry schemes." TRUE)
option(iAMAS_DOUBLE_PRECISION "Use double precision 64-bit Floating point." TRUE)
option(iAMAS_PROFILE "Enable GPTL profiling" OFF)
option(iAMAS_OPENMP "Enable OpenMP" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(USE_PIO2 "Use PIO 2 library" TRUE)
option(USE_DNN "Use DNN" FALSE)
option(USE_KPP "Use KPP" FALSE)
option(USE_SW "USE_SW" FALSE)
option(USE_SERIALIZATION "USE_SERIALIZATION" TRUE)
option(USE_DEBUG "Use Debug" FALSE)

message(STATUS "[OPTION] iAMAS_CORES: ${iAMAS_CORES}")
message(STATUS "[OPTION] iAMAS_DOUBLE_PRECISION: ${iAMAS_DOUBLE_PRECISION}")
message(STATUS "[OPTION] iAMAS_PROFILE: ${iAMAS_PROFILE}")
message(STATUS "[OPTION] iAMAS_OPENMP: ${iAMAS_OPENMP}")
message(STATUS "[OPTION] BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")

# Build product output locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Set default build type to RelWithDebInfo
if(NOT CMAKE_BUILD_TYPE)
  message(
    STATUS
      "Setting default build type to Release.  Specify CMAKE_BUILD_TYPE to override."
  )
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "CMake Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# Detect iAMAS git version
if(NOT iAMAS_GIT_VERSION)
  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --dirty
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      OUTPUT_VARIABLE _iamas_git_version
      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    set(_iamas_git_version "Unknown")
  endif()
  set(iAMAS_GIT_VERSION
      ${_iamas_git_version}
      CACHE STRING "iAMAS-Model git version")
endif()
# Dependencies
find_package(OpenMP REQUIRED COMPONENTS Fortran CXX)
find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(NetCDF REQUIRED COMPONENTS Fortran C)
find_package(PnetCDF REQUIRED COMPONENTS Fortran)
find_package(PIO REQUIRED COMPONENTS Fortran C)
if(iAMAS_PROFILE)
  find_package(GPTL REQUIRED)
endif()
if(USE_DNN)
  # search for torch. @zhxia
  enable_language(CXX)
  find_package(Torch REQUIRED)
  if(Torch_FOUND)
    message(STATUS "Found Torch: ${Torch_DIR}")
  else()
    message(STATUS "Torch not found! Please install Torch!")
  endif()
endif()
if(USE_KPP)
  # search for Flex & Bison. @zhxia
  find_package(FLEX REQUIRED)
  if(NOT FLEX_FOUND)
    message(FATAL_ERROR "FLEX not found! It is required to build the internal KPP pre-processor.")
  endif()
  message(STATUS "Found FLEX: ${FLEX_EXECUTABLE}")

  find_package(BISON REQUIRED)
  if(NOT BISON_FOUND)
    message(FATAL_ERROR "BISON not found! It is required to build the internal KPP pre-processor.")
  endif()
  message(STATUS "Found BISON: ${BISON_EXECUTABLE}")
endif()

set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

find_package(Kokkos CONFIG)
if(Kokkos_FOUND)
  message(STATUS "Found Kokkos: ${Kokkos_DIR} (version \"${Kokkos_VERSION}\")")
else()
  message(STATUS "Kokkos not found externally. Fetching via FetchContent.")
  include(FetchContent)
  find_package(Patch REQUIRED)
  FetchContent_Declare(
    Kokkos
    URL https://github.com/kokkos/kokkos/releases/download/4.7.00/kokkos-4.7.00.tar.gz
    URL_HASH MD5=24cd603e2a047fc8d67d814f33769f54
    DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/downloaded #
    DOWNLOAD_EXTRACT_TIMESTAMP ON #
    PATCH_COMMAND ${Patch_EXECUTABLE} -p1 -i
                  ${PROJECT_SOURCE_DIR}/cmake/kokkos.patch FIND_PACKAGE_ARGS)
  FetchContent_MakeAvailable(Kokkos)
endif()

find_package(flcl CONFIG)
if(flcl_FOUND)
  message(STATUS "Found flcl: ${flcl_DIR} (version \"${flcl_VERSION}\")")
else()
  message(STATUS "flcl not found externally. Fetching via FetchContent.")
  include(FetchContent)
  set(FLCL_BUILD_EXAMPLES OFF)
  set(FLCL_BUILD_TESTS OFF)
  FetchContent_Declare(
    flcl
    GIT_REPOSITORY https://github.com/waldonchen/kokkos-fortran-interop/
    GIT_TAG develop
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/downloaded/flcl FIND_PACKAGE_ARGS)
  FetchContent_MakeAvailable(flcl)
endif()

# Find C pre-processor
if(CMAKE_C_COMPILER_ID MATCHES GNU)
  find_program(CPP_EXECUTABLE NAMES cpp REQUIRED)
  set(CPP_EXTRA_FLAGS -traditional)
elseif(CMAKE_C_COMPILER_ID MATCHES "(Apple)?Clang")
  find_program(CPP_EXECUTABLE NAMES cpp REQUIRED)
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  find_program(CPP_EXECUTABLE NAMES cpp REQUIRED)
  set(CPP_EXTRA_FLAGS -traditional -D_MPI)
else()
  message(STATUS "Unknown compiler: ${CMAKE_C_COMPILER_ID}")
  set(CPP_EXECUTABLE ${CMAKE_C_COMPILER})
endif()

# Common Variables

# Fortran module output directory for build interface
set(iAMAS_MODULE_DIR
    ${PROJECT_NAME}/module/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION}
)
# Install Fortran module directory
install(DIRECTORY ${CMAKE_BINARY_DIR}/${iAMAS_MODULE_DIR}/
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/${iAMAS_MODULE_DIR}/)

# Location of common subdriver module compiled by each cores
set(iAMAS_MAIN_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/iamas.F)
set(iAMAS_SUBDRIVER_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/iamas_subdriver.F)

# Create targets
if(USE_DNN)
  if(NOT "${iAMAS_CORES}" STREQUAL "atmosphere")
    message(ERROR "DNN is only supported for atmosphere core! Current core: ${iAMAS_CORES}")
  endif()
endif()

if(USE_KPP)
  if(NOT "${iAMAS_CORES}" STREQUAL "atmosphere")
    message(ERROR "KPP is only supported for atmosphere core! Current core: ${iAMAS_CORES}")
  endif()
endif()

add_subdirectory(src/external/ezxml) # Target: iAMAS::external::ezxml
if(ESMF_FOUND)
  message(STATUS "Configure iAMAS for external ESMF")
  add_definitions(-DiAMAS_EXTERNAL_ESMF_LIB -DiAMAS_NO_ESMF_INIT)
  add_library(${PROJECT_NAME}::external::esmf ALIAS esmf)
else()
  message(STATUS "Configure iAMAS for internal ESMF")
  add_subdirectory(src/external/esmf_time_f90) # Target:
  # iAMAS::external::esmf_time
endif()
add_subdirectory(src/tools/input_gen) # Targets: namelist_gen, streams_gen
add_subdirectory(src/tools/registry) # Targets: iamas_parse_<core_name>
add_subdirectory(src/framework) # Target: iAMAS::framework
add_subdirectory(src/operators) # Target: iAMAS::operators

foreach(_core IN LISTS iAMAS_CORES)
  add_subdirectory(src/${_core}) # Target: iAMAS::core::<core_name>
endforeach()

# Package config
include(CMakePackageConfigHelpers)

# Build-tree target exports
export(
  EXPORT ${PROJECT_NAME}ExportsExternal
  NAMESPACE ${PROJECT_NAME}::external::
  FILE ${PROJECT_NAME}-targets-external.cmake)
export(
  EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake)
export(
  EXPORT ${PROJECT_NAME}ExportsCore
  NAMESPACE ${PROJECT_NAME}::core::
  FILE ${PROJECT_NAME}-targets-core.cmake)

# CMake Config file install location
set(CONFIG_INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
# Install iAMAS-supplied Find<Pkg>.cmake modules for use by downstream CMake
# dependencies
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules
        DESTINATION ${CONFIG_INSTALL_DESTINATION})

# <pkgname>-config.cmake: build-tree Variables to export for use from build-tree
set(BINDIR ${CMAKE_BINARY_DIR}/bin)
set(CORE_DATADIR_ROOT ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
set(CMAKE_MODULE_INSTALL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
configure_package_config_file(
  cmake/PackageConfig.cmake.in ${PROJECT_NAME_LOWER}-config.cmake
  INSTALL_DESTINATION .
  INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}
  PATH_VARS BINDIR CORE_DATADIR_ROOT CMAKE_MODULE_INSTALL_PATH)

# <pkgname>-config.cmake: install-tree Variables to export for use from
# install-tree
set(BINDIR ${CMAKE_INSTALL_BINDIR})
set(CORE_DATADIR_ROOT ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})
set(CMAKE_MODULE_INSTALL_PATH ${CONFIG_INSTALL_DESTINATION}/Modules)
configure_package_config_file(
  cmake/PackageConfig.cmake.in install/${PROJECT_NAME_LOWER}-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION}
  PATH_VARS BINDIR CORE_DATADIR_ROOT CMAKE_MODULE_INSTALL_PATH)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/install/${PROJECT_NAME_LOWER}-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# <pkgname>-config-version.cmake
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME_LOWER}-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

# package-targets.cmake and package-targets-<build-type>.cmake
install(
  EXPORT ${PROJECT_NAME}ExportsExternal
  NAMESPACE ${PROJECT_NAME}::external::
  FILE ${PROJECT_NAME_LOWER}-targets-external.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(
  EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME_LOWER}-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(
  EXPORT ${PROJECT_NAME}ExportsCore
  NAMESPACE ${PROJECT_NAME}::core::
  FILE ${PROJECT_NAME_LOWER}-targets-core.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})
