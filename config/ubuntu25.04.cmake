### Ubuntu 25.04, 64-bit - Dell XPS 13 (9305) ###

if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE {LITTLE})
endif()

# Compiler for parallel build
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build
   set(ENV{FC} mpifort)
   set(CMAKE_Fortran_COMPILER mpifort)
   set(USER_Fortran_FLAGS_RELEASE "-cpp -O3 -fconvert=little-endian -ffast-math -mtune=native -march=native")
   #   set(USER_Fortran_FLAGS_RELEASE "-O0 -ggdb")
   add_definitions(-DUSE_MPI -DUSE_MPI_IO)
   set(CMAKE_BUILD_TYPE RELEASE)

# Compiler for serial build
else()
   set(ENV{FC} gfortran)
   set(CMAKE_Fortran_COMPILER gfortran)
   
   if    ( ${BUILD_TYPE} STREQUAL "BIG" )
      set(USER_Fortran_FLAGS_RELEASE "-fconvert=big-endian -ffpe-summary=none -O3 -ffast-math -mtune=native -march=native")
      set(CMAKE_BUILD_TYPE RELEASE)
      
   elseif( ${BUILD_TYPE} STREQUAL "LITTLE" )
      set(USER_Fortran_FLAGS_RELEASE "-cpp -fconvert=little-endian -ffpe-summary=none -O3 -ffast-math -mtune=native -march=native")
      set(CMAKE_BUILD_TYPE RELEASE)
      
   elseif( ${BUILD_TYPE} STREQUAL "PROFILE" )
      set(USER_Fortran_FLAGS_DEBUG "-cpp -fconvert=little-endian -O0 -pg -ggdb -ffpe-summary=none")
      set(CMAKE_BUILD_TYPE DEBUG)
      
   else()
      set(USER_Fortran_FLAGS_DEBUG "-cpp -fconvert=little-endian -O0 -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow")#,underflow,precision,denormal")
      set(USER_Fortran_FLAGS_DEBUG "-cpp -fconvert=little-endian -O0 -ggdb -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow")
      add_definitions(-D_DEBUG)
      set(CMAKE_BUILD_TYPE DEBUG)

   endif()
    
endif()

set(GNU_SED "gsed")

# Libraries - FFTW
add_definitions(-DUSE_FFTW)
#set(FFTW_INCLUDE_DIR   "/usr/local/include")
#set(FFTW_LIB           "/usr/local/lib/libfftw3.a")
set(FFTW_LIB           "-lfftw3")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

# Libraries - NETCDF
add_definitions(-DUSE_NETCDF)
set(NC_INCLUDE_DIR     "/usr/include")
set(NC_LIB             "-L/usr/lib -lnetcdff -lnetcdf")
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${NC_INCLUDE_DIR})
set(LIBS ${LIBS} ${NC_LIB})
