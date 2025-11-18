if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE PARALLEL)
endif()

set(USER_Fortran_FLAGS  " -fpp ${USER_profile_FLAGS} -nbs -save-temps -heap-arrays -simd -vec-threshold50 -unroll-aggressive ${USER_omp_FLAGS} " )

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" ) # compiler for parallel build
  set(ENV{FC} mpifort)
  set(CMAKE_Fortran_COMPILER mpifort)
  set(USER_Fortran_FLAGS_RELEASE  " -march=core-avx2 -mtune=core-avx2 -qopt-prefetch -O3 -ipo" )
  #set(USER_Fortran_FLAGS_RELEASE  " -march=skylake-avx512 -axcommon-avx512,SSE4.2 -qopt-prefetch -O3 -ipo" )
  #set(USER_Fortran_FLAGS_RELEASE  " -axCORE-AVX2 -qopt-prefetch -O3 -ipo" )
  add_definitions(-DUSE_MPI)
  set(CMAKE_BUILD_TYPE RELEASE)

elseif ( ${BUILD_TYPE} STREQUAL "GPU" ) # compiler for GPU build
  set(ENV{FC} nvfortran)
  set(CMAKE_Fortran_COMPILER nvfortran)
  set(USER_Fortran_FLAGS          " -cpp -g -O2 -Mnofma " ) #DEBUG
  #set(USER_Fortran_FLAGS          " -Mpreprocess -Minfo=accel -O0 -acc -Kieee -g " ) #DEBUG
  #set(USER_Fortran_FLAGS          " -Mpreprocess -Minfo=accel -O3 -acc -Kieee " ) #BASIC - Note: -Mcpp does not work... and dont use -mp (OpenMP)
  #set(USER_Fortran_FLAGS          " -Mpreprocess -Mnostack_arrays -Mvect=simd -tp=haswell -Mprefetch -O3 -Kieee " ) #-Munroll=?? -Mipa=fast,inline " ) #TEST-OPT

  add_definitions(-DNO_ASSUMED_RANKS)
  ##add_definitions(-D_DEBUG)
  set(CMAKE_BUILD_TYPE RELEASE)

else() # compiler for serial build
  set(ENV{FC} ifort)
  set(CMAKE_Fortran_COMPILER ifort)

  if    ( ${BUILD_TYPE} STREQUAL "SERIAL" )
   set(USER_Fortran_FLAGS_RELEASE  " -march=core-avx2 -mtune=core-avx2 -qopt-prefetch -O3 -ipo" )
   #set(USER_Fortran_FLAGS_RELEASE  " -march=skylake-avx512 -axcommon-avx512,SSE4.2 -qopt-prefetch -O3 -ipo" )
   set(CMAKE_BUILD_TYPE RELEASE)

  else()
    set(USER_Fortran_FLAGS_DEBUG    " -g -traceback -debug all ")
    add_definitions(-D_DEBUG)
    set(CMAKE_BUILD_TYPE DEBUG)

  endif()

endif()

add_definitions(-DUSE_FFTW)
set(FFTW_INCLUDE_DIR   "/sw/spack-levante/fftw-3.3.10-fnfhvr/include")
set(FFTW_LIB           "/sw/spack-levante/fftw-3.3.10-fnfhvr/lib/libfftw3.a")
set(INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(LIBS ${FFTW_LIB})

add_definitions(-DUSE_NETCDF)
if ( ${BUILD_TYPE} STREQUAL "GPU" ) # compiler for GPU build
  set(NC_INCLUDE_DIR     "/sw/spack-levante/netcdf-fortran-4.6.1-m5p7yn/include")
  set(NC_LIB             "-L/sw/spack-levante/netcdf-fortran-4.6.1-m5p7yn/lib -lnetcdff -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.6.1-m5p7yn/lib")
else()
  set(NC_INCLUDE_DIR     "/sw/spack-levante/netcdf-fortran-4.5.3-pvmcx6/include")
  set(NC_LIB             "-L/sw/spack-levante/netcdf-fortran-4.5.3-pvmcx6/lib -lnetcdff -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.3-pvmcx6/lib")
endif()
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${NC_INCLUDE_DIR})
set(LIBS ${LIBS} ${NC_LIB})
