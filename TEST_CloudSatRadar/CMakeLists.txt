
cmake_minimum_required(VERSION 3.16)
project(crtm_app LANGUAGES Fortran)

# CMake modules
set( crtm_src_root /home/yma/data/CRTM/Repos/CRTMv3)
list( APPEND CMAKE_MODULE_PATH  ${crtm_src_root}/cmake )

#--
set(CMAKE_BUILD_TYPE Debug )
include( crtm_compiler_flags )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR} )
add_compile_options( -warn all,nounused)

# CRTM lib
set( crtm_lib_file /home/yma/data/CRTM/crtm_lib/CherryTree_noOpenMP/lib64/libcrtm.a)
set( crtm_mod_root /home/yma/data/CRTM/crtm_lib/CherryTree_noOpenMP/module/crtm/Intel/2021.5.0.20211109/)

# OpenMP
find_package( OpenMP COMPONENTS Fortran )

# Find NetCDF package
list(APPEND CMAKE_PREFIX_PATH $ENV{SSEC_NETCDF_DIR})
find_package(NetCDF REQUIRED Fortran)


# Compile the project
set( exeFile  run_crtm_Z.exe )
add_executable( ${exeFile}
                run_crtm_Z.f90  )
target_link_libraries(${exeFile} PUBLIC OpenMP::OpenMP_Fortran)
target_link_libraries(${exeFile} PUBLIC NetCDF::NetCDF_Fortran)
target_link_libraries(${exeFile} PUBLIC ${crtm_lib_file})
target_include_directories(${exeFile} PUBLIC  ${crtm_mod_root})

