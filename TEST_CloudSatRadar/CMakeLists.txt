
cmake_minimum_required(VERSION 3.16)
project(crtm_app LANGUAGES Fortran)

# CRTM root
set( crtm_root ${CMAKE_CURRENT_SOURCE_DIR}/../)

# Use CRTM stuff
list( APPEND CMAKE_MODULE_PATH  ${crtm_root}/cmake )
include( crtm_compiler_flags )
set( CMAKE_BUILD_TYPE Debug )
#add_compile_options( -warn all,nounused)

# Path for the compiled excutable
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR} )


# Find CRTM lib and mod 
MACRO(subDirList result curdir)
  FILE(GLOB_RECURSE children LIST_DIRECTORIES true RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist ${curdir})
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()
subDirList( subDirs_lib ${crtm_root}/build/src)
subDirList( subDirs_mod ${crtm_root}/build/module)
#set( crtm_lib_file ${crtm_root}/build/src/libcrtm.a)
#set( crtm_mod_root ${crtm_root}/build/module/crtm/GNU/9.2.0)
find_library( crtm_libFile NAMES "libcrtm.a" PATHS ${subDirs_lib} REQUIRED)
find_path( crtm_modPath NAMES "crtm_module.mod" PATHS ${subDirs_mod} REQUIRED)


# OpenMP
find_package( OpenMP COMPONENTS Fortran )

# Find NetCDF package
if(DEFINED ENV{NETCDF_PATH})
  list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF_PATH})
elseif(DEFINED ENV{NETCDF})
  list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF})
elseif(DEFINED ENV{SSEC_NETCDF_DIR})
  list(APPEND CMAKE_PREFIX_PATH $ENV{SSEC_NETCDF_DIR})
endif()
find_package(NetCDF REQUIRED Fortran)


# Compile the project
set( exeFile  run_crtm_Z.exe )
add_executable( ${exeFile}
                run_crtm_Z.f90  )
target_link_libraries(${exeFile} PRIVATE OpenMP::OpenMP_Fortran)
target_link_libraries(${exeFile} PRIVATE NetCDF::NetCDF_Fortran)
#target_link_libraries(${exeFile} PRIVATE crtm )
target_link_libraries(${exeFile} PRIVATE ${crtm_libFile})
target_include_directories(${exeFile} PRIVATE  ${crtm_modPath})

