#==============================================================================
#
# Makefile for users to simulate radiance 
#
#==============================================================================

#-----------------------------------------------------------------------------
#                          -- Define macros --
#-----------------------------------------------------------------------------

include ./make.macros

# -------------
# This makefile
# -------------

MAKE_FILE = makefile


# ---------------
# Executable file
# ---------------

#EXE_FILE = run_crtm_Z
EXE_FILE = CRTM_radar_test_FWTLADK
#EXE_FILE = CRTM_Forward_CloudSAT

OBJ_FILES = read_cloudsat_data.o $(EXE_FILE).o
# ------------
# Object files
#  v3.0
##INCLUDES = -I/data/smcd8/qliu/Tropics/CRTM_Sim/REL-3.0.1/libsrc -I/data/starfs1/libs/netcdf-4.2-ifort/include 
##LIBRARIES = -qopenmp -L/data/smcd8/qliu/Tropics/CRTM_Sim/REL-3.0.1/libsrc -lcrtm -L/data/starfs1/libs/netcdf-4.2-ifort/lib -lnetcdff

#INCLUDES = -I../src/Build/libsrc -I/data/starfs1/libs/netcdf-4.2-ifort/include 
#LIBRARIES = -qopenmp -L../src/Build/libsrc -lcrtm -L/data/starfs1/libs/netcdf-4.2-ifort/lib -lnetcdff

INCLUDES = -I../build/module/crtm/Intel/19.0.5.20190815 -I/data/starfs1/libs/netcdf-4.2-ifort/include 
LIBRARIES = -qopenmp -L../build/src -lcrtm -L/data/starfs1/libs/netcdf-4.2-ifort/lib -lnetcdff

#INCLUDES = -I../../../CRTM_09202023/libsrc -I/data/starfs1/libs/netcdf-4.2-ifort/include 
#LIBRARIES = -qopenmp -L../../../CRTM_09202023/libsrc -lcrtm -L/data/starfs1/libs/netcdf-4.2-ifort/lib -lnetcdff



all:
	@echo "OS type detected: "`uname -s`
	@case `uname -s` in \
	  "SunOS")   make -f $(MAKE_FILE) test_program $(SUNOS_FLAGS) ;; \
	  "AIX")     make -f $(MAKE_FILE) test_program $(AIX_FLAGS) ;; \
	  "IRIX64" ) make -f $(MAKE_FILE) test_program $(IRIX64_FLAGS) ;; \
	  "Linux" )  make -f $(MAKE_FILE) test_program $(LINUX_FLAGS) ;; \
	  *) echo "This system is not supported" ;; \
       esac


# -- Targets for specific Linux compilers.
#
# *** NOTE: The PGI compiler must be v6 or later but even ***
# ***       that may not work due to compiler bugs        ***
# IBM AIX Compiler
ibm_debug:
	make -f $(MAKE_FILE) test_program $(AIX_FLAGS_DEBUG)

ibm:
	make -f $(MAKE_FILE) test_program $(AIX_FLAGS_PROD)

intel_debug:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_INTEL_DEBUG)

intel:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_INTEL_PROD)

lahey_debug:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_LAHEY_DEBUG)

lahey:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_LAHEY_PROD)

pgi_debug:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_PGI_DEBUG)

pgi:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_PGI_PROD)

g95_debug:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_G95_DEBUG)

g95:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_G95_PROD)

gfortran_debug:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_GFORTRAN_DEBUG)
gfortran:
	make -f $(MAKE_FILE) test_program $(LINUX_FLAGS_GFORTRAN_PROD)


# ----------------
# Make the program
# ----------------

test_program: $(OBJ_FILES)
	$(FL) $(OBJ_FILES) $(EXTRA_FL_FLAGS) $(FL_FLAGS) $(EXE_FILE)



# --------
# Clean up
# --------

clean:
	$(REMOVE) $(OBJ_FILES) $(EXE_FILE) *.mod *.MOD *.stb


# ---------------
# Dependency list
# ---------------
$read_cloudsat_data.o : read_cloudsat_data.f90
$(EXE_FILE).o : $(EXE_FILE).f90
#-----------------------------------------------------------------------------
#                          -- Define default rules --
#-----------------------------------------------------------------------------

include ./make.rules
