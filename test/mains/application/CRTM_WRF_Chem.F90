!
! CRTM_WRF_Chem
!
! Check/example program for the CRTM Forward and K-Matrix functions.
!
! NOTE: No results are output or compared in this program.
!       At this stage, it is included as an example only,
!       and to determine that the CRTM library built
!       correctly such that it can be linked to create
!       an executable.
!

! Aerosol-Index mapping:
!  --- CRTM ---
!  1 = Dust
!  2 = Sea salt-SSAM
!  3 = Sea salt-SSCM1
!  4 = Sea salt-SSCM2
!  5 = Sea salt-SSCM3
!  6 = Organic carbon
!  7 = Black carbon
!  8 = Sulfate
!  --- CMAQ ---
!  1 = Dust
!  2 = Soot
!  3 = Water soluble
!  4 = Sulfate
!  5 = Sea salt
!  6 = Water
!  7 = Insoluble
!  8 = dust-like
!  --- GOCART-GEOS5 ---
!  1, 2, 3, 4, 5  = Dust 1, 2, 3, 4, 5
!  6, 7, 8, 9, 10 = Sea salt 1, 2, 3, 4, 5
!  11, 12 = Organic carbon 1, 2
!  13, 14 = Black carbon 1, 2
!  15, 16 = Sulfate 1, 2
!  17, 18, 19 = Nitrate 1, 2, 3
!  --- NAAPS ---
!  1 = Dust
!  2 = Smoke
!  3 = Sea Salt
!  4 = Anthropogenic and Biogenic Fine Particles

PROGRAM CRTM_WRF_Chem

  ! ============================================================================
  ! STEP 1. **** ENVIRONMENT SETUP FOR CRTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  USE CRTM_RTSolution_Define, ONLY: CRTM_RTSolution_WriteFile
  USE netcdf
  USE CRTM_SpcCoeff, ONLY: SC
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================



  ! --------------------------
  ! Some non-CRTM-y Parameters
  ! --------------------------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'CRTM_WRF_Chem'

  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  ! Directory location of coefficients
#ifdef LITTLE_ENDIAN
  CHARACTER(*), PARAMETER :: ENDIAN_TYPE='little_endian'
#else
  CHARACTER(*), PARAMETER :: ENDIAN_TYPE='big_endian'
#endif

  CHARACTER(*), PARAMETER :: COEFFICIENT_PATH='./testinput/'
  CHARACTER(*), PARAMETER :: NC_COEFFICIENT_PATH='./testinput/'

  ! Aerosol/Cloud coefficient format
  CHARACTER(*), PARAMETER :: Coeff_Format = 'Binary'
  !CHARACTER(*), PARAMETER :: Coeff_Format = 'netCDF'

  ! Aerosol/Cloud coefficient scheme
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'CRTM'
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'CMAQ'
  CHARACTER(*), PARAMETER :: Aerosol_Model = 'GOCART-GEOS5'
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'NAAPS'
  CHARACTER(*), PARAMETER :: Cloud_Model   = 'CRTM'

  ! Directory location of results for comparison [NOT USED YET]
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/'

  ! Profile dimensions
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_LAYERS    = 30
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 31  !tbd
  INTEGER, PARAMETER :: N_TIME = 1
  INTEGER, PARAMETER :: N_LAT = 1
  INTEGER, PARAMETER :: N_LON = 1
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'v.modis_aqua'/)

  ! Some pretend geometry angles. The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! WRF-Chem input files
  INTEGER :: FileId
  INTEGER :: NF90_Status
  INTEGER :: VarId
  REAL(Double), DIMENSION(N_LAYERS)   :: rh, Reff_i, Reff_j, Reff_k, &
                                         PresLayer, Temperature
  REAL(Double), DIMENSION(0:N_LAYERS) :: PresLevel
  CHARACTER(*), PARAMETER :: File_Ancillary = '/Users/dangch/Documents/WRF/WRF_output/WRFChem_map2_CRTM.nc4'
  CHARACTER(*), PARAMETER :: File_WRF_Rhodz = '/Users/dangch/Documents/WRF/WRF_output/seoul_only/rhodz.Seoul.nc'
  CHARACTER(*), PARAMETER :: File_WRF_RH    = '/Users/dangch/Documents/WRF/WRF_output/seoul_only/rh.wrfout_d02_2023-04-13_00:00:00.Seoul.nc'
  CHARACTER(*), PARAMETER :: File_WRF_Aero  = '/Users/dangch/Documents/WRF/WRF_output/seoul_only/wrfout_d02_2023-04-13_00:00:00.Seoul.nc'
  REAL(Double), DIMENSION(N_LAYERS)   :: rhodz
  REAL(Double), DIMENSION(N_LON, N_LAT, N_LAYERS, N_TIME) :: antha, asoa1i, asoa1j, asoa2i, asoa2j,      &
                                       asoa3i, asoa3j, asoa4i, asoa4j, bsoa1i, bsoa1j, bsoa2i, bsoa2j,   &
                                       bsoa3i, bsoa3j, bsoa4i, bsoa4j, clai, claj, eci, ecj, naai, naaj, &
                                       nh4ai, nh4aj, no3ai, no3aj, orgpai, orgpaj, p25i, p25j, seas,     &
                                       so4ai, so4aj, soila

  ! ============================================================================

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  CHARACTER(256) :: AerosolCoeff_File
  CHARACTER(256) :: AerosolCoeff_Format
  CHARACTER(256) :: CloudCoeff_File
  CHARACTER(256) :: CloudCoeff_Format
  CHARACTER(256) :: Aerosol_Scheme
  CHARACTER(256) :: Cloud_Scheme
  CHARACTER(256) :: output_nc_file
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels
  INTEGER :: l, m, n, nc
  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: geo(N_PROFILES)


  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type)              :: atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)

  ! 3c. Define the K-MATRIX variables
  ! ---------------------------------
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)
  ! ============================================================================


  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Check/example program for the CRTM Forward and K-Matrix functions using '//&
    ENDIAN_TYPE//' coefficient datafiles', &
    'CRTM Version: '//TRIM(Version) )

  ! ============================================================================
  ! STEP 4. **** INITIALIZE THE CRTM ****
  !
  ! 4a. Initialise all the sensors at once
  ! --------------------------------------
  ! ... Cloud coefficient information
  IF ( Cloud_Model /= 'CRTM' ) THEN
      Cloud_Scheme = Cloud_Model//'.'
  ELSE
      Cloud_Scheme = ' '
  END IF
  ! ... Aerosol coefficient information
  IF ( Aerosol_Model /= 'CRTM' ) THEN
      Aerosol_Scheme = Aerosol_Model//'.'
  ELSE
      Aerosol_Scheme = ' '
  END IF
  ! ... Coefficient table format
  IF ( Coeff_Format == 'Binary' ) THEN
    AerosolCoeff_Format = 'Binary'
    AerosolCoeff_File   = 'AerosolCoeff.'//TRIM(Aerosol_Scheme)//'bin'
    CloudCoeff_Format   = 'Binary'
    CloudCoeff_File     = 'CloudCoeff.'//TRIM(Cloud_Scheme)//'bin'
  ELSE IF ( Coeff_Format == 'netCDF' ) THEN
    AerosolCoeff_Format = 'netCDF'
    AerosolCoeff_File   = 'AerosolCoeff.'//TRIM(Aerosol_Scheme)//'nc4'
    CloudCoeff_Format   = 'netCDF'
    CloudCoeff_File     = 'CloudCoeff.'//TRIM(Cloud_Scheme)//'nc4'
  END IF

  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  err_stat = CRTM_Init( SENSOR_ID, &
                        chinfo, &
                        Aerosol_Model, &
                        AerosolCoeff_Format, &
                        AerosolCoeff_File, &
                        Cloud_Model, &
                        CloudCoeff_Format, &
                        CloudCoeff_File, &
                        File_Path=COEFFICIENT_PATH, &
                        NC_File_Path=NC_COEFFICIENT_PATH, &
                        Quiet=.TRUE.)
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  DO n = 1, N_SENSORS
    WRITE( *,'(7x,i0," from ",a)' ) &
      CRTM_ChannelInfo_n_Channels(chinfo(n)), TRIM(SENSOR_ID(n))
  END DO


  ! ============================================================================



  ! Begin loop over sensors
  ! ----------------------
  Sensor_Loop: DO n = 1, N_SENSORS

    ! ==========================================================================
    ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
    !
    ! 5a. Determine the number of channels
    !     for the current sensor
    ! ------------------------------------
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

    ! 5b. Allocate the ARRAYS
    ! -----------------------
    ALLOCATE( rts( n_channels, N_PROFILES ), &
              atm_K( n_channels, N_PROFILES ), &
              sfc_K( n_channels, N_PROFILES ), &
              rts_K( n_channels, N_PROFILES ), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF


    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF

    ! The output K-MATRIX structure
    CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF

    ! The output structure
    CALL CRTM_RTSolution_Create( rts, N_LAYERS )
    IF ( ANY(.NOT. CRTM_RTSolution_Associated(rts)) ) THEN
      Message = 'Error allocating CRTM RTSolution structures'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
    ! ==========================================================================

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    CALL Load_Atm_Data()
    CALL Load_Sfc_Data()

    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    !  The Sensor_Scan_Angle is optional.
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle = ZENITH_ANGLE, &
                                 Sensor_Scan_Angle   = SCAN_ANGLE )
    ! ==========================================================================




    ! ==========================================================================
    ! STEP 7. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
    !
    ! 7a. Zero the K-matrix OUTPUT structures
    ! ---------------------------------------
    CALL CRTM_Atmosphere_Zero( atm_K )
    CALL CRTM_Surface_Zero( sfc_K )


    ! 7b. Inintialize the K-matrix INPUT so
    !     that the results are dTb/dx
    ! -------------------------------------
    rts_K%Radiance               = ZERO
    rts_K%Brightness_Temperature = ONE
    ! ==========================================================================

    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    WRITE( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(SENSOR_ID(n))

    ! 8a. The forward model
    ! ---------------------
    err_stat = CRTM_Forward( atm        , &  ! Input
                             sfc        , &  ! Input
                             geo        , &  ! Input
                             chinfo(n:n), &  ! Input
                             rts          )  ! Output
    IF ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM Forward Model for '//TRIM(SENSOR_ID(n))
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF


    ! 8b. The K-matrix model
    ! ----------------------
    err_stat = CRTM_K_Matrix( atm        , &  ! FORWARD  Input
                              sfc        , &  ! FORWARD  Input
                              rts_K      , &  ! K-MATRIX Input
                              geo        , &  ! Input
                              chinfo(n:n), &  ! Input
                              atm_K      , &  ! K-MATRIX Output
                              sfc_K      , &  ! K-MATRIX Output
                              rts          )  ! FORWARD  Output
    IF ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM K-Matrix Model for '//TRIM(SENSOR_ID(n))
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF

    ! 8b. The AOD model
    ! ----------------------
    err_stat = CRTM_AOD( atm        , &
                         chinfo(n:n), &
                         rts )
    IF ( err_stat /= SUCCESS ) THEN
      Message = 'Error in CRTM AOD Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
    ! ==========================================================================

   ! ============================================================================
   ! 8c. **** OUTPUT THE RESULTS TO SCREEN ****
   !
   ! User should read the user guide or the source code of the routine
   ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
   ! select the needed variables for outputs.  These variables are contained
   ! in the structure RTSolution.
   ! DO m = 1, N_PROFILES
   !   WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(Sensor_Id(n))
   !   DO l = 1, n_Channels
   !     WRITE( *, '(/5x,"Channel ",i0," results")') chinfo(n)%Sensor_Channel(l)
   !     CALL CRTM_RTSolution_Inspect(rts(l,m))
   !     CALL CRTM_RTSolution_Inspect(rts_K(l,m))
   !     CALL CRTM_Atmosphere_Inspect(atm_K(l,m))
   !     CALL CRTM_Surface_Inspect(sfc_K(l,m))
   !
   !   END DO
   !   CALL CRTM_Atmosphere_Inspect(atm(m))
   !   CALL CRTM_Surface_Inspect(sfc(m))
   ! END DO

     ! 8d. **** OUTPUT THE RESULTS TO NetCDF files ****
     PRINT *, '  channel     wavelength'
     DO l = 1, n_Channels
       WRITE( *,'(I5,8x,f12.5,2ES14.5,3ES17.5)') l,10000.0/SC(1)%wavenumber(l)
     END DO

     output_nc_file = TRIM('WRF_Chem.RTS.'//TRIM(SENSOR_ID(n))//'.nc')
     err_stat = CRTM_RTSolution_WriteFile( output_nc_file, rts, NetCDF=.TRUE., Quiet=.TRUE. )
     PRINT *, 'CRTM RTS are saved in file: ', output_nc_file
     IF ( err_stat /= SUCCESS ) THEN
       message = 'Error calling CRTM_RTSolution_WriteFile for '//TRIM(SENSOR_ID(n))
       CALL Display_Message( PROGRAM_NAME, message, FAILURE )
       STOP
     END IF

    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT SENSOR ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------
    CALL CRTM_Atmosphere_Destroy(atm_K)
    CALL CRTM_Atmosphere_Destroy(atm)


    ! 9b. Deallocate the arrays
    ! -------------------------
    DEALLOCATE(rts, rts_K, sfc_k, atm_k, STAT = alloc_stat)
    ! ==========================================================================

  END DO Sensor_Loop


  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF
  ! ==========================================================================




  ! ==========================================================================
  ! 11. **** CREATE A SIGNAL FILE FOR TESTING SUCCESS ****
  !
  ! This step is just to allow the CRTM library build process
  ! to detect success or failure at the shell level
  CALL SignalFile_Create()
  ! ==========================================================================


CONTAINS


  ! ==========================================================================
  !                Below are some internal procedures that load the
  !                necessary input structures with some pretend data
  ! ==========================================================================
  ! Internal subprogam to load some test profile data
  !
  SUBROUTINE Load_Atm_Data()
    ! Local variables
    INTEGER :: nc
    INTEGER :: k1, k2

    CALL Load_WRF_Data()
    CALL Load_Ancillary_Data()

    ! 4a.1 Profile #1
    ! ---------------
    ! ...Profile and absorber definitions
    atm(1)%Climatology         = US_STANDARD_ATMOSPHERE
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID                 , O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    ! ...Profile data
    ! a. Dummy variables that won't impact WRF AOD but required for CRTM
    ! ... If reading from a file
    atm(1)%Level_Pressure = PresLevel
    atm(1)%Pressure = PresLayer
    atm(1)%Temperature = Temperature
    ! ... If assigned
    ! atm(1)%Level_Pressure = &
    ! (/ 7.14000000e-01, 3.51721379e+01, 6.96302759e+01, 1.04088414e+02, &
    !    1.38546552e+02, 1.73004690e+02, 2.07462828e+02, 2.41920966e+02, &
    !    2.76379103e+02, 3.10837241e+02, 3.45295379e+02, 3.79753517e+02, &
    !    4.14211655e+02, 4.48669793e+02, 4.83127931e+02, 5.17586069e+02, &
    !    5.52044207e+02, 5.86502345e+02, 6.20960483e+02, 6.55418621e+02, &
    !    6.89876759e+02, 7.24334897e+02, 7.58793034e+02, 7.93251172e+02, &
    !    8.27709310e+02, 8.62167448e+02, 8.96625586e+02, 9.31083724e+02, &
    !    9.65541862e+02, 1.00000000e+03/)
    ! atm(1)%Pressure = &
    ! (/ 7.14000000e-01, 3.73568667e+01, 7.39997333e+01, 1.10642600e+02, &
    !    1.47285467e+02, 1.83928333e+02, 2.20571200e+02, 2.57214067e+02, &
    !    2.93856933e+02, 3.30499800e+02, 3.67142667e+02, 4.03785533e+02, &
    !    4.40428400e+02, 4.77071267e+02, 5.13714133e+02, 5.50357000e+02, &
    !    5.86999867e+02, 6.23642733e+02, 6.60285600e+02, 6.96928466e+02, &
    !    7.33571333e+02, 7.70214200e+02, 8.06857067e+02, 8.43499933e+02, &
    !    8.80142800e+02, 9.16785667e+02, 9.53428533e+02, 9.90071400e+02, &
    !    1.02671427e+03, 1.063357130e+3, 1.10000000e+03/)
    ! atm(1)%Temperature = &
    !    (/ 2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, &
    !       2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, &
    !       2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, &
    !       2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, &
    !       2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02, 2.38e+02 /)

    ! Relative_Humidity, rh from WRF, convert to fraction
    atm(1)%Relative_Humidity = rh * 1e-2  !

    ! Gases
    atm(1)%Absorber(:,1) = &
    (/4.187E-03_fp,4.401E-03_fp,4.250E-03_fp,3.688E-03_fp,3.516E-03_fp,3.739E-03_fp,3.694E-03_fp,3.449E-03_fp, &
      3.228E-03_fp,3.212E-03_fp,3.245E-03_fp,3.067E-03_fp,2.886E-03_fp,2.796E-03_fp,2.704E-03_fp,2.617E-03_fp, &
      2.568E-03_fp,2.536E-03_fp,2.506E-03_fp,2.468E-03_fp,2.427E-03_fp,2.438E-03_fp,2.493E-03_fp,2.543E-03_fp, &
      2.586E-03_fp,2.632E-03_fp,2.681E-03_fp,2.703E-03_fp,2.636E-03_fp,2.512E-03_fp/)

    atm(1)%Absorber(:,2) = &
    (/3.035E+00_fp,3.943E+00_fp,4.889E+00_fp,5.812E+00_fp,6.654E+00_fp,7.308E+00_fp,7.660E+00_fp,7.745E+00_fp, &
      7.696E+00_fp,7.573E+00_fp,7.413E+00_fp,7.246E+00_fp,7.097E+00_fp,6.959E+00_fp,6.797E+00_fp,6.593E+00_fp, &
      6.359E+00_fp,6.110E+00_fp,5.860E+00_fp,5.573E+00_fp,5.253E+00_fp,4.937E+00_fp,4.625E+00_fp,4.308E+00_fp, &
      3.986E+00_fp,3.642E+00_fp,3.261E+00_fp,2.874E+00_fp,2.486E+00_fp,2.102E+00_fp/)


    ! !Load CO2 absorber data if there are three absorrbers
    ! IF ( atm(1)%n_Absorbers > 2 ) THEN
    !   atm(1)%Absorber_Id(3)    = CO2_ID
    !   atm(1)%Absorber_Units(3) = VOLUME_MIXING_RATIO_UNITS
    !   atm(1)%Absorber(:,3)     = 380.0_fp
    ! END IF


    ! Load aerosol data
    ! * The aerosol index (idx_aerosol) is ranked based on
    ! antha,
    ! asoa1i, asoa1j, asoa2i, asoa2j, asoa3i, asoa3j, asoa4i, asoa4j,
    ! bsoa1i, bsoa1j, bsoa2i, bsoa2j, bsoa3i, bsoa3j, bsoa4i, bsoa4j, &
    ! clai, claj, eci, ecj, naai, naaj, &
    ! nh4ai, nh4aj, no3ai, no3aj, orgpai, orgpaj, seas,  &
    ! so4ai, so4aj, soila
    ! * Total of 31 aerosols --> max(idx_aerosol) = 31
    ! * SUBROUTINE Map_Aerosol(idx_profile, idx_aerosol, aerotype, radius, conc)

    ! Dust, aerotype = 3
    CALL Map_Aerosol(1, 1,  3, Reff_k, antha )
    CALL Map_Aerosol(1, 31, 3, Reff_k, soila)

    ! Sea salt, i,j,k mode --> aerotype = 6,7,9
    CALL Map_Aerosol(1, 18, 6, Reff_i, clai)
    CALL Map_Aerosol(1, 19, 7, Reff_j, claj)
    CALL Map_Aerosol(1, 22, 6, Reff_i, naai)
    CALL Map_Aerosol(1, 23, 7, Reff_j, naaj)
    CALL Map_Aerosol(1, 28, 9, Reff_k, seas)

    ! Organic carbon hydrophilic, aerotype = 12
    CALL Map_Aerosol(1, 2, 12, Reff_i, asoa1i)
    CALL Map_Aerosol(1, 3, 12, Reff_j, asoa1j)
    CALL Map_Aerosol(1, 4, 12, Reff_i, asoa2i)
    CALL Map_Aerosol(1, 5, 12, Reff_j, asoa2j)
    CALL Map_Aerosol(1, 6, 12, Reff_i, asoa3i)
    CALL Map_Aerosol(1, 7, 12, Reff_j, asoa3j)
    CALL Map_Aerosol(1, 8, 12, Reff_i, asoa4i)
    CALL Map_Aerosol(1, 9, 12, Reff_j, asoa4j)

    CALL Map_Aerosol(1, 10, 12, Reff_i, bsoa1i)
    CALL Map_Aerosol(1, 11, 12, Reff_j, bsoa1j)
    CALL Map_Aerosol(1, 12, 12, Reff_i, bsoa2i)
    CALL Map_Aerosol(1, 13, 12, Reff_j, bsoa2j)
    CALL Map_Aerosol(1, 14, 12, Reff_i, bsoa3i)
    CALL Map_Aerosol(1, 15, 12, Reff_j, bsoa3j)
    CALL Map_Aerosol(1, 16, 12, Reff_i, bsoa4i)
    CALL Map_Aerosol(1, 17, 12, Reff_j, bsoa4j)

    CALL Map_Aerosol(1, 26, 12, Reff_i, orgpai)
    CALL Map_Aerosol(1, 27, 12, Reff_j, orgpaj)

    ! Black carbon hydrophilic, aerotype = 14
    CALL Map_Aerosol(1, 20, 14, Reff_i, eci)
    CALL Map_Aerosol(1, 21, 14, Reff_j, ecj)

    ! Sulfate, aerotype = 15, 16
    CALL Map_Aerosol(1, 29, 15, Reff_i, so4ai)
    CALL Map_Aerosol(1, 30, 16, Reff_j, so4aj)

    ! Nitrate, aerotype = 15, 16
    CALL Map_Aerosol(1, 24, 15, Reff_i, nh4ai)
    CALL Map_Aerosol(1, 25, 16, Reff_j, nh4aj)

  END SUBROUTINE Load_Atm_Data

  !
  ! Internal subprogam to load some test surface data
  !
  SUBROUTINE Load_Sfc_Data()


    ! 4a.0 Surface type definitions for default SfcOptics definitions
    !      For IR and VIS, this is the NPOESS reflectivities.
    ! ---------------------------------------------------------------
    INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
    INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
    INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
    INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
    INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
    INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics



    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc(1)%Land_Coverage     = 0.1_fp
    sfc(1)%Land_Type         = TUNDRA_SURFACE_TYPE
    sfc(1)%Land_Temperature  = 272.0_fp
    sfc(1)%Lai               = 0.17_fp
    sfc(1)%Soil_Type         = COARSE_SOIL_TYPE
    sfc(1)%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc(1)%Water_Coverage    = 0.5_fp
    sfc(1)%Water_Type        = SEA_WATER_TYPE
    sfc(1)%Water_Temperature = 275.0_fp
    ! ...Snow coverage characteristics
    sfc(1)%Snow_Coverage    = 0.25_fp
    sfc(1)%Snow_Type        = FRESH_SNOW_TYPE
    sfc(1)%Snow_Temperature = 265.0_fp
    ! ...Ice surface characteristics
    sfc(1)%Ice_Coverage    = 0.15_fp
    sfc(1)%Ice_Type        = FRESH_ICE_TYPE
    sfc(1)%Ice_Temperature = 269.0_fp

  END SUBROUTINE Load_Sfc_Data

  !
  ! Internal subprogam to load WRF-Chem data
  !
  SUBROUTINE Load_WRF_Data()
    ! Read WRF-Chem data

    ! 1. Read Rhodz from file File_WRF_Rhodz
    PRINT *, 'Reading WRF-Chem Data: ', File_WRF_Rhodz
    IF ( .NOT. File_Exists(File_WRF_Rhodz) ) THEN
       WRITE(*,*) 'File does not exist (File_WRF_Rhodz)'
    END IF
    CALL check_nc(NF90_OPEN(File_WRF_Rhodz, NF90_NOWRITE, FileId))
    CALL Read_NetCDF_Var_1D(FileId, 'rhodz', rhodz)
    CALL check_nc(NF90_CLOSE(FileId))

    ! 2. Read RH from file File_WRF_RH
    PRINT *, 'Reading WRF-Chem Data: ', File_WRF_RH
    IF ( .NOT. File_Exists(File_WRF_RH) ) THEN
       WRITE(*,*) 'File does not exist (File_WRF_RH)'
    END IF
    CALL check_nc(NF90_OPEN(File_WRF_RH, NF90_NOWRITE, FileId))
    CALL Read_NetCDF_Var_1D(FileId, 'rh', rh)
    CALL check_nc(NF90_CLOSE(FileId))

    ! 3. Read aerosol concentrations from File_WRF_Aero
    ! All 33 aerosols but p25i, p25j are included, PARAMETER N_AEROSOLS = 31
    ! antha, asoa1i, asoa1j, asoa2i, asoa2j,      &
    ! asoa3i, asoa3j, asoa4i, asoa4j, bsoa1i, bsoa1j, bsoa2i, bsoa2j,   &
    ! bsoa3i, bsoa3j, bsoa4i, bsoa4j, clai, claj, eci, ecj, naai, naaj, &
    ! nh4ai, nh4aj, no3ai, no3aj, orgpai, orgpaj, seas,     &
    ! so4ai, so4aj, soila
    PRINT *, 'Reading WRF-Chem Data: ', File_WRF_Aero
    IF ( .NOT. File_Exists(File_WRF_Aero)) THEN
       WRITE(*,*) 'File does not exist (File_WRF_Aero)'
    END IF
    CALL check_nc(NF90_OPEN(File_WRF_Aero, NF90_NOWRITE, FileId))
    CALL Read_NetCDF_Var_4D(FileId, 'antha' , antha )
    CALL Read_NetCDF_Var_4D(FileId, 'asoa1i', asoa1i)
    CALL Read_NetCDF_Var_4D(FileId, 'asoa1j', asoa1j)
    CALL Read_NetCDF_Var_4D(FileId, 'asoa2i', asoa2i)
    CALL Read_NetCDF_Var_4D(FileId, 'asoa2j', asoa2j)
    CALL Read_NetCDF_Var_4D(FileId, 'asoa3i', asoa3i)
    CALL Read_NetCDF_Var_4D(FileId, 'asoa3j', asoa3j)
    CALL Read_NetCDF_Var_4D(FileId, 'asoa4i', asoa4i)
    CALL Read_NetCDF_Var_4D(FileId, 'asoa4j', asoa4j)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa1i', bsoa1i)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa1j', bsoa1j)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa2i', bsoa2i)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa2j', bsoa2j)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa3i', bsoa3i)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa3j', bsoa3j)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa4i', bsoa4i)
    CALL Read_NetCDF_Var_4D(FileId, 'bsoa4j', bsoa4j)
    CALL Read_NetCDF_Var_4D(FileId, 'clai'  , clai  )
    CALL Read_NetCDF_Var_4D(FileId, 'claj'  , claj  )
    CALL Read_NetCDF_Var_4D(FileId, 'eci'   , eci   )
    CALL Read_NetCDF_Var_4D(FileId, 'ecj'   , ecj   )
    CALL Read_NetCDF_Var_4D(FileId, 'naai'  , naai  )
    CALL Read_NetCDF_Var_4D(FileId, 'naaj'  , naaj  )
    CALL Read_NetCDF_Var_4D(FileId, 'nh4ai' , nh4ai )
    CALL Read_NetCDF_Var_4D(FileId, 'nh4aj' , nh4aj )
    CALL Read_NetCDF_Var_4D(FileId, 'no3ai' , no3ai )
    CALL Read_NetCDF_Var_4D(FileId, 'no3aj' , no3aj )
    CALL Read_NetCDF_Var_4D(FileId, 'orgpai', orgpai)
    CALL Read_NetCDF_Var_4D(FileId, 'orgpaj', orgpaj)
    ! CALL Read_NetCDF_Var_4D(FileId, 'p25i'  , p25i  )
    ! CALL Read_NetCDF_Var_4D(FileId, 'p25j'  , p25j  )
    CALL Read_NetCDF_Var_4D(FileId, 'seas'  , seas  )
    CALL Read_NetCDF_Var_4D(FileId, 'so4ai' , so4ai )
    CALL Read_NetCDF_Var_4D(FileId, 'so4aj' , so4aj )
    CALL Read_NetCDF_Var_4D(FileId, 'soila' , soila )
    CALL check_nc(NF90_CLOSE(FileId))

    ! Convert aerosol concentration from ugkg to kgm-2
    ! C [kg/m^2] = r [ug/kg-dryair] * air_density[kg/m3] * dz [m] * 1e-9
    CALL ugkg_to_kgm2(rhodz, antha)
    CALL ugkg_to_kgm2(rhodz, asoa1i)
    CALL ugkg_to_kgm2(rhodz, asoa1j)
    CALL ugkg_to_kgm2(rhodz, asoa2i)
    CALL ugkg_to_kgm2(rhodz, asoa2j)
    CALL ugkg_to_kgm2(rhodz, asoa3i)
    CALL ugkg_to_kgm2(rhodz, asoa3j)
    CALL ugkg_to_kgm2(rhodz, asoa4i)
    CALL ugkg_to_kgm2(rhodz, asoa4j)
    CALL ugkg_to_kgm2(rhodz, bsoa1i)
    CALL ugkg_to_kgm2(rhodz, bsoa1j)
    CALL ugkg_to_kgm2(rhodz, bsoa2i)
    CALL ugkg_to_kgm2(rhodz, bsoa2j)
    CALL ugkg_to_kgm2(rhodz, bsoa3i)
    CALL ugkg_to_kgm2(rhodz, bsoa3j)
    CALL ugkg_to_kgm2(rhodz, bsoa4i)
    CALL ugkg_to_kgm2(rhodz, bsoa4j)
    CALL ugkg_to_kgm2(rhodz, clai)
    CALL ugkg_to_kgm2(rhodz, claj)
    CALL ugkg_to_kgm2(rhodz, eci)
    CALL ugkg_to_kgm2(rhodz, ecj)
    CALL ugkg_to_kgm2(rhodz, naai)
    CALL ugkg_to_kgm2(rhodz, naaj)
    CALL ugkg_to_kgm2(rhodz, nh4ai)
    CALL ugkg_to_kgm2(rhodz, nh4aj)
    CALL ugkg_to_kgm2(rhodz, no3ai)
    CALL ugkg_to_kgm2(rhodz, no3aj)
    CALL ugkg_to_kgm2(rhodz, orgpai)
    CALL ugkg_to_kgm2(rhodz, orgpaj)
    ! CALL ugkg_to_kgm2(rhodz, p25i)
    ! CALL ugkg_to_kgm2(rhodz, p25j)
    CALL ugkg_to_kgm2(rhodz, seas)
    CALL ugkg_to_kgm2(rhodz, so4ai)
    CALL ugkg_to_kgm2(rhodz, so4aj)
    CALL ugkg_to_kgm2(rhodz, soila)
  END SUBROUTINE Load_WRF_Data

  !
  ! Subroutine for processing WRF-Chem Aerosol Data
  !
  SUBROUTINE ugkg_to_kgm2(rhodz, C)
    REAL(DOUBLE), INTENT(IN    ) :: rhodz(:)
    REAL(DOUBLE), INTENT(IN OUT) :: C(:,:,:,:)
    INTEGER :: nlay, nlat, nlon, ntime
    DO ntime = 1, N_TIME
      DO nlat = 1, N_LAT
      DO nlon = 1, N_LON
        DO nlay = 1, N_LAYERS
          C(nlon,nlat,nlay,ntime) = C(nlon,nlat,nlay,ntime) * rhodz(nlay) * 1.0e-9
        END DO
      END DO
      END DO
    END DO
  END SUBROUTINE ugkg_to_kgm2

  !
  ! Subroutine for loading Ancillary_Data
  !
  SUBROUTINE Load_Ancillary_Data()
    PRINT *, 'Reading Ancillary Data: ', File_Ancillary
    IF ( .NOT. File_Exists(File_Ancillary) ) THEN
       WRITE(*,*) 'File does not exist (File_Ancillary)'
    END IF
    CALL check_nc(NF90_OPEN(File_Ancillary, NF90_NOWRITE, FileId))
    CALL Read_NetCDF_Var_1D(FileId, 'PresLayer', PresLayer)
    CALL Read_NetCDF_Var_1D(FileId, 'PresLevel', PresLevel)
    CALL Read_NetCDF_Var_1D(FileId, 'Temperature', Temperature)
    CALL Read_NetCDF_Var_1D(FileId, 'Reff_i', Reff_i)
    CALL Read_NetCDF_Var_1D(FileId, 'Reff_j', Reff_j)
    CALL Read_NetCDF_Var_1D(FileId, 'Reff_j', Reff_k)
    CALL check_nc(NF90_CLOSE(FileId))
  END SUBROUTINE Load_Ancillary_Data

  !
  ! Subroutine for mapping aerosols
  !
  SUBROUTINE Map_Aerosol(idx_profile, idx_aerosol, aerotype, radius, conc)
    INTEGER, INTENT(IN)      :: idx_profile, idx_aerosol, aerotype
    REAL(DOUBLE), INTENT(IN) :: radius(:), conc(:,:,:,:)
    atm(idx_profile)%Aerosol(idx_aerosol)%Type             = aerotype
    atm(idx_profile)%Aerosol(idx_aerosol)%Effective_Radius = radius
    atm(idx_profile)%Aerosol(idx_aerosol)%Concentration    = conc(1,1,:,1)
  END SUBROUTINE Map_Aerosol

  !
  ! Internal subprogam to create a signal file
  !
  SUBROUTINE SignalFile_Create()
    CHARACTER(256) :: Filename
    INTEGER :: fid
    Filename = '.signal'
    fid = Get_Lun()
    OPEN( fid, FILE = Filename )
    WRITE( fid,* ) TRIM(Filename)
    CLOSE( fid )
  END SUBROUTINE SignalFile_Create

  !
  ! Internal subprogram for NetCDF I/O
  !
  SUBROUTINE check_nc(NF90_Status)
    INTEGER, INTENT(in) :: NF90_Status
    IF(NF90_Status /= NF90_NOERR) THEN
      PRINT *, 'netCDF error!:'
      PRINT *, TRIM(NF90_STRERROR(NF90_Status))
      STOP "Stopped"
    END IF
  END SUBROUTINE check_nc

  SUBROUTINE Read_NetCDF_Var_1D(FileId, VarName, VarArray)
    INTEGER,      INTENT(IN)  :: FileId
    CHARACTER(*), INTENT(IN)  :: VarName
    REAL(Double), INTENT(OUT) :: VarArray(:)
    ! Local variable
    INTEGER :: VarId
    NF90_Status = NF90_INQ_VARID( FileId, VarName, VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID for vairable: ', VarName
    NF90_Status = NF90_GET_VAR( FileId, VarId, VarArray)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR for vairable: ', VarName
  END SUBROUTINE Read_NetCDF_Var_1D

  SUBROUTINE Read_NetCDF_Var_4D(FileId, VarName, VarArray)
    INTEGER,      INTENT(IN)  :: FileId
    CHARACTER(*), INTENT(IN)  :: VarName
    REAL(Double), INTENT(OUT) :: VarArray(:,:,:,:)
    ! Local variable
    INTEGER :: VarId
    NF90_Status = NF90_INQ_VARID( FileId, VarName, VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID for vairable: ', VarName
    NF90_Status = NF90_GET_VAR( FileId, VarId, VarArray)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR for vairable: ', VarName
  END SUBROUTINE Read_NetCDF_Var_4D


END PROGRAM CRTM_WRF_Chem
