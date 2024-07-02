!
! CRTM_ActiveSensor_Module
!
! This module calculates attenuated reflectivity for active sensors within the CRTM. 
! The module is based on the CRTM_AOD_Module. This module designed for radar reflectivity,
! lidar reflectivity and scatterometer for surface. The first application is for radar 
! attenuated reflectivity in dbZ here.
!
!
! CREATION HISTORY:  September 22, 2023
!       Contributed by:  Quanhua Liu
!                        Benjamin Johnson
!                        Isaac Moradi
!                        Yingtao Ma
!
! Inputs are the same as the CRTM model: sensor ID, geometryInfo (e.g. zenith and azimuth angles),
!   Atmosphere and Surface data types.
!  The output will be 
!     RTSolution(channel_idx, profile_idx)%Reflectivity_Attenuated(1:n_Layers).
!     Missing value or zero value is set to -9999.0 for the dbZ values.
!
MODULE CRTM_ActiveSensor_Module


  ! ------------
  ! Module usage
  ! ------------
  USE Type_Kinds,                 ONLY: fp, LLong
  USE ODPS_CoordinateMapping,     ONLY: Geopotential_Height
  USE Message_Handler,            ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,            ONLY: SET,NOT_SET,ZERO,ONE, RT_VMOM, &
                                        MAX_N_LAYERS        , &
                                        MAX_N_PHASE_ELEMENTS, &
                                        MAX_N_LEGENDRE_TERMS, &
                                        MAX_N_STOKES        , &
                                        MAX_N_ANGLES        , &
                                        MAX_N_AZIMUTH_FOURIER, &
                                        MAX_SOURCE_ZENITH_ANGLE, &
                                        MAX_N_STREAMS, &
                                        AIRCRAFT_PRESSURE_THRESHOLD, &
                                        MIN_COVERAGE_THRESHOLD, &
                                        SCATTERING_ALBEDO_THRESHOLD
  USE CRTM_Atmosphere_Define,   ONLY: CRTM_Atmosphere_type, &
                                      CRTM_Atmosphere_IsValid
  USE CRTM_Surface_Define,        ONLY: CRTM_Surface_type, &
                                        CRTM_Surface_IsValid
  USE CRTM_ChannelInfo_Define,  ONLY: CRTM_ChannelInfo_type, &
                                      CRTM_ChannelInfo_n_Channels
  USE CRTM_Options_Define,      ONLY: CRTM_Options_type
  USE CRTM_AtmOptics_Define,    ONLY: CRTM_AtmOptics_type      , &
                                      CRTM_AtmOptics_Associated, &
                                      CRTM_AtmOptics_Create    , &
                                      CRTM_AtmOptics_Destroy   , &
                                      CRTM_AtmOptics_Zero
  USE CRTM_GeometryInfo_Define,   ONLY: CRTM_GeometryInfo_type, &
                                        CRTM_GeometryInfo_SetValue, &
                                        CRTM_GeometryInfo_GetValue
  USE CRTM_GeometryInfo,          ONLY: CRTM_GeometryInfo_Compute
  USE CRTM_Geometry_Define,       ONLY: CRTM_Geometry_type, &
                                        CRTM_Geometry_IsValid
  USE CRTM_Predictor_Define,      ONLY: CRTM_Predictor_type      , &
                                        CRTM_Predictor_Associated, &
                                        CRTM_Predictor_Destroy   , &
                                        CRTM_Predictor_Create
  USE CRTM_Predictor,             ONLY: CRTM_PVar_type => iVar_type, &
       CRTM_Compute_Predictors,CRTM_Compute_Predictors_TL,CRTM_Compute_Predictors_AD
  USE CRTM_AtmAbsorption,         ONLY: CRTM_AAvar_type => iVar_type, &
           CRTM_Compute_AtmAbsorption, CRTM_Compute_AtmAbsorption_TL,CRTM_Compute_AtmAbsorption_AD
  USE CRTM_AerosolScatter,      ONLY: CRTM_Compute_AerosolScatter, &
                                      CRTM_Compute_AerosolScatter_TL, &
                                      CRTM_Compute_AerosolScatter_AD
  USE CRTM_SfcOptics_Define,      ONLY: CRTM_SfcOptics_type      , &
                                        CRTM_SfcOptics_Associated, &
                                        CRTM_SfcOptics_Create    , &
                                        CRTM_SfcOptics_Destroy
  USE CRTM_SfcOptics,             ONLY: CRTM_Compute_SurfaceT,CRTM_Compute_SurfaceT_TL, &
                                        CRTM_Compute_SurfaceT_AD
  USE CRTM_RTSolution_Define,     ONLY: CRTM_RTSolution_type   , &
                                        CRTM_RTSolution_Destroy, &
                                        CRTM_RTSolution_Zero,    &
                                        CRTM_RTSolution_Inspect
  USE CRTM_AncillaryInput_Define, ONLY: CRTM_AncillaryInput_type
  USE CRTM_AerosolCoeff,        ONLY: CRTM_AerosolCoeff_IsLoaded
  USE CRTM_CloudCoeff,            ONLY: CRTM_CloudCoeff_IsLoaded
  USE CRTM_TauCoeff,              ONLY: TC
  USE CRTM_SpcCoeff,              ONLY: SC, &
                                        SpcCoeff_IsVisibleSensor, &
                                        SpcCoeff_IsMicrowaveSensor, &
                                        SpcCoeff_IsInfraredSensor, &
                                        SpcCoeff_IsUltravioletSensor
  USE CRTM_CloudScatter, ONLY: CRTM_Compute_CloudScatter, CRTM_Compute_CloudScatter_TL, &
        CRTM_Compute_CloudScatter_AD

  ! ...CloudScatter
  USE CSvar_Define, ONLY: CSvar_type, &
                          CSvar_Associated, &
                          CSvar_Destroy   , &
                          CSvar_Create  
  ! Internal variable definition modules
  ! ...AerosolScatter
  USE ASvar_Define, ONLY: ASvar_type, &
                          ASvar_Associated, &
                          ASvar_Destroy   , &
                          ASvar_Create

  ! ...AtmOptics
  USE AOvar_Define, ONLY: AOvar_type, &
                          AOvar_Associated, &
                          AOvar_Destroy   , &
                          AOvar_Create
  ! ...Radiative transfer
  USE RTV_Define,   ONLY: RTV_type, &
                          RTV_Associated, &
                          RTV_Destroy   , &
                          RTV_Create
  ! Active Sensors
  USE ActiveSensor_Model, ONLY: Radar_Solution, Radar_Solution_TL, Radar_Solution_AD

  
  ! -----------------------
  ! Disable implicit typing
  ! -----------------------
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Public procedures
  PUBLIC :: CRTM_ActiveSensor
  PUBLIC :: CRTM_ActiveSensor_TL
  PUBLIC :: CRTM_ActiveSensor_AD
  PUBLIC :: CRTM_ActiveSensor_K
  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


CONTAINS


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_ActiveSensor
!
! PURPOSE:
!       Function that calculates layer total optical depth profile at nadir.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_ActiveSensor( Atmosphere       , &
!                                         ChannelInfo      , &
!                                         RTSolution       , &
!                                         Options = Options  )
!
! INPUTS:
!       Atmosphere:     Structure containing the Atmosphere data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Rank-1 (n_Profiles)
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelInfo:    Structure returned from the CRTM_Init() function
!                       that contains the satellite/sensor channel index
!                       information.
!                       UNITS:      N/A
!                       TYPE:       CRTM_ChannelInfo_type
!                       DIMENSION:  Rank-1 (n_Sensors)
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       RTSolution:     Structure containing the layer aerosol optical
!                       profile for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Rank-2 (n_Channels x n_Profiles)
!                       ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       Options:        Options structure containing the optional arguments
!                       for the CRTM.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Options_type
!                       DIMENSION:  Same as input Atmosphere structure
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       - Many of the components of the Options optional input structure
!         are not used in this function. Consult the CRTM User Guide for
!         which Options components are usable for AOD calculations.
!
!:sdoc-:
!--------------------------------------------------------------------------------
  FUNCTION CRTM_ActiveSensor( &
    Atmosphere , &  ! Input, M
    Surface    , &  ! Input, M
    Geometry   , &  ! Input, M
    ChannelInfo, &  ! Input, n_Sensors
    RTSolution , &  ! Output, L x M
    Options    ) &  ! Optional input, M
  RESULT( Error_Status )

    ! Arguments
    TYPE(CRTM_Surface_type),           INTENT(IN)     :: Surface(:)        ! M
    TYPE(CRTM_Geometry_type),          INTENT(IN)     :: Geometry(:)       ! M
    TYPE(CRTM_SfcOptics_type)                         :: SfcOptics
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT)     :: Atmosphere(:)     ! M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)    ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)   ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)        ! M
    TYPE(CRTM_AncillaryInput_type) :: AncillaryInput
    TYPE(CRTM_GeometryInfo_type) :: GeometryInfo
    TYPE(CRTM_Predictor_type)    :: Predictor
    TYPE(CSvar_type)      :: CSvar  ! CloudScatter
    TYPE(AOvar_type)      :: AOvar  ! AtmOptics
    TYPE(RTV_type)        :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_ActiveSensor'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: ln, H2O_idx
    REAL(fp), Allocatable :: Height(:), deltaZ(:)
    REAL(fp) :: transmittance1, transmittance2
    ! Component variables
    TYPE(CRTM_AtmOptics_type) :: AtmOptics
    TYPE(CRTM_PVar_type)  :: PVar   ! Predictor
    TYPE(ASVar_type) :: ASvar
    TYPE(CRTM_AAvar_type) :: AAvar  ! AtmAbsorption
!
    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS

    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN


    ! Check the number of channels
    IF ( SIZE(RTSolution,DIM=1) < n_Channels ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure array too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    ! ...Number of atmospheric profiles.
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(RTSolution,DIM=2) /= n_Profiles ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for RTSolution argument.'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF


    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles

    CALL CRTM_SfcOptics_Create( SfcOptics  , MAX_N_ANGLES, MAX_N_STOKES )

    Allocate( Height(0:Atmosphere(m)%n_Layers), deltaZ(Atmosphere(m)%n_Layers) )
      H2O_idx = 1  ! default
    CALL Geopotential_Height(Atmosphere(m)%Level_Pressure      , & ! Input
                            Atmosphere(m)%Temperature         , & ! Input
                            Atmosphere(m)%Absorber(:, H2O_idx), & ! Input
                            ZERO                    , & ! Input - surface height
                            Height                ) ! Output in km
!    write(6,'(6f12.5)') Height
    deltaZ(:) = Height(0:Atmosphere(m)%n_Layers-1)-Height(1:Atmosphere(m)%n_Layers)
    Atmosphere(m)%Height(:) = Height(0:Atmosphere(m)%n_Layers)
    CALL CRTM_RTSolution_Zero(RTSolution(:,m))
      ! ...Compute derived geometry
    CALL CRTM_GeometryInfo_SetValue( GeometryInfo, Geometry=Geometry(m) )
    CALL CRTM_GeometryInfo_Compute( GeometryInfo )
    CALL CRTM_Compute_SurfaceT( Surface(m), SfcOptics )
      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF
!
      CALL RTV_Create( RTV, MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, Atmosphere(m)%n_Layers )
      
      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
        ! Check whether to skip this profile
        IF ( Options(m)%Skip_Profile ) CYCLE Profile_Loop
      END IF
      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF
      ! Check the RTSolution layer dimension
      IF ( ANY(RTSolution(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! Allocate AtmOptics based on Atmosphere dimension
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS, &
                                  MAX_N_PHASE_ELEMENTS  )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structure for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set the scattering switch
      AtmOptics%Include_Scattering = Options(m)%Include_Scattering
      ! ...Allocate the atmospheric optics internal structure
      CALL AOvar_Create( AOvar, Atmosphere(m)%n_Layers )

      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms = 16

      ! Allocate the scattering internal variables if necessary
      ! ...Cloud
      IF ( Atmosphere(m)%n_Clouds > 0 ) THEN
        CALL CSvar_Create( CSvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Clouds    )
      END IF
      ! ...Aerosol
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Aerosols  )
      END IF
!
      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for channel(l)/sensor(n) count
      ln = 0
      Sensor_Loop: DO n = 1, n_Sensors

        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index

          CALL CRTM_Predictor_Create( &
                   Predictor, &
                   Atmosphere(m)%n_Layers,  &
                   SensorIndex    )

          IF ( .NOT. CRTM_Predictor_Associated(Predictor) ) THEN
            Error_Status=FAILURE
            WRITE( Message,'("Error allocating predictor structure for profile #",i0, &
                   &" and ",a," sensor.")' ) m, SC(SensorIndex)%Sensor_Id
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          END IF

          ! ...And now fill them
          CALL CRTM_Compute_Predictors( SensorIndex   , &  ! Input
                                        Atmosphere(m) , &  ! Input
                                        GeometryInfo  , &  ! Input
                                        AncillaryInput, &  ! Input
                                        Predictor , &  ! Output
                                        PVar        )  ! Internal variable output
        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
            ! Compute the gas absorption
          CALL CRTM_Compute_AtmAbsorption( SensorIndex   , &  ! Input
                                           ChannelIndex  , &  ! Input
                                           AncillaryInput, &  ! Input
                                           Predictor , &  ! Input
                                           AtmOptics , &  ! Output
                                           AAvar       )  ! Internal variable output

          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status = CRTM_Compute_AerosolScatter( Atmosphere(m), &  ! Input
                                                        SensorIndex  , &  ! Input
                                                        ChannelIndex , &  ! Input
                                                        AtmOptics    , &  ! In/Output
                                                        ASVar          )  ! Internal variable output
            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF
          END IF
          ! Compute the cloud particle absorption/scattering properties
          IF( Atmosphere(m)%n_Clouds > 0 ) THEN
            Error_Status = CRTM_Compute_CloudScatter( Atmosphere(m) , &  ! Input
                                                      SensorIndex , &  ! Input
                                                      ChannelIndex, &  ! Input
                                                      AtmOptics   , &  ! Output
                                                      CSvar         )  ! Internal variable output

            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing CloudScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            END IF
          END IF
!       
          IF ( Options(m)%Use_Emissivity ) THEN
            ! ...Cloudy/all-sky case
            SfcOptics%Compute = .FALSE.
            SfcOptics%Emissivity(1,1)       = Options(m)%Emissivity(ln)
            SfcOptics%Reflectivity(1,1,1,1) = ONE - Options(m)%Emissivity(ln)
            IF ( Options(m)%Use_Direct_Reflectivity ) THEN
              SfcOptics%Direct_Reflectivity(1,1) = Options(m)%Direct_Reflectivity(ln)
            ELSE
              SfcOptics%Direct_Reflectivity(1,1) = SfcOptics%Reflectivity(1,1,1,1)
            END IF
          END IF
          ! Save the nadir optical depth
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth
!
          Error_Status = Radar_Solution( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    RTSolution(ln,m), &  ! Output
                    RTV               )  ! Internal variable output
        END DO Channel_Loop

      END DO Sensor_Loop
      ! Deallocate local sensor independent data structures
      CALL CRTM_Predictor_Destroy( Predictor )
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_SfcOptics_Destroy( SfcOptics )
      ! ...Internal variables
      CALL AOvar_Destroy( AOvar )
      CALL CSvar_Destroy( CSvar )
      CALL ASvar_Destroy( ASvar )
      CALL RTV_Destroy( RTV )
      DEALLOCATE( Height, deltaZ )
    END DO Profile_Loop

  END FUNCTION CRTM_ActiveSensor
!
!
  FUNCTION CRTM_ActiveSensor_TL( &
    Atmosphere , &  ! Input, M
    Surface    , &  ! Input, M
    Atmosphere_TL, &  ! TL  Input, M
    Surface_TL   , &  ! TL  Input, M
    Geometry   , &  ! Input, M
    ChannelInfo, &  ! Input, n_Sensors
    RTSolution , &  ! Output, L x M
    RTSolution_TL, &  ! Output, L x M
    Options    ) &  ! Optional input, M
  RESULT( Error_Status )

    ! Arguments
    TYPE(CRTM_Surface_type),           INTENT(IN)     :: Surface(:), Surface_TL(:)  ! M
    TYPE(CRTM_Geometry_type),          INTENT(IN)     :: Geometry(:)       ! M
    TYPE(CRTM_SfcOptics_type)                         :: SfcOptics, SfcOptics_TL
    TYPE(CRTM_Atmosphere_type),        INTENT(IN)     :: Atmosphere(:), Atmosphere_TL(:)     ! M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)    ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)   ! L x M
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution_TL(:,:)   ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)        ! M
    TYPE(CRTM_AncillaryInput_type) :: AncillaryInput
    TYPE(CRTM_GeometryInfo_type) :: GeometryInfo
    TYPE(CRTM_Predictor_type)    :: Predictor, Predictor_TL
    TYPE(CSvar_type)      :: CSvar  ! CloudScatter
    TYPE(AOvar_type)      :: AOvar  ! AtmOptics
    TYPE(RTV_type)        :: RTV
    ! Function result
    INTEGER :: Error_Status, Error_Status_FWD, Error_Status_TL
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_ActiveSensor_TL'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: ln, H2O_idx
    REAL(fp), Allocatable :: Height(:), deltaZ(:)
    REAL(fp) :: transmittance1, transmittance2
    ! Component variables
    TYPE(CRTM_AtmOptics_type) :: AtmOptics, AtmOptics_TL
    TYPE(CRTM_PVar_type)  :: PVar   ! Predictor
    TYPE(ASVar_type) :: ASvar
    TYPE(CRTM_AAvar_type) :: AAvar  ! AtmAbsorption
!
    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS

    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN

    ! Check the number of channels
    IF ( SIZE(RTSolution,DIM=1) < n_Channels ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure array too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    ! ...Number of atmospheric profiles.
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(RTSolution,DIM=2) /= n_Profiles ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for RTSolution argument.'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF

    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles

    CALL CRTM_SfcOptics_Create( SfcOptics  , MAX_N_ANGLES, MAX_N_STOKES )
    CALL CRTM_SfcOptics_Create( SfcOptics_TL  , MAX_N_ANGLES, MAX_N_STOKES )
    Allocate( Height(0:Atmosphere(m)%n_Layers), deltaZ(Atmosphere(m)%n_Layers) )
      H2O_idx = 1  ! default
    CALL Geopotential_Height(Atmosphere(m)%Level_Pressure      , & ! Input
                            Atmosphere(m)%Temperature         , & ! Input
                            Atmosphere(m)%Absorber(:, H2O_idx), & ! Input
                            ZERO                    , & ! Input - surface height
                            Height                ) ! Output in km
!    write(6,'(6f12.5)') Height
    deltaZ(:) = Height(0:Atmosphere(m)%n_Layers-1)-Height(1:Atmosphere(m)%n_Layers)

    CALL CRTM_RTSolution_Zero(RTSolution(:,m))
    CALL CRTM_RTSolution_Zero(RTSolution_TL(:,m))
      ! ...Compute derived geometry
    CALL CRTM_GeometryInfo_SetValue( GeometryInfo, Geometry=Geometry(m) )
    CALL CRTM_GeometryInfo_Compute( GeometryInfo )
    CALL CRTM_Compute_SurfaceT( Surface(m), SfcOptics )
    CALL CRTM_Compute_SurfaceT_TL( Surface(m),Surface_TL(m), SfcOptics_TL )

      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF
!
      CALL RTV_Create( RTV, MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, Atmosphere(m)%n_Layers )
      
      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
        ! Check whether to skip this profile
        IF ( Options(m)%Skip_Profile ) CYCLE Profile_Loop
      END IF
      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF
      ! Check the RTSolution layer dimension
      IF ( ANY(RTSolution(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! Allocate AtmOptics based on Atmosphere dimension
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS, &
                                  MAX_N_PHASE_ELEMENTS  )
      CALL CRTM_AtmOptics_Create( AtmOptics_TL, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS  , &
                                  MAX_N_PHASE_ELEMENTS   )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structure for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set the scattering switch
      AtmOptics%Include_Scattering = Options(m)%Include_Scattering
      ! ...Allocate the atmospheric optics internal structure
      CALL AOvar_Create( AOvar, Atmosphere(m)%n_Layers )

      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms = 16

      ! Allocate the scattering internal variables if necessary
      ! ...Cloud
      IF ( Atmosphere(m)%n_Clouds > 0 ) THEN
        CALL CSvar_Create( CSvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Clouds    )
      END IF
      ! ...Aerosol
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Aerosols  )
      END IF
!

      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for channel(l)/sensor(n) count
      ln = 0
      Sensor_Loop: DO n = 1, n_Sensors

        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index
        ! Compute predictors for AtmAbsorption calcs
        ! ...Allocate the predictor structures
        CALL CRTM_Predictor_Create( &
               Predictor   , &
               atmosphere(m)%n_Layers, &
               SensorIndex , &
               SaveFWV = 1   )
        CALL CRTM_Predictor_Create( &
               Predictor_TL, &
               atmosphere(m)%n_Layers, &
               SensorIndex   )
        IF ( (.NOT. CRTM_Predictor_Associated(Predictor)) .OR. &
             (.NOT. CRTM_Predictor_Associated(Predictor_TL)) ) THEN
          Error_Status=FAILURE
          WRITE( Message,'("Error allocating predictor structures for profile #",i0, &
                 &" and ",a," sensor.")' ) m, SC(SensorIndex)%Sensor_Id
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF

          IF ( .NOT. CRTM_Predictor_Associated(Predictor) ) THEN
            Error_Status=FAILURE
            WRITE( Message,'("Error allocating predictor structure for profile #",i0, &
                   &" and ",a," sensor.")' ) m, SC(SensorIndex)%Sensor_Id
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          END IF

          ! ...And now fill them
          CALL CRTM_Compute_Predictors( SensorIndex   , &  ! Input
                                        Atmosphere(m) , &  ! Input
                                        GeometryInfo  , &  ! Input
                                        AncillaryInput, &  ! Input
                                        Predictor , &  ! Output
                                        PVar        )  ! Internal variable output

          CALL CRTM_Compute_Predictors_TL( SensorIndex   , &  ! Input
                                         Atmosphere(m) , &  ! Input
                                         Predictor     , &  ! Input
                                         Atmosphere_TL(m), &  ! Input
                                         AncillaryInput, &  ! Input
                                         Predictor_TL  , &  ! Output
                                         PVar           )  ! Internal variable input
        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels
          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
          CALL CRTM_AtmOptics_Zero( AtmOptics_TL )
            ! Compute the gas absorption
          CALL CRTM_Compute_AtmAbsorption( SensorIndex   , &  ! Input
                                           ChannelIndex  , &  ! Input
                                           AncillaryInput, &  ! Input
                                           Predictor , &  ! Input
                                           AtmOptics , &  ! Output
                                           AAvar       )  ! Internal variable output
          CALL CRTM_Compute_AtmAbsorption_TL( SensorIndex     , &  ! Input
                                              ChannelIndex    , &  ! Input
                                              Predictor       , &  ! Input
                                              Predictor_TL    , &  ! Input
                                              AtmOptics_TL    , &  ! Output
                                              AAvar            )  ! Internal variable input

          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status_FWD = CRTM_Compute_AerosolScatter( Atmosphere(m), &  ! Input
                                                        SensorIndex  , &  ! Input
                                                        ChannelIndex , &  ! Input
                                                        AtmOptics    , &  ! In/Output
                                                        ASVar          )  ! Internal variable output
!
            Error_Status_TL  = CRTM_Compute_AerosolScatter_TL( Atmosphere(m) , &  ! FWD Input
                                                         AtmOptics   , &  ! FWD Input
                                                         Atmosphere_TL(m), &  ! TL  Input
                                                         SensorIndex , &  ! Input
                                                         ChannelIndex, &  ! Input
                                                         AtmOptics_TL, &  ! TL  Output
                                                         ASvar        )  ! Internal variable input
!
            IF ( Error_Status_FWD /= SUCCESS .OR. Error_Status_TL /= SUCCESS) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              CYCLE  !RETURN
            END IF
          END IF

          ! Compute the cloud particle absorption/scattering properties
          IF( Atmosphere(m)%n_Clouds > 0 ) THEN
            Error_Status_FWD = CRTM_Compute_CloudScatter( Atmosphere(m) , &  ! Input
                                                      SensorIndex , &  ! Input
                                                      ChannelIndex, &  ! Input
                                                      AtmOptics   , &  ! Output
                                                      CSvar         )  ! Internal variable output
            Error_Status_TL = CRTM_Compute_CloudScatter_TL( Atmosphere(m), &  ! FWD Input
                                                      AtmOptics   , &  ! FWD Input
                                                      Atmosphere_TL(m) , &  ! TL  Input
!                                                      GeometryInfo, &  ! Input
                                                      SensorIndex , &  ! Input
                                                      ChannelIndex, &  ! Input
                                                      AtmOptics_TL, &  ! TL  Output
                                                      CSvar        )  ! Internal variable input
            IF ( Error_Status_FWD /= SUCCESS .OR. Error_Status_TL /= SUCCESS) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing CloudScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              CYCLE  !RETURN
            END IF
          END IF
!
          IF ( Options(m)%Use_Emissivity ) THEN
            ! ...Cloudy/all-sky case
            SfcOptics%Compute = .FALSE.
            SfcOptics%Emissivity(1,1)       = Options(m)%Emissivity(ln)
            SfcOptics%Reflectivity(1,1,1,1) = ONE - Options(m)%Emissivity(ln)
            IF ( Options(m)%Use_Direct_Reflectivity ) THEN
              SfcOptics%Direct_Reflectivity(1,1) = Options(m)%Direct_Reflectivity(ln)
            ELSE
              SfcOptics%Direct_Reflectivity(1,1) = SfcOptics%Reflectivity(1,1,1,1)
            END IF
          END IF
          ! Save the nadir optical depth
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth
!
          Error_Status_FWD = Radar_Solution( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    RTSolution(ln,m), &  ! Output
                    RTV               )  ! Internal variable output
!
          Error_Status_FWD = Radar_Solution_TL( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    Atmosphere_TL(m)   , &  ! Input
                    Surface_TL(m)      , &  ! Input
                    AtmOptics_TL    , &  ! Input
                    SfcOptics_TL    , &  ! Input                    
                    RTSolution(ln,m), &  ! Output
                    RTSolution_TL(ln,m), &  ! Output                
                    RTV               )  ! Internal variable output
        END DO Channel_Loop

      END DO Sensor_Loop

      ! Deallocate local sensor independent data structures
      CALL CRTM_Predictor_Destroy( Predictor )
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_SfcOptics_Destroy( SfcOptics )
      ! ...Internal variables
      CALL AOvar_Destroy( AOvar )
      CALL CSvar_Destroy( CSvar )
      CALL ASvar_Destroy( ASvar )
      CALL RTV_Destroy( RTV )
      DEALLOCATE( Height, deltaZ )
    END DO Profile_Loop

  END FUNCTION CRTM_ActiveSensor_TL
!
!
  FUNCTION CRTM_ActiveSensor_AD( &
    Atmosphere , &  ! Input, M
    Surface    , &  ! Input, M
    RTSolution_AD, &  ! AD  Input, L x M
    Geometry   , &  ! Input, M
    ChannelInfo, &  ! Input, n_Sensors
    Atmosphere_AD, &  ! AD  Output, M
    Surface_AD   , &  ! AD  Output, M
    RTSolution , &  ! Output, L x M
    Options    ) &  ! Optional input, M
  RESULT( Error_Status )

    ! Arguments
    TYPE(CRTM_Surface_type),           INTENT(IN)     :: Surface(:) ! M
    TYPE(CRTM_Surface_type),           INTENT(IN OUT) :: Surface_AD(:)  ! M
    TYPE(CRTM_Geometry_type),          INTENT(IN)     :: Geometry(:)       ! M
    TYPE(CRTM_SfcOptics_type)                         :: SfcOptics
    TYPE(CRTM_SfcOptics_type)                         :: SfcOptics_AD
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT)     :: Atmosphere(:)  ! M
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT)     :: Atmosphere_AD(:)     ! M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)    ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)   ! L x M
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution_AD(:,:)   ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)        ! M
    TYPE(CRTM_AncillaryInput_type) :: AncillaryInput
    TYPE(CRTM_GeometryInfo_type) :: GeometryInfo
    TYPE(CRTM_Predictor_type)    :: Predictor, Predictor_AD
    TYPE(CSvar_type)      :: CSvar  ! CloudScatter
    TYPE(AOvar_type)      :: AOvar  ! AtmOptics
    TYPE(RTV_type)        :: RTV
    ! Function result
    INTEGER :: Error_Status, Error_Status_FWD, Error_Status_AD
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_ActiveSensor_TL'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: Check_Input
    INTEGER :: n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles
    INTEGER :: ln, H2O_idx
    REAL(fp), Allocatable :: Height(:), deltaZ(:)
    REAL(fp) :: transmittance1, transmittance2
    ! Component variables
    TYPE(CRTM_AtmOptics_type) :: AtmOptics, AtmOptics_AD
    TYPE(CRTM_PVar_type)  :: PVar   ! Predictor
    TYPE(ASVar_type) :: ASvar
    TYPE(CRTM_AAvar_type) :: AAvar  ! AtmAbsorption
!
    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS

    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN
    ! Check the number of channels
    IF ( SIZE(RTSolution,DIM=1) < n_Channels ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure array too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF


    ! Check the number of profiles
    ! ...Number of atmospheric profiles.
    n_Profiles = SIZE(Atmosphere)
    ! ...Check the profile dimensionality of the other mandatory arguments
    IF ( SIZE(RTSolution,DIM=2) /= n_Profiles ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for RTSolution argument.'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF

    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles

    CALL CRTM_SfcOptics_Create( SfcOptics  , MAX_N_ANGLES, MAX_N_STOKES )
    CALL CRTM_SfcOptics_Create( SfcOptics_AD  , MAX_N_ANGLES, MAX_N_STOKES )
    Allocate( Height(0:Atmosphere(m)%n_Layers), deltaZ(Atmosphere(m)%n_Layers) )
      H2O_idx = 1  ! default
    CALL Geopotential_Height(Atmosphere(m)%Level_Pressure      , & ! Input
                            Atmosphere(m)%Temperature         , & ! Input
                            Atmosphere(m)%Absorber(:, H2O_idx), & ! Input
                            ZERO                    , & ! Input - surface height
                            Height                ) ! Output in km
!    write(6,'(6f12.5)') Height
    deltaZ(:) = Height(0:Atmosphere(m)%n_Layers-1)-Height(1:Atmosphere(m)%n_Layers)

    CALL CRTM_RTSolution_Zero(RTSolution(:,m))
!    CALL CRTM_RTSolution_Zero(RTSolution_TL(:,m))
      ! ...Compute derived geometry
    CALL CRTM_GeometryInfo_SetValue( GeometryInfo, Geometry=Geometry(m) )
    CALL CRTM_GeometryInfo_Compute( GeometryInfo )
    CALL CRTM_Compute_SurfaceT( Surface(m), SfcOptics )

      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF
!
      CALL RTV_Create( RTV, MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, Atmosphere(m)%n_Layers )

      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
        ! Check whether to skip this profile
        IF ( Options(m)%Skip_Profile ) CYCLE Profile_Loop
      END IF
      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF
      ! Check the RTSolution layer dimension
      IF ( ANY(RTSolution(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF

      ! Allocate AtmOptics based on Atmosphere dimension
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers, &
                                  MAX_N_LEGENDRE_TERMS, &
                                  MAX_N_PHASE_ELEMENTS  )
      CALL CRTM_AtmOptics_Create( AtmOptics_AD, &
                                  Atmosphere(m)%n_Layers        , &
                                  MAX_N_LEGENDRE_TERMS, &
                                  MAX_N_PHASE_ELEMENTS  )
      !AtmOptics%Optical_Depth  = ZERO
      !AtmOptics_AD%Optical_Depth = ZERO      
!
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structure for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
      ! ...Set the scattering switch
      AtmOptics%Include_Scattering = Options(m)%Include_Scattering
      ! ...Allocate the atmospheric optics internal structure
      CALL AOvar_Create( AOvar, Atmosphere(m)%n_Layers )

      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms = 16

      ! Allocate the scattering internal variables if necessary
      ! ...Cloud
      IF ( Atmosphere(m)%n_Clouds > 0 ) THEN
        CALL CSvar_Create( CSvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Clouds    )
      END IF
      ! ...Aerosol
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Aerosols  )
      END IF
!
      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for channel(l)/sensor(n) count
      ln = 0
      Sensor_Loop: DO n = 1, n_Sensors
        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index
        ! Allocate the AtmAbsorption predictor structures
        CALL CRTM_Predictor_Create( &
               Predictor   , &
               atmosphere(m)%n_Layers, &
               SensorIndex , &
               SaveFWV = 1   )
        CALL CRTM_Predictor_Create( &
               Predictor_AD, &
               atmosphere(m)%n_Layers, &
               SensorIndex   )
        IF ( (.NOT. CRTM_Predictor_Associated(Predictor)) .OR. &
             (.NOT. CRTM_Predictor_Associated(Predictor_AD)) ) THEN
          Error_Status=FAILURE
          WRITE( Message,'("Error allocating predictor structures for profile #",i0, &
                 &" and ",a," sensor.")' ) m, SC(SensorIndex)%Sensor_Id
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF


          ! ...And now fill them
          CALL CRTM_Compute_Predictors( SensorIndex   , &  ! Input
                                        Atmosphere(m) , &  ! Input
                                        GeometryInfo  , &  ! Input
                                        AncillaryInput, &  ! Input
                                        Predictor , &  ! Output
                                        PVar        )  ! Internal variable output
        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)


          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
          CALL CRTM_AtmOptics_Zero( AtmOptics_AD )
            ! Compute the gas absorption
          CALL CRTM_Compute_AtmAbsorption( SensorIndex   , &  ! Input
                                           ChannelIndex  , &  ! Input
                                           AncillaryInput, &  ! Input
                                           Predictor , &  ! Input
                                           AtmOptics , &  ! Output
                                           AAvar       )  ! Internal variable output

          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status_FWD = CRTM_Compute_AerosolScatter( Atmosphere(m), &  ! Input
                                                        SensorIndex  , &  ! Input
                                                        ChannelIndex , &  ! Input
                                                        AtmOptics    , &  ! In/Output
                                                        ASVar          )  ! Internal variable output
!
            IF ( Error_Status_FWD /= SUCCESS ) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              CYCLE  !RETURN
            END IF
          END IF

          ! Compute the cloud particle absorption/scattering properties
          IF( Atmosphere(m)%n_Clouds > 0 ) THEN
            Error_Status_FWD = CRTM_Compute_CloudScatter( Atmosphere(m) , &  ! Input
                                                      SensorIndex , &  ! Input
                                                      ChannelIndex, &  ! Input
                                                      AtmOptics   , &  ! Output
                                                      CSvar         )  ! Internal variable output
!
            IF ( Error_Status_FWD /= SUCCESS ) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing CloudScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              CYCLE  !RETURN
            END IF
          END IF
!
          IF ( Options(m)%Use_Emissivity ) THEN
            ! ...Cloudy/all-sky case
            SfcOptics%Compute = .FALSE.
            SfcOptics%Emissivity(1,1)       = Options(m)%Emissivity(ln)
            SfcOptics%Reflectivity(1,1,1,1) = ONE - Options(m)%Emissivity(ln)
            IF ( Options(m)%Use_Direct_Reflectivity ) THEN
              SfcOptics%Direct_Reflectivity(1,1) = Options(m)%Direct_Reflectivity(ln)
            ELSE
              SfcOptics%Direct_Reflectivity(1,1) = SfcOptics%Reflectivity(1,1,1,1)
            END IF
          END IF
          ! Save the nadir optical depth
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth
!
          Error_Status_FWD = Radar_Solution( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    RTSolution(ln,m), &  ! Output
                    RTV               )  ! Internal variable output
!
          Error_Status_AD = Radar_Solution_AD( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    RTSolution_AD(ln,m), &  ! Output
                    Atmosphere_AD(m)   , &  ! Input
                    Surface_AD(m)      , &  ! Input
                    AtmOptics_AD    , &  ! Input
                    SfcOptics_AD    , &  ! Input                    
                    RTSolution(ln,m), &  ! Output
                    RTV               )  ! Internal variable output
!
            IF ( Error_Status_AD /= SUCCESS ) THEN
              WRITE( Message,'( "Error computing Radar_Solution_ADD for ", a, &
                     &", channel ", i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status_AD )
              RETURN
            END IF
          ! Compute the adjoint cloud absorption/scattering properties
          IF ( Atmosphere(m)%n_Clouds > 0 ) THEN
            Error_Status_AD = CRTM_Compute_CloudScatter_AD( Atmosphere(m) , &  ! FWD Input
                                                         AtmOptics   , &  ! FWD Input
                                                         AtmOptics_AD, &  ! AD  Input
                                                         SensorIndex , &  ! Input
                                                         ChannelIndex, &  ! Input
                                                         Atmosphere_AD(m) , &  ! AD  Output
                                                         CSvar         )  ! Internal variable input
            IF ( Error_Status_AD /= SUCCESS ) THEN
              WRITE( Message,'("Error computing CloudScatter_AD for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status_AD )
              RETURN
            END IF
          END IF
!
          ! Compute the adjoint aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status_AD = CRTM_Compute_AerosolScatter_AD( Atmosphere(m) , &  ! FWD Input
                                                           AtmOptics   , &  ! FWD Input
                                                           AtmOptics_AD, &  ! AD  Input
                                                           SensorIndex , &  ! Input
                                                           ChannelIndex, &  ! Input
                                                           Atmosphere_AD(m) , &  ! AD  Output
                                                           ASvar         )  ! Internal variable input
            IF ( Error_Status_AD /= SUCCESS ) THEN
              WRITE( Message,'("Error computing AerosolScatter_AD for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status_AD )
              RETURN
            END IF
          END IF
!
          ! Compute the adjoint molecular scattering properties
!!          IF( RTV%Visible_Flag_true ) THEN
!!            Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex)
!!            Error_Status = CRTM_Compute_MoleculeScatter_AD( &
!!                             Wavenumber  , &
!!                             AtmOptics_AD, &
!!                             Atm_AD        )
!!            IF ( Error_Status /= SUCCESS ) THEN
!!              WRITE( Message,'("Error computing MoleculeScatter_AD for ",a,&
!!                     &", channel ",i0,", profile #",i0)' ) &
!!                     TRIM(ChannelInfo(n)%Sensor_ID), &
!!                     ChannelInfo(n)%Sensor_Channel(l), &
!!                     m
!!              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
!!              RETURN
!!            END IF
!!          END IF
!
          ! Compute the adjoint gaseous absorption
          CALL CRTM_Compute_AtmAbsorption_AD( SensorIndex     , &  ! Input
                                              ChannelIndex    , &  ! Input
                                              Predictor       , &  ! FWD Input
                                              AtmOptics_AD    , &  ! AD  Input
                                              Predictor_AD    , &  ! AD  Output
                                              AAVar             )  ! Internal variable input

        END DO Channel_Loop
        ! Adjoint of the predictor calculations
        CALL CRTM_Compute_Predictors_AD( SensorIndex   , &  ! Input
                                         Atmosphere(m) , &  ! FWD Input
                                         Predictor     , &  ! FWD Input
                                         Predictor_AD  , &  ! AD  Input
                                         AncillaryInput, &  ! Input
                                         Atmosphere_AD(m), &  ! AD  Output
                                         PVar            )  ! Internal variable input


        ! Deallocate local sensor dependent data structures
        ! ...RTV structure
        IF ( RTV_Associated(RTV) ) CALL RTV_Destroy(RTV)
        ! ...Predictor structures
        CALL CRTM_Predictor_Destroy( Predictor )
        CALL CRTM_Predictor_Destroy( Predictor_AD )

      END DO Sensor_Loop

      ! ...Adjoint of average surface skin temperature for multi-surface types
      CALL CRTM_Compute_SurfaceT_AD( Surface(m), SfcOptics_AD, Surface_AD(m) )
      
      ! Deallocate local sensor independent data structures
      CALL CRTM_Predictor_Destroy( Predictor )
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_SfcOptics_Destroy( SfcOptics )
      ! ...Internal variables
      CALL AOvar_Destroy( AOvar )
      CALL CSvar_Destroy( CSvar )
      CALL ASvar_Destroy( ASvar )
      CALL RTV_Destroy( RTV )
      DEALLOCATE( Height, deltaZ )
    END DO Profile_Loop

  END FUNCTION CRTM_ActiveSensor_AD
!
!
  FUNCTION CRTM_ActiveSensor_K( &
    Atmosphere , &  ! Input, M
    Surface    , &  ! Input, M
    RTSolution_K, &  ! AD  Input, L x M
    Geometry   , &  ! Input, M
    ChannelInfo, &  ! Input, n_Sensors
    Atmosphere_K, &  ! AD  Output, M
    Surface_K  , &  ! AD  Output, M
    RTSolution , &  ! Output, L x M
    Options    ) &  ! Optional input, M
  RESULT( Error_Status )

    ! Arguments
    TYPE(CRTM_Surface_type),           INTENT(IN)     :: Surface(:) ! M
    TYPE(CRTM_Surface_type),           INTENT(IN OUT) :: Surface_K(:,:)  ! M
    TYPE(CRTM_Geometry_type),          INTENT(IN)     :: Geometry(:)       ! M
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT) :: Atmosphere(:)  ! M
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT) :: Atmosphere_K(:,:)     ! M
    TYPE(CRTM_ChannelInfo_type),       INTENT(IN)     :: ChannelInfo(:)    ! n_Sensors
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:,:)   ! L x M
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution_K(:,:)   ! L x M
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options(:)        ! M
    TYPE(CRTM_AncillaryInput_type) :: AncillaryInput
    TYPE(CRTM_GeometryInfo_type) :: GeometryInfo
    TYPE(CRTM_Predictor_type)    :: Predictor, Predictor_K
    TYPE(CRTM_SfcOptics_type)    :: SfcOptics, SfcOptics_K
    TYPE(CRTM_AtmOptics_type)    :: AtmOptics, AtmOptics_K
    TYPE(CSvar_type)      :: CSvar  ! CloudScatter
    TYPE(AOvar_type)      :: AOvar  ! AtmOptics
    TYPE(RTV_type)        :: RTV
    ! Function result
    INTEGER :: Error_Status, Error_Status_FWD, Error_Status_K
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_ActiveSensor_TL'
    ! Local variables
    CHARACTER(ML) :: Message
    LOGICAL :: Options_Present
    LOGICAL :: User_Emissivity, User_Direct_Reflectivity, User_N_Streams
    LOGICAL :: Check_Input
    INTEGER :: j, n, n_Sensors,  SensorIndex
    INTEGER :: l, n_Channels, ChannelIndex
    INTEGER :: m, n_Profiles, nc, na
    INTEGER :: ln, H2O_idx
    INTEGER :: RT_Algorithm_Id
    REAL(fp), Allocatable :: Height(:), deltaZ(:)
    REAL(fp) :: transmittance1, transmittance2
    ! Component variables
    TYPE(CRTM_Options_type) :: Default_Options
    TYPE(CRTM_PVar_type)  :: PVar   ! Predictor
    TYPE(ASVar_type) :: ASvar
    TYPE(CRTM_AAvar_type) :: AAvar  ! AtmAbsorption
!
    ! ------
    ! SET UP
    ! ------
    Error_Status = SUCCESS

    ! If no sensors or channels, simply return
    n_Sensors  = SIZE(ChannelInfo)
    n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    IF ( n_Sensors == 0 .OR. n_Channels == 0 ) RETURN


    ! Check the number of channels
    IF ( SIZE(RTSolution,DIM=1) < n_Channels ) THEN
      Error_Status = FAILURE
      WRITE( Message,'("Output RTSolution structure array too small (",i0,&
             &") to hold results for the number of requested channels (",i0,")")') &
             SIZE(RTSolution,DIM=1), n_Channels
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF

    ! Check the number of profiles
    ! ...Number of atmospheric profiles.
    n_Profiles = SIZE(Atmosphere)
    IF ( SIZE(Surface)            /= n_Profiles .OR. &
         SIZE(RTSolution_K,DIM=2) /= n_Profiles .OR. &
         SIZE(Geometry)           /= n_Profiles .OR. &
         SIZE(Atmosphere_K,DIM=2) /= n_Profiles .OR. &
         SIZE(Surface_K   ,DIM=2) /= n_Profiles .OR. &
         SIZE(RTSolution  ,DIM=2) /= n_Profiles      ) THEN
      Error_Status = FAILURE
      Message = 'Inconsistent profile dimensionality for input arguments.'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ...Check the profile dimensionality of the other optional arguments
    Options_Present = PRESENT(Options)
    IF ( Options_Present ) THEN
      IF ( SIZE(Options) /= n_Profiles ) THEN
        Error_Status = FAILURE
        Message = 'Inconsistent profile dimensionality for Options optional input argument.'
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF
    END IF

    ! Allocate the profile independent surface optics local structure
    CALL CRTM_SfcOptics_Create( SfcOptics  , MAX_N_ANGLES, MAX_N_STOKES )
    CALL CRTM_SfcOptics_Create( SfcOptics_K, MAX_N_ANGLES, MAX_N_STOKES )
    IF ( (.NOT. CRTM_SfcOptics_Associated(SfcOptics  )) .OR. &
         (.NOT. CRTM_SfcOptics_Associated(SfcOptics_K)) ) THEN
      Error_Status = FAILURE
      Message = 'Error allocating SfcOptics data structures'
      CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
      RETURN
    END IF
    ! ------------
    ! PROFILE LOOP
    ! ------------
    Profile_Loop: DO m = 1, n_Profiles
    
    Allocate( Height(0:Atmosphere(m)%n_Layers), deltaZ(Atmosphere(m)%n_Layers) )
      H2O_idx = 1  ! default
    CALL Geopotential_Height(Atmosphere(m)%Level_Pressure      , & ! Input
                            Atmosphere(m)%Temperature         , & ! Input
                            Atmosphere(m)%Absorber(:, H2O_idx), & ! Input
                            ZERO                    , & ! Input - surface height
                            Height                ) ! Output in km
!    write(6,'(6f12.5)') Height
    deltaZ(:) = Height(0:Atmosphere(m)%n_Layers-1)-Height(1:Atmosphere(m)%n_Layers)

    CALL CRTM_RTSolution_Zero(RTSolution(:,m))
!    CALL CRTM_RTSolution_Zero(RTSolution_TL(:,m))
      ! ...Compute derived geometry
    CALL CRTM_GeometryInfo_SetValue( GeometryInfo, Geometry=Geometry(m) )
    CALL CRTM_GeometryInfo_Compute( GeometryInfo )
    CALL CRTM_Compute_SurfaceT( Surface(m), SfcOptics )

      ! Check the aerosol coeff. data for cases with aerosols
      IF( Atmosphere(m)%n_Aerosols > 0 .AND. .NOT. CRTM_AerosolCoeff_IsLoaded() )THEN
         Error_Status = FAILURE
         WRITE( Message,'("The AerosolCoeff data must be loaded (with CRTM_Init routine) ", &
                &"for the aerosol case profile #",i0)' ) m
         CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
         RETURN
      END IF
!
      CALL RTV_Create( RTV, MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, Atmosphere(m)%n_Layers )
!

      ! Copy over forward "non-variable" inputs to K-matrix outputs
      DO l = 1, n_Channels
        ! ...Atmosphere
        Atmosphere_K(l,m)%Climatology = Atmosphere(m)%Climatology
        ! Loop over absorbers
        DO j = 1, Atmosphere(m)%n_Absorbers
          Atmosphere_K(l,m)%Absorber_ID(j)    = Atmosphere(m)%Absorber_ID(j)
          Atmosphere_K(l,m)%Absorber_Units(j) = Atmosphere(m)%Absorber_Units(j)
        END DO
        ! Loop over and assign cloud types
        DO nc = 1, Atmosphere(m)%n_Clouds
          Atmosphere_K(l,m)%Cloud(nc)%Type = Atmosphere(m)%Cloud(nc)%Type
        END DO
        ! Loop over and assign aerosol types
        DO na = 1, Atmosphere(m)%n_Aerosols
          Atmosphere_K(l,m)%Aerosol(na)%Type = Atmosphere(m)%Aerosol(na)%Type
        END DO
        ! ...Surface
        Surface_K(l,m)%Land_Coverage  = Surface(m)%Land_Coverage
        Surface_K(l,m)%Water_Coverage = Surface(m)%Water_Coverage
        Surface_K(l,m)%Snow_Coverage  = Surface(m)%Snow_Coverage
        Surface_K(l,m)%Ice_Coverage   = Surface(m)%Ice_Coverage
        Surface_K(l,m)%Land_Type  = Surface(m)%Land_Type
        Surface_K(l,m)%Water_Type = Surface(m)%Water_Type
        Surface_K(l,m)%Snow_Type  = Surface(m)%Snow_Type
        Surface_K(l,m)%Ice_Type   = Surface(m)%Ice_Type
      END DO

      ! Check the optional Options structure argument
      Check_Input = .TRUE.
      IF (Options_Present) THEN
        Check_Input = Options(m)%Check_Input
        ! Check whether to skip this profile
        IF ( Options(m)%Skip_Profile ) CYCLE Profile_Loop
      END IF
      ! Check the input atmosphere if required
      IF ( Check_Input ) THEN
        IF ( .NOT. CRTM_Atmosphere_IsValid( Atmosphere(m) ) ) THEN
          Error_Status = FAILURE
          WRITE( Message,'("Input data check failed for profile #",i0)' ) m
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
      END IF
!
      User_Emissivity       = Default_Options%Use_Emissivity
      RT_Algorithm_Id       = Default_Options%RT_Algorithm_Id
      User_N_Streams        = Default_Options%Use_N_Streams
      ! ...Check the Options argument
      IF (Options_Present) THEN
        ! Override input checker with option
        Check_Input = Options(m)%Check_Input
        ! Check if the supplied emissivity should be used
        User_Emissivity = Options(m)%Use_Emissivity
        IF ( Options(m)%Use_Emissivity ) THEN
          ! Are the channel dimensions consistent
          IF ( Options(m)%n_Channels < n_Channels ) THEN
            Error_Status = FAILURE
            WRITE( Message,'( "Input Options channel dimension (", i0, ") is less ", &
                   &"than the number of requested channels (",i0, ")" )' ) &
                   Options(m)%n_Channels, n_Channels
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
          ! Check if the supplied direct reflectivity should be used
          User_Direct_Reflectivity = Options(m)%Use_Direct_Reflectivity
        END IF
        ! Copy over surface optics input
        SfcOptics%Use_New_MWSSEM = .NOT. Options(m)%Use_Old_MWSSEM
        ! Specify the RT algorithm
        RT_Algorithm_Id = Options(m)%RT_Algorithm_ID
        ! Check if n_Streams should be used
        User_N_Streams = Options(m)%Use_N_Streams
        IF ( User_N_Streams ) THEN
          IF ( Options(m)%n_Streams <= 0 .OR. MOD(Options(m)%n_Streams,2) /= 0 .OR. &
               Options(m)%n_Streams > MAX_N_STREAMS ) THEN
              Error_Status = FAILURE
              WRITE( Message,'( "Input Options n_Streams (", i0, ") is invalid" )' ) &
                     Options(m)%n_Streams
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
          END IF
        END IF
      END IF

      ! Check the RTSolution layer dimension
      IF ( ANY(RTSolution(:,m)%n_Layers < Atmosphere(m)%n_Layers) ) THEN
        Error_Status=FAILURE
        WRITE( Message,'("Number of RTSolution layers < Atmosphere for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF

      ! ...Allocate the atmospheric optics structures based on Atm extension
      CALL CRTM_AtmOptics_Create( AtmOptics, &
                                  Atmosphere(m)%n_Layers        , &
                                  MAX_N_LEGENDRE_TERMS, &
                                  MAX_N_PHASE_ELEMENTS  )
      CALL CRTM_AtmOptics_Create( AtmOptics_K, &
                                  Atmosphere(m)%n_Layers        , &
                                  MAX_N_LEGENDRE_TERMS, &
                                  MAX_N_PHASE_ELEMENTS  )
      IF ( .NOT. CRTM_AtmOptics_Associated( Atmoptics ) .OR. &
           .NOT. CRTM_AtmOptics_Associated( Atmoptics_K ) ) THEN
        Error_Status = FAILURE
        WRITE( Message,'("Error allocating AtmOptics data structures for profile #",i0)' ) m
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN
      END IF

      ! ...Set the scattering switch
      AtmOptics%Include_Scattering = Options(m)%Include_Scattering
      ! ...Allocate the atmospheric optics internal structure
      CALL AOvar_Create( AOvar, Atmosphere(m)%n_Layers )

      ! ...Set default number of streams
      AtmOptics%n_Legendre_Terms = 16

      ! Allocate the scattering internal variables if necessary
      ! ...Cloud
      IF ( Atmosphere(m)%n_Clouds > 0 ) THEN
        CALL CSvar_Create( CSvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Clouds    )
      END IF
      ! ...Aerosol
      IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
        CALL ASvar_Create( ASvar, &
                           MAX_N_LEGENDRE_TERMS   , &
                           MAX_N_PHASE_ELEMENTS   , &
                           Atmosphere(m)%n_Layers , &
                           Atmosphere(m)%n_Aerosols  )
      END IF
!
      ! -----------
      ! SENSOR LOOP
      ! -----------
      ! Initialise channel counter for channel(l)/sensor(n) count
      ln = 0
      Sensor_Loop: DO n = 1, n_Sensors

        ! Shorter name
        SensorIndex = ChannelInfo(n)%Sensor_Index
        ! Allocate the AtmAbsorption predictor structures
        CALL CRTM_Predictor_Create( &
               Predictor   , &
               atmosphere(m)%n_Layers, &
               SensorIndex , &
               SaveFWV = 1   )
        CALL CRTM_Predictor_Create( &
               Predictor_K , &
               atmosphere(m)%n_Layers, &
               SensorIndex   )
        IF ( (.NOT. CRTM_Predictor_Associated(Predictor)) .OR. &
             (.NOT. CRTM_Predictor_Associated(Predictor_K)) ) THEN
          Error_Status=FAILURE
          WRITE( Message,'("Error allocating predictor structures for profile #",i0, &
                 &" and ",a," sensor.")' ) m, SC(SensorIndex)%Sensor_Id
          CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
          RETURN
        END IF
          ! ...And now fill them
          CALL CRTM_Compute_Predictors( SensorIndex   , &  ! Input
                                        Atmosphere(m) , &  ! Input
                                        GeometryInfo  , &  ! Input
                                        AncillaryInput, &  ! Input
                                        Predictor , &  ! Output
                                        PVar        )  ! Internal variable output
        ! ------------
        ! CHANNEL LOOP
        ! ------------
        Channel_Loop: DO l = 1, ChannelInfo(n)%n_Channels

          ! Channel setup
          ! ...Skip channel if requested
          IF ( .NOT. ChannelInfo(n)%Process_Channel(l) ) CYCLE Channel_Loop
          ! ...Shorter name
          ChannelIndex = ChannelInfo(n)%Channel_Index(l)
          ! ...Increment the processed channel counter
          ln = ln + 1
          ! ...Assign sensor+channel information to output
          RTSolution(ln,m)%Sensor_Id        = ChannelInfo(n)%Sensor_Id
          RTSolution(ln,m)%WMO_Satellite_Id = ChannelInfo(n)%WMO_Satellite_Id
          RTSolution(ln,m)%WMO_Sensor_Id    = ChannelInfo(n)%WMO_Sensor_Id
          RTSolution(ln,m)%Sensor_Channel   = ChannelInfo(n)%Sensor_Channel(l)

          ! Initialisations
          CALL CRTM_AtmOptics_Zero( AtmOptics )
          CALL CRTM_AtmOptics_Zero( AtmOptics_K )
            ! Compute the gas absorption
          CALL CRTM_Compute_AtmAbsorption( SensorIndex   , &  ! Input
                                           ChannelIndex  , &  ! Input
                                           AncillaryInput, &  ! Input
                                           Predictor , &  ! Input
                                           AtmOptics , &  ! Output
                                           AAvar       )  ! Internal variable output

          ! Compute the aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status_FWD = CRTM_Compute_AerosolScatter( Atmosphere(m), &  ! Input
                                                        SensorIndex  , &  ! Input
                                                        ChannelIndex , &  ! Input
                                                        AtmOptics    , &  ! In/Output
                                                        ASVar          )  ! Internal variable output
!
            IF ( Error_Status_FWD /= SUCCESS ) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing AerosolScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              CYCLE  !RETURN
            END IF
          END IF
          ! Compute the cloud particle absorption/scattering properties
          IF( Atmosphere(m)%n_Clouds > 0 ) THEN
            Error_Status_FWD = CRTM_Compute_CloudScatter( Atmosphere(m) , &  ! Input
                                                      SensorIndex , &  ! Input
                                                      ChannelIndex, &  ! Input
                                                      AtmOptics   , &  ! Output
                                                      CSvar         )  ! Internal variable output
            IF ( Error_Status_FWD /= SUCCESS ) THEN
              Error_Status = FAILURE
              WRITE( Message,'("Error computing CloudScatter for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              CYCLE  !RETURN
            END IF
          END IF
!
          IF ( Options(m)%Use_Emissivity ) THEN
            ! ...Cloudy/all-sky case
            SfcOptics%Compute = .FALSE.
            SfcOptics%Emissivity(1,1)       = Options(m)%Emissivity(ln)
            SfcOptics%Reflectivity(1,1,1,1) = ONE - Options(m)%Emissivity(ln)
            IF ( Options(m)%Use_Direct_Reflectivity ) THEN
              SfcOptics%Direct_Reflectivity(1,1) = Options(m)%Direct_Reflectivity(ln)
            ELSE
              SfcOptics%Direct_Reflectivity(1,1) = SfcOptics%Reflectivity(1,1,1,1)
            END IF
          END IF
          ! Save the nadir optical depth
          RTSolution(ln,m)%Layer_Optical_Depth(1:Atmosphere(m)%n_Layers) = AtmOptics%Optical_Depth
!
          Error_Status_FWD = Radar_Solution( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    RTSolution(ln,m), &  ! Output
                    RTV               )  ! Internal variable output
!
!
          Error_Status_K = Radar_Solution_AD( &
                    Atmosphere(m)   , &  ! Input
                    Surface(m)      , &  ! Input
                    AtmOptics       , &  ! Input
                    SfcOptics       , &  ! Input
                    GeometryInfo    , &  ! Input
                    SensorIndex     , &  ! Input
                    ChannelIndex    , &  ! Input
                    deltaZ          , &  ! Input
                    SC(SensorIndex)%Wavenumber(ChannelIndex), & ! Input
                    RTSolution_K(ln,m), &  ! Output
                    Atmosphere_K(ln,m)   , &  ! Input
                    Surface_K(ln,m)      , &  ! Input
                    AtmOptics_K    , &  ! Input
                    SfcOptics_K    , &  ! Input                    
                    RTSolution(ln,m), &  ! Output
                    RTV               )  ! Internal variable output
!
            IF ( Error_Status_K /= SUCCESS ) THEN
              WRITE( Message,'( "Error computing Radar_Solution_ADD for ", a, &
                     &", channel ", i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status_K )
              RETURN
            END IF
          ! Compute the adjoint cloud absorption/scattering properties
          IF ( Atmosphere(m)%n_Clouds > 0 ) THEN
            Error_Status_K = CRTM_Compute_CloudScatter_AD( Atmosphere(m) , &  ! FWD Input
                                                         AtmOptics   , &  ! FWD Input
                                                         AtmOptics_K, &  ! AD  Input
                                                         SensorIndex , &  ! Input
                                                         ChannelIndex, &  ! Input
                                                         Atmosphere_K(ln,m) , &  ! AD  Output
                                                         CSvar         )  ! Internal variable input
            IF ( Error_Status_K /= SUCCESS ) THEN
              WRITE( Message,'("Error computing CloudScatter_AD for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status_K )
              RETURN
            END IF
          END IF
!
          ! Compute the adjoint aerosol absorption/scattering properties
          IF ( Atmosphere(m)%n_Aerosols > 0 ) THEN
            Error_Status_K = CRTM_Compute_AerosolScatter_AD( Atmosphere(m) , &  ! FWD Input
                                                           AtmOptics   , &  ! FWD Input
                                                           AtmOptics_K, &  ! AD  Input
                                                           SensorIndex , &  ! Input
                                                           ChannelIndex, &  ! Input
                                                           Atmosphere_K(ln,m) , &  ! AD  Output
                                                           ASvar         )  ! Internal variable input
            IF ( Error_Status_K /= SUCCESS ) THEN
              WRITE( Message,'("Error computing AerosolScatter_AD for ",a,&
                     &", channel ",i0,", profile #",i0)' ) &
                     TRIM(ChannelInfo(n)%Sensor_ID), ChannelInfo(n)%Sensor_Channel(l), m
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status_K )
              RETURN
            END IF
          END IF
!
          ! Compute the adjoint molecular scattering properties
!!          IF( RTV%Visible_Flag_true ) THEN
!!            Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex)
!!            Error_Status = CRTM_Compute_MoleculeScatter_AD( &
!!                             Wavenumber  , &
!!                             AtmOptics_AD, &
!!                             Atm_AD        )
!!            IF ( Error_Status /= SUCCESS ) THEN
!!              WRITE( Message,'("Error computing MoleculeScatter_AD for ",a,&
!!                     &", channel ",i0,", profile #",i0)' ) &
!!                     TRIM(ChannelInfo(n)%Sensor_ID), &
!!                     ChannelInfo(n)%Sensor_Channel(l), &
!!                     m
!!              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
!!              RETURN
!!            END IF
!!          END IF
!
          ! Compute the adjoint gaseous absorption
          CALL CRTM_Compute_AtmAbsorption_AD( SensorIndex     , &  ! Input
                                              ChannelIndex    , &  ! Input
                                              Predictor       , &  ! FWD Input
                                              AtmOptics_K    , &  ! AD  Input
                                              Predictor_K    , &  ! AD  Output
                                              AAVar             )  ! Internal variable input
        ! Adjoint of the predictor calculations
        CALL CRTM_Compute_Predictors_AD( SensorIndex   , &  ! Input
                                         Atmosphere(m) , &  ! FWD Input
                                         Predictor     , &  ! FWD Input
                                         Predictor_K  , &  ! AD  Input
                                         AncillaryInput, &  ! Input
                                         Atmosphere_K(ln,m), &  ! AD  Output
                                         PVar            )  ! Internal variable input

          ! Postprocess some input data
          ! ...K-matrix of average surface skin temperature for multi-surface types
          CALL CRTM_Compute_SurfaceT_AD( Surface(m), SfcOptics_K, Surface_K(ln,m) )
          
        END DO Channel_Loop

        ! Deallocate local sensor dependent data structures
        ! ...RTV structure
        IF ( RTV_Associated(RTV) ) CALL RTV_Destroy(RTV)
        ! ...Predictor structures
        CALL CRTM_Predictor_Destroy( Predictor )
        CALL CRTM_Predictor_Destroy( Predictor_K )

      END DO Sensor_Loop

      ! Deallocate local sensor independent data structures
      CALL CRTM_Predictor_Destroy( Predictor )
      CALL CRTM_AtmOptics_Destroy( AtmOptics )
      CALL CRTM_SfcOptics_Destroy( SfcOptics )
      ! ...Internal variables
      CALL AOvar_Destroy( AOvar )
      CALL CSvar_Destroy( CSvar )
      CALL ASvar_Destroy( ASvar )
      CALL RTV_Destroy( RTV )
      DEALLOCATE( Height, deltaZ )
    END DO Profile_Loop

  END FUNCTION CRTM_ActiveSensor_K


END MODULE CRTM_ActiveSensor_Module
!
