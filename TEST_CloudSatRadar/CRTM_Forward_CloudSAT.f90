!-------------------------------------------------------------------------------
! PURPOSE: An example code for simulating CloudSat attenuated reflectivity.
!        
! Quanhua Liu, NOAA/NESDIS Center for Satellite Applications and Research
!    quanhua.liu@noaa.gov        September 26, 2022
!-------------------------------------------------------------------------------
PROGRAM CRTM_Forward_CloudSAT
  USE CRTM_Module
  USE CRTM_SpcCoeff
  USE SpcCoeff_Define
  USE CRTM_ActiveSensor_module, ONLY: CRTM_ActiveSensor
  USE read_cloudsat_module, ONLY: load_atm_sfc_geo_2BCLDGPM
  IMPLICIT NONE
  CHARACTER( * ), PARAMETER :: PROGRAM_NAME   = 'CRTM_Forward_CloudSAT'
  INTEGER :: Error_Status, L, i, j, k, m, kk, Allocate_Status, n_Channels
  INTEGER, PARAMETER :: N_SENSORS = 1
  CHARACTER( 256 ) :: Sensor_ID(N_SENSORS)
  CHARACTER( 256 ) :: AerosolCoeff_File,EmisCoeff_File,CloudCoeff_File,fixcrtm
  TYPE( CRTM_Atmosphere_type ), ALLOCATABLE :: Atmosphere(:)
  TYPE( CRTM_Surface_type ), ALLOCATABLE :: Surface(:)
  TYPE(CRTM_Geometry_type), ALLOCATABLE  :: Geometry(:)
      
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE( CRTM_RTSolution_type ), ALLOCATABLE, DIMENSION( :,: ) :: RTSolution
  TYPE( CRTM_Options_type ), ALLOCATABLE  :: Options(:)
  INTEGER :: n_PROFILES, n_Stokes
  INTEGER :: n_Layers, n_Absorbers, n_Clouds, n_Aerosols, chan_index
  CHARACTER( 256 ) :: CRTM_INPUT_FILE
!
  OPEN(11, file='crtm_cloudsat.txt')
  READ(11,'(A)') Sensor_ID(1)
  print *,trim(Sensor_ID(1))
  READ(11,'(A)') CRTM_INPUT_FILE
!
  Error_Status = CRTM_Init( Sensor_Id, &  ! Input... must be an array, hence the (/../)
         ChannelInfo  , &  ! Output
         CloudCoeff_File='Cloud_V3.bin', &
         AerosolCoeff_File='Aerosol_V3.bin', & 
         File_Path='../CRTM_coeff/') 
!
  n_Channels = sum(CRTM_ChannelInfo_n_Channels(ChannelInfo(:))  ) 
  n_Stokes = 1   ! Using scalar RT model

  print *,' number of channels = ',n_Channels,' for ',Sensor_ID(1)
!
!  Read CRTM input data
   CALL load_atm_sfc_geo_2BCLDGPM( atmosphere, surface, Geometry, trim(CRTM_INPUT_FILE) )
!
   n_Profiles  = size( atmosphere )
   print *,' n_Profiles = ',n_Profiles
   !n_Layers    = atmosphere(1)%n_Layers      !Don't use this n_Layers, some profile may have different number of layers than atm(1)
   n_Absorbers = atmosphere(1)%n_Absorbers    ! H2O & O3
   n_Clouds    = atmosphere(1)%n_Clouds   
   n_Aerosols  = 0                     ! for radar
   print *,' n_Clouds = ',n_Clouds, ' n_Absorbers = ',n_Absorbers
!
   Atmosphere(:)%Climatology = US_STANDARD_ATMOSPHERE
   ALLOCATE( RTSolution( n_Channels, n_PROFILES ), STAT = Allocate_Status )  
   print *,' RTSolution allocate ',Allocate_Status
   DO m = 1, n_Profiles
         CALL CRTM_RTSolution_Create( RTSolution(:,m), atmosphere(m)%n_Layers )
   END DO    
  ALLOCATE(Options(n_Profiles), STAT=Allocate_Status ) 
   print *,' Options allocate ',Allocate_Status
  CALL CRTM_Options_Create( Options, n_Channels )
  Options(:)%n_Stokes = n_Stokes
!    Options(1)%RT_Algorithm_Id = RT_ADA !RT_VMOM  !RT_SOI  !RT_ADA 
   DO m = 1, n_Profiles
     Options(m)%Use_Emissivity = .TRUE.
     Options(m)%Use_Direct_Reflectivity = .TRUE.
     Options(m)%Emissivity(:) = (ONE-0.1)
     Options(m)%Direct_Reflectivity(:) = (ONE-Options(m)%Emissivity(:)) * PI  
   END DO
!
!    Geometry(1)%Sensor_Zenith_Angle = 0.0
!

   m = 1368
   Options(m)%Use_n_Streams = .true.
   Options(m)%n_Streams = 8
   Error_Status = CRTM_Forward( Atmosphere(m:m) , &  
                               Surface(m:m)    , & 
                               Geometry(m:m)   , &  
                               ChannelInfo, &  
                               RTSolution(:,m:m),  &
                               Options = Options(m:m) )

  print *,'  channel     Frequency    Radiance    Brightness Temperature '
    DO k = 1, n_Channels
      write(6,'(I5,6x,f12.5,ES14.5,f12.5)') k,SC(1)%Frequency(k),RTSolution(k,1368)%Radiance, &
      RTSolution(k,1368)%Brightness_Temperature
    END DO

!  IF(m > 0) STOP

   Error_Status = CRTM_ActiveSensor( Atmosphere , &  
                               Surface    , & 
                               Geometry   , &  
                               ChannelInfo, &  
                               RTSolution,  &
                               Options = Options )  

 print *,' Error_Status ',Error_Status




 CLOSE(60)
 DO m = 1, n_Profiles
!   IF(maxval(RTSolution(1,m)%Reflectivity_Attenuated) > 0.0 ) THEN
!     print *,m, maxval(RTSolution(1,m)%Reflectivity_Attenuated), &
!       maxval(RTSolution(1,m)%Reflectivity)
     write(60,'(I10,10f10.3)') m, maxval(RTSolution(1,m)%Reflectivity_Attenuated), &
       maxval(RTSolution(1,m)%Reflectivity), &
       (maxval(Atmosphere(m)%Cloud(j)%Water_Content),j=1,4), &
       (maxval(Atmosphere(m)%Cloud(j)%Effective_Radius),j=1,4)
!   END IF
 END DO
 i = 13
 print *,Atmosphere(i)%Cloud(:)%Type
 DO k = 1, Atmosphere(i)%n_layers
   write(6,511) k,Atmosphere(i)%Level_Pressure(k),Atmosphere(i)%Pressure(k), &
   Atmosphere(i)%Temperature(k),Atmosphere(i)%Absorber(k,1), &
   (Atmosphere(i)%Cloud(j)%Water_Content(k),j=1,4)
 END DO
 511 FORMAT(I5,3f9.3,6E12.3)
 CLOSE(60)
END PROGRAM CRTM_Forward_CloudSAT
!   
