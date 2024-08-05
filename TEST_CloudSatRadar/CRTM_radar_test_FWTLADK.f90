!-------------------------------------------------------------------------------
! PURPOSE: A simple example for testing CRTM operator for active sensors.
!    1. forward simulation
!    2. finite difference for derivative
!    3. Tangent-linear
!    4. comparison between finite derivative and tangent-linear
!    5. Adjoint
!    6. consistency between adjoint and tangent-linear
!    7. K-matrix
!    8. K-matrix (jacobian) vs tangent linear
! 
! Quanhua Liu, NOAA/NESDIS Center for Satellite Applications and Research
!    quanhua.liu@noaa.gov        June 30, 2024
!-------------------------------------------------------------------------------
PROGRAM CRTM_radar_test_FWTLADK
  USE CRTM_Module
  USE CRTM_Model_Profiles   , ONLY: MODEL_LEVEL_PRESSURE, CRTM_Get_Model_Profile 
  USE CRTM_SpcCoeff
  USE SpcCoeff_Define
  USE CRTM_ActiveSensor_module, ONLY: CRTM_ActiveSensor,CRTM_ActiveSensor_TL,  &
           CRTM_ActiveSensor_AD,CRTM_ActiveSensor_K 
  IMPLICIT NONE
  CHARACTER( * ), PARAMETER :: PROGRAM_NAME   = 'CRTM_radar_test_FWTLADK'
  INTEGER :: Error_Status, L, i, j, k, kk, Allocate_Status, n_Channels
  INTEGER, PARAMETER :: N_PROFILES=1, n_Sensors = 1, M_nz = 97  !number of atmospheric levels  
  CHARACTER( 256 ) :: Sensor_ID(n_Sensors)
  CHARACTER( 256 ) :: AerosolCoeff_File,EmisCoeff_File,CloudCoeff_File,fixcrtm
  TYPE( CRTM_Atmosphere_type )   :: Atmosphere(1), Atmosphere_TL(1), Atmosphere_AD(1)
  TYPE( CRTM_Surface_type )      :: Surface(1), Surface_TL(1), Surface_AD(1)
    
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE( CRTM_RTSolution_type ), ALLOCATABLE, DIMENSION( :,: ) :: RTSolution, RTSolution_TL, &
        RTSolution_AD
  TYPE( CRTM_RTSolution_type ), ALLOCATABLE, DIMENSION( :,: ) :: RTSolution1
  TYPE( CRTM_Options_type )  :: Options(N_PROFILES)
  REAL(fp) :: level_p(0:M_nz),level_t(0:M_nz),le_absorber(0:M_nz,2),total_OD,u,tt

  REAL(fp), Allocatable :: transmittance(:)
  REAL(fp) :: deltaX, dy, tTL, tAD
  INTEGER :: mProfile
  INTEGER :: n_Layers, n_Absorbers, n_Clouds, n_Aerosols, chan_index

!  INTEGER, PARAMETER :: absorb_id(3)=(/ 1,3,2/)

  INTEGER :: Err_Atm, Err_Sfc, Ki
  INTEGER :: used_n_channel, n_channel, channel_index(1000), used_channel_index(1000)
!
  TYPE(CRTM_Atmosphere_type)  , DIMENSION(:,:), ALLOCATABLE       :: Atmosphere_K
  TYPE(CRTM_Surface_type)     , DIMENSION(:,:), ALLOCATABLE       :: Surface_K
  TYPE(CRTM_RTSolution_type)  , DIMENSION(:,:), ALLOCATABLE       :: RTSolution_K 
!
!
   Sensor_ID(1)='cpr_cloudsat'
   
   k = 0
   
   IF(k >0) THEN
   Error_Status = CRTM_Init( Sensor_Id, &
                             ChannelInfo, &
                             Aerosol_Model       = 'CRTM', &
                             AerosolCoeff_Format = 'netCDF', &
                             AerosolCoeff_File   = './testinput/AerosolCoeff.nc4', &
!!                             AerosolCoeff_File   = './testinput/Aerosol_V3.bin', &
                             Cloud_Model         = 'CRTM', & 
                             CloudCoeff_Format   = "netCDF", &
!!                             CloudCoeff_File     = './testinput/CloudCoeff.nc', &
                             CloudCoeff_File     = './testinput/Cloud_V3.nc', &
!!                             CloudCoeff_File     = './testinput/CloudCoeff.nc4', &
                             SpcCoeff_Format     = 'netCDF', &
                             TauCoeff_Format     = 'netCDF', &
         File_Path='./testinput/') 
  END IF

   Error_Status = CRTM_Init( Sensor_Id, &
                             ChannelInfo, &
                             Aerosol_Model       = 'CRTM', &
                             AerosolCoeff_Format = 'netCDF', &
                             AerosolCoeff_File   = &
       '../CRTM_Coeff_BigEndian/AerosolCoeff_liu.nc', &
       Cloud_Model         = 'CRTM', & 
       CloudCoeff_Format   = "netCDF", &
       CloudCoeff_File     = '../CRTM_Coeff_BigEndian/CloudCoeff_liu.nc', &
                             SpcCoeff_Format     = 'netCDF', &
                             TauCoeff_Format     = 'netCDF', &
       File_Path='../CRTM_Coeff_BigEndian/') 



   SC(1)%Is_Active_Sensor  = .TRUE.
  n_Channels = sum(CRTM_ChannelInfo_n_Channels(ChannelInfo(:))  ) 
  CALL CRTM_Options_Create( Options, n_Channels )

  print *,' number of channels = ',n_Channels,' for ',Sensor_ID(1)

  IF( trim(Sensor_ID(1)) == 'u.omps-npAllFOV_j2' .or. trim(Sensor_ID(1)) == 'u.omps-tcAllFOV_j2' ) THEN
    PRINT *,' input nFOV, if enter  0  then nadir FOV will be used '
!    READ(5,*)  Options(1)%nFOV
!    IF(  Options(1)%nFOV < 0 )  Options(1)%nFOV = 0
  END IF
!
   !  for a model profile under clear-sky condition.
   n_Layers = 97
   n_Clouds = 1
   n_Aerosols = 0
   n_Absorbers = 2             ! two absorbers, H2O & O3

  ! Allocate atmospheric arrays for atmospheric structure
  Allocate( transmittance(0:n_Layers) )
  CALL CRTM_Atmosphere_Create( Atmosphere,n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atmosphere_TL,n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atmosphere_AD,n_Layers, n_Absorbers, n_Clouds, n_Aerosols )  

  ALLOCATE( Atmosphere_K( n_Channels, n_Profiles ), &
            Surface_K( n_Channels, n_Profiles ), &
            STAT = Allocate_Status )

  IF ( Allocate_Status /= 0 ) THEN                                                       
    Allocate_Status = FAILURE                                                               
    CALL Display_Message( PROGRAM_NAME, &                                                
                          'Error allocating structure arrays Atm_K, Sfc_K & RTSolution_K', &                         
                           Allocate_Status)                                                 
    STOP                                                                                 
  END IF
  ! Allocate atmospheric arrays for atmospheric structure
  CALL CRTM_Atmosphere_Create( Atmosphere,n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atmosphere_TL,n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atmosphere_AD,n_Layers, n_Absorbers, n_Clouds, n_Aerosols )  
  CALL CRTM_Atmosphere_Create( Atmosphere_K,n_Layers, n_Absorbers, n_Clouds, n_Aerosols ) 

   print *,' input a number 1 to 6 for following standard profile '
   print *,' 1 :  TROPICAL '
   print *,' 2 :  MIDLATITUDE_SUMMER '
   print *,' 3 :  MIDLATITUDE_WINTER '
   print *,' 4 :  SUBARCTIC_SUMMER '
   print *,' 5 :  SUBARCTIC_WINTER '
   print *,' 6 :  US_STANDARD_ATMOSPHERE '

   mProfile = 1
  GET_Poption: SELECT CASE (mProfile)

  CASE (TROPICAL)
    Atmosphere(1)%Climatology = TROPICAL      ! using tropical atmosphere

  CASE (MIDLATITUDE_SUMMER)
    Atmosphere(1)%Climatology = MIDLATITUDE_SUMMER

  CASE (MIDLATITUDE_WINTER)
    Atmosphere(1)%Climatology = MIDLATITUDE_WINTER
   
  CASE (SUBARCTIC_SUMMER)
    Atmosphere(1)%Climatology = SUBARCTIC_SUMMER
   
  CASE (SUBARCTIC_WINTER)
    Atmosphere(1)%Climatology = SUBARCTIC_WINTER
   
  CASE (US_STANDARD_ATMOSPHERE)
    Atmosphere(1)%Climatology = US_STANDARD_ATMOSPHERE

  CASE DEFAULT
    Atmosphere(1)%Climatology = US_STANDARD_ATMOSPHERE  
  END SELECT GET_Poption

    Atmosphere(1)%Absorber_Id    = (/ H2O_ID, O3_ID /)
    Atmosphere(1)%Absorber_Units = (/ MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS /)

    Atmosphere(1)%level_pressure(0) = 0.005_fp

   !  -- CALL CRTM ROUTINE TO GET MODEL ATMOSPHERE --
   CALL CRTM_Get_Model_Profile(Atmosphere(1)%Absorber_Id, level_p,level_t,le_Absorber, &
                                Model=Atmosphere(1)%Climatology)

   ! -- SUM OF CHANNELS FROM ALL SENSORS AND ALLOCATE OUTPUT RTSolution --
   ALLOCATE( RTSolution( n_Channels,N_PROFILES ), STAT = Allocate_Status )  
   ALLOCATE( RTSolution1( n_Channels,N_PROFILES ), STAT = Allocate_Status )  
   ! open an array for holding optical depth
   CALL CRTM_RTSolution_Create( RTSolution, Atmosphere(1)%n_Layers )
   CALL CRTM_RTSolution_Create( RTSolution1, Atmosphere(1)%n_Layers )
    
  ALLOCATE( RTSolution_TL( n_Channels,N_PROFILES ), STAT = Allocate_Status ) 
  CALL CRTM_RTSolution_Create( RTSolution_TL, Atmosphere(1)%n_Layers )
  ALLOCATE( RTSolution_AD( n_Channels,N_PROFILES ), STAT = Allocate_Status ) 
  CALL CRTM_RTSolution_Create( RTSolution_AD, Atmosphere(1)%n_Layers )
  ALLOCATE( RTSolution_K( n_Channels,N_PROFILES ), STAT = Allocate_Status ) 
  CALL CRTM_RTSolution_Create( RTSolution_K, Atmosphere(1)%n_Layers )


   !  -- SETTING GEOMETRY PARAMETERS --
   Geometry(1)%Source_Zenith_Angle  = 30.0_fp  !150.00_fp   ! >= 90 at night
   Geometry(1)%Sensor_Zenith_Angle  = ZERO 
   Geometry(1)%Sensor_Scan_Angle  = ZERO   ! 400 km above the height

   ! needed for Cross-scan MW sensors
   Geometry(1)%Sensor_Azimuth_Angle = 0.0_fp   ! needed for WINDSAT sensor
   Geometry(1)%Source_Azimuth_Angle = 150.0_fp ! needed by IR sensors for day time

   !  -- CONSTRUCT LAYER VALUES FROM LEVEL VALUES --
   Atmosphere(1)%Level_Pressure(0) = 0.005_fp
   Layer_Loop: DO kk = 1, Atmosphere(1)%n_Layers
     i = kk 
     Atmosphere(1)%Level_Pressure(kk) = level_p(i)
     Atmosphere(1)%Temperature(kk) = (level_t(i-1)+level_t(i))/2.0_fp
     Atmosphere(1)%Pressure(kk)    = (level_p(i-1)+level_p(i))/2.0_fp
     
     Absorber_Loop: DO k = 1, Atmosphere(1)%n_Absorbers
       IF( k <= 2 ) THEN
         Atmosphere(1)%Absorber(kk,k)=(le_absorber(i-1,k)+le_absorber(i,k))/2.0_fp
       ELSE
         Atmosphere(1)%Absorber(kk,k)= 382.5  !1.75 
       END IF
     END DO Absorber_Loop
   END DO Layer_Loop
  
   !  -- SETTING SURFACE PARAMETERS --
   Surface(1)%Water_Type = 1 
   Surface(1)%Land_Coverage = 0.0_fp
   Surface(1)%Water_Coverage = 1.0_fp
   Surface(1)%Snow_Coverage = 0.0_fp
   Surface(1)%Ice_Coverage = 0.0_fp
   Surface(1)%Water_Temperature = level_t( n_Layers )
   Surface(1)%Land_Temperature = level_t( n_Layers )
   Surface(1)%Wind_Speed = 7.1_fp
   Surface(1)%Wind_Direction = 120.5018_fp
   Surface(1)%Salinity = 33.0_fp


!    Options(1)%RT_Algorithm_Id = RT_ADA !RT_VMOM  !RT_SOI  !RT_ADA 

    Options(1)%Use_Emissivity = .TRUE.
    Options(1)%Use_Direct_Reflectivity = .TRUE.
    Options(1)%Emissivity(:) = (ONE-0.4)
    Options(1)%Direct_Reflectivity(:) = (ONE-Options(1)%Emissivity(:)) !* PI  


    Options(1)%Aircraft_Pressure = Atmosphere(1)%level_pressure(80)

!    Options(1)%obs_4_downward_P = ZERO
    
     IF( Atmosphere(1)%n_Clouds > 0 ) THEN
       Atmosphere(1)%Cloud_Fraction(:) = ONE
       Atmosphere(1)%Cloud(1)%Type = RAIN_CLOUD
       Atmosphere(1)%Cloud(1)%Effective_Radius(60:70) = 500.0_fp
       Atmosphere(1)%Cloud(1)%Effective_Variance(60:70) = 2.0_fp
       Atmosphere(1)%Cloud(1)%Water_Content(60:70) = 0.3
       Atmosphere(1)%Cloud(1)%Water_Content(65) = 0.1
       Atmosphere(1)%Cloud(1)%Water_Content(68) = 0.4
     END IF

     IF( Atmosphere(1)%n_Clouds > 1 ) THEN
       Atmosphere(1)%Cloud(2)%Type = SNOW_CLOUD
       Atmosphere(1)%Cloud(2)%Effective_Radius(60) = 300.0_fp
       Atmosphere(1)%Cloud(2)%Effective_Variance(60) = 2.0_fp
       Atmosphere(1)%Cloud(2)%Water_Content(60) = 0.3
     END IF     

   tt = 28.9644/47.9982 * 0.021415
    total_od = 0.0
    DO K = 1, Atmosphere(1)%n_Layers
    total_od = total_od + Atmosphere(1)%Absorber(k,2)/tt* &
        (Atmosphere(1)%Level_Pressure(k)-Atmosphere(1)%Level_Pressure(k-1))/9.8/10.0
    END DO
! normalize O3 to 300 DU    
    Atmosphere(1)%Absorber(:,2) = Atmosphere(1)%Absorber(:,2)*300.d0/total_od
    
    Geometry(1)%Sensor_Zenith_Angle = 50.0
    print *,' column ozone = ',total_od
    print *,' n_Clouds ',n_Clouds,Geometry(1)%Sensor_Zenith_Angle
   u = cos( Geometry(1)%Sensor_Zenith_Angle*3.14159/180.0)

  ! ---------------------------------------------------------
  ! ***********  nFOV is only used for OMPS   **************

  write(60,'(2I10)') n_Channels, Atmosphere(1)%n_layers
  DO k = 1,  Atmosphere(1)%n_layers
    write(60,'(I5,3f10.3,E15.6)') k, Atmosphere(1)%level_pressure(k),Atmosphere(1)%pressure(k), &
        Atmosphere(1)%Temperature(k),Atmosphere(1)%Absorber(k,1)
  END DO

         CALL CRTM_RTSolution_Zero( RTSolution )
         CALL CRTM_RTSolution_Zero( RTSolution1 )
  print *,' call Forward model '
   Error_Status = CRTM_ActiveSensor( Atmosphere , &  
                               Surface    , & 
                               Geometry   , &  
                               ChannelInfo, &  
                               RTSolution,  &
                               Options = Options )  
  print *,' after Forward model ',Error_Status
   RTSolution1(1,1)%Reflectivity_Attenuated = RTSolution(1,1)%Reflectivity_Attenuated
   RTSolution1(1,1)%Reflectivity = RTSolution(1,1)%Reflectivity
    
   print *,' Forward simulation Reflectivity_Attenuated -9999.0 (default)'
   write(6,'(6ES14.5)') RTSolution1(1,1)%Reflectivity_Attenuated(65)
!  Puturbation
   deltaX = 0.00001_fp
   Atmosphere(1)%Cloud(1)%Water_Content(65) = Atmosphere(1)%Cloud(1)%Water_Content(65) + deltaX
   Error_Status = CRTM_ActiveSensor( Atmosphere , &  
                               Surface    , & 
                               Geometry   , &  
                               ChannelInfo, &  
                               RTSolution,  &
                               Options = Options )  
   print *,' finite difference or derivative '
   dy = (RTSolution(1,1)%Reflectivity_Attenuated(65)-RTSolution1(1,1)%Reflectivity_Attenuated(65))/deltaX
   print *,' dy = ',dy
   Atmosphere(1)%Cloud(1)%Water_Content(65) = Atmosphere(1)%Cloud(1)%Water_Content(65) - deltaX
!  
   CALL CRTM_Atmosphere_Zero( Atmosphere_TL )   
   CALL CRTM_Surface_Zero( Surface_TL )   
   Atmosphere_TL(1)%Cloud(1)%Water_Content(65) = ONE
   
   CALL CRTM_RTSolution_Zero( RTSolution )
   CALL CRTM_RTSolution_Zero( RTSolution_TL )
 
   print *,' TL '  
   Error_Status = CRTM_ActiveSensor_TL( Atmosphere , &  
                               Surface    , & 
                               Atmosphere_TL, &  ! TL  Input, M
                               Surface_TL   , &  ! TL  Input, M                               
                               Geometry   , &  
                               ChannelInfo, &  
                               RTSolution,  &
                               RTSolution_TL, &  ! Output, L x M
                               Options = Options ) 
   print *,' tangent-linear vs finite difference TL, dy, (TL-dy)/TL'
   print *, 'dy ',dy, RTSolution_TL(1,1)%Reflectivity_Attenuated(65), &
             (RTSolution_TL(1,1)%Reflectivity_Attenuated(65)-dy)/RTSolution_TL(1,1)%Reflectivity_Attenuated(65)
   print *,' Forward simulation by calling CRTM_ActiveSensor_TL'
   write(6,'(6ES14.5)') RTSolution1(1,1)%Reflectivity_Attenuated(65)
!  write(6,'(6E14.5)') RTSolution_TL(1,1)%Reflectivity_Attenuated
        
   print *,' Adjoint ' 
   CALL CRTM_RTSolution_Zero( RTSolution )
   CALL CRTM_RTSolution_Zero( RTSolution_AD )
   CALL CRTM_Atmosphere_Zero( Atmosphere_AD )   
   CALL CRTM_Surface_Zero( Surface_AD )   
   RTSolution_AD(1,1)%Reflectivity_Attenuated(:) = ZERO
   RTSolution_AD(1,1)%Reflectivity(:) = ZERO
!   RTSolution_AD(1,1)%Reflectivity_Attenuated(:) = RTSolution_TL(1,1)%Reflectivity_Attenuated(:)
   RTSolution_AD(1,1)%Reflectivity_Attenuated(65) = RTSolution_TL(1,1)%Reflectivity_Attenuated(65)
   Error_Status = CRTM_ActiveSensor_AD( Atmosphere , &  
                               Surface    , & 
                               RTSolution_AD, &  ! Output, L x M
                               Geometry   , &  
                               ChannelInfo, &  
                               Atmosphere_AD, &  ! TL  Input, M
                               Surface_AD   , &  ! TL  Input, M    
                               RTSolution,  &
                               Options = Options ) 
!        tAD = sum(Atmosphere_AD(1)%Cloud(1)%Water_Content(:) * Atmosphere_TL(1)%Cloud(1)%Water_Content(:))
!        tTL = sum(RTSolution_TL(1,1)%Reflectivity_Attenuated(:)*RTSolution_AD(1,1)%Reflectivity_Attenuated(:))

        tAD = (Atmosphere_AD(1)%Cloud(1)%Water_Content(65) * Atmosphere_TL(1)%Cloud(1)%Water_Content(65))
        tTL = (RTSolution_TL(1,1)%Reflectivity_Attenuated(65)*RTSolution_AD(1,1)%Reflectivity_Attenuated(65))        
        print *,' TL, AD  (TL-AD)/TL*100.0 (5)',tTL, tAD, (tTL-tAD)/tAD * 100.0_fp
   print *,' Forward simulation by calling CRTM_ActiveSensor_AD'
   write(6,'(6ES14.5)') RTSolution1(1,1)%Reflectivity_Attenuated(65)
        print *,' K-matrix ' 
   CALL CRTM_RTSolution_Zero( RTSolution )
   CALL CRTM_RTSolution_Zero( RTSolution_K )
   CALL CRTM_Atmosphere_Zero( Atmosphere_K )   
   CALL CRTM_Surface_Zero( Surface_K )   
   
   RTSolution_K(1,1)%Reflectivity_Attenuated(65) = ONE
   Error_Status = CRTM_ActiveSensor_K( Atmosphere , &  
                               Surface    , & 
                               RTSolution_K, &  ! Output, L x M
                               Geometry   , &  
                               ChannelInfo, &  
                               Atmosphere_K, &  ! TL  Input, M
                               Surface_K   , &  ! TL  Input, M    
                               RTSolution,  &
                               Options = Options ) 
!   write(6,'(6ES13.4)') Atmosphere_K(1,1)%Cloud(1)%Water_Content
   dy = RTSolution_TL(1,1)%Reflectivity_Attenuated(65)
    print *,' TL vs K  dy, K, (dy-K)/dy*100.0 ',dy, Atmosphere_K(1,1)%Cloud(1)%Water_Content(65), &
    (dy-Atmosphere_K(1,1)%Cloud(1)%Water_Content(65))/dy*100.0
   print *,' Forward simulation by calling CRTM_ActiveSensor_K'
   write(6,'(6ES14.5)') RTSolution1(1,1)%Reflectivity_Attenuated(65)
   Error_Status = CRTM_Destroy( ChannelInfo )
END PROGRAM CRTM_radar_test_FWTLADK
!   
