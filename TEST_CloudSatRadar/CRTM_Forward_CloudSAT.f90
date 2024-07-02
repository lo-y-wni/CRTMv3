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
  INTEGER :: n_PROFILES, n_Stokes, FOUTbin
  INTEGER :: n_Layers, n_Absorbers, n_Clouds, n_Aerosols, chan_index
  CHARACTER( 256 ) :: CRTM_INPUT_FILE
!
  OPEN(11, file='crtm_cloudsat.txt')
  READ(11,'(A)') Sensor_ID(1)
  print *,trim(Sensor_ID(1))
  READ(11,'(A)') CRTM_INPUT_FILE
!
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


   SC(1)%Is_Active_Sensor  = .TRUE.


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
   print *,' m = ',m,Atmosphere(m:m)%n_Clouds
   Options(m)%Use_n_Streams = .true.
   Options(m)%n_Streams = 8
!   Error_Status = CRTM_Forward( Atmosphere(m:m) , &  
!                               Surface(m:m)    , & 
!                               Geometry(m:m)   , &  
!                               ChannelInfo, &  
!                               RTSolution(:,m:m),  &
!                               Options = Options(m:m) )

  print *,'  channel     Frequency    Radiance    Brightness Temperature '
    DO k = 1, n_Channels
      write(6,'(I5,6x,f12.5,ES14.5,f12.5)') k,SC(1)%Frequency(k),RTSolution(k,1368)%Radiance, &
      RTSolution(k,1368)%Brightness_Temperature
    END DO

!  IF(m > 0) STOP
   print *,'Prf 1 ',atmosphere(1)%n_Clouds,atmosphere(1)%n_layers, &
    atmosphere(1)%level_pressure(atmosphere(1)%n_layers)

   Error_Status = CRTM_ActiveSensor( Atmosphere , &  
                               Surface    , & 
                               Geometry   , &  
                               ChannelInfo, &  
                               RTSolution,  &
                               Options = Options )  

 print *,' Error_Status ',Error_Status

!
 FOUTbin = 61
 CLOSE(FOUTbin)
 OPEN(FOUTbin,file='cpr_ref.bin', FORM='unformatted', STATUS = 'REPLACE', ACCESS='STREAM', &
     CONVERT='LITTLE_ENDIAN')
   print *,' height ',maxval(Atmosphere(1)%Height(:)), minval(Atmosphere(1)%Height(:))
 CALL Dump_Result(FOUTbin, ChannelInfo(1), Atmosphere, RTSolution)     
 CLOSE(FOUTbin)
 OPEN(FOUTbin,file='cpr_ref.txt', FORM='formatted', STATUS = 'REPLACE')
 CALL Print_Result(FOUTbin, ChannelInfo(1), Atmosphere, RTSolution)
 CLOSE(60)
 DO m = 1, n_Profiles
!   IF(maxval(RTSolution(1,m)%Reflectivity_Attenuated) > 0.0 ) THEN
!     print *,m, maxval(RTSolution(1,m)%Reflectivity_Attenuated), &
!       maxval(RTSolution(1,m)%Reflectivity)
     write(60,'(I10,10f12.5)') m, maxval(RTSolution(1,m)%Reflectivity_Attenuated), &
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
 
 CONTAINS

   !----------------------------------------------------------------------------
   SUBROUTINE Print_Result(fid, ChannelInfo, Atm, RTSolution)
   !----------------------------------------------------------------------------
      INTEGER,                     INTENT( IN )  :: fid
      TYPE(CRTM_ChannelInfo_type), INTENT( IN )  :: ChannelInfo
      TYPE(CRTM_Atmosphere_type),  INTENT( IN )  :: Atm(:)
      TYPE(CRTM_RTSolution_type),  INTENT( IN )  :: RTSolution(:,:)

      ! Local
      INTEGER  :: l, k, m, n_Profiles
      REAL(fp) :: data_out(100)
      
      n_Channels = SIZE(RTSolution, DIM=1)
      n_Profiles = SIZE(RTSolution, DIM=2)
      
      DO m = 1, n_Profiles
      DO l = 1, n_Channels                                                                                  
                                                                                                               
            WRITE( fid, '(/7x,"Profile ",i0," Reflectivity for ",a,&
                        &" channel ",i0, ", n_Layers = ", i0)') m, TRIM(ChannelInfo%Sensor_ID),& 
                         ChannelInfo%Sensor_Channel(l), Atm(m)%n_Layers
            
            WRITE( fid,'(/6x, a, i0)')" Pressure(mb)    Reflectivity      Reflectivity_Attenuated"
            
            DO k = 1, Atm(m)%n_Layers                                                                           
              data_out(1:3) = [ Atm(m)%Pressure(k), &
                                RTSolution(l, m)%Reflectivity(k), &                  
                                RTSolution(l, m)%Reflectivity_Attenuated(k)]                                   
              WRITE( fid,'(7x,f8.3,2x,es13.6,2x,es13.6)' ) data_out(1:3)                                        
            END DO                                                                                              
      END DO                                                                                                
      END DO                                                                                                  
   END SUBROUTINE Print_Result

    !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   SUBROUTINE Dump_Result(fid, ChannelInfo, Atm, RTSolution)
   !----------------------------------------------------------------------------
      INTEGER,                     INTENT( IN )  :: fid
      TYPE(CRTM_ChannelInfo_type), INTENT( IN )  :: ChannelInfo
      TYPE(CRTM_Atmosphere_type),  INTENT( IN )  :: Atm(:)
      TYPE(CRTM_RTSolution_type),  INTENT( IN )  :: RTSolution(:,:)
      ! Local
      INTEGER  :: l, k, m, n_Profiles, n_Layers
      REAL(fp) :: data_out(100)
      
      n_Channels = SIZE(RTSolution, DIM=1)
      n_Profiles = SIZE(RTSolution, DIM=2)

    print *,' dump to a binary file ',n_Channels, n_Profiles
    print *,Atm(1)%n_Layers
    DO k = 1, Atm(1)%n_Layers
      write(6,'(I5,3f12.4,4E15.6)') k, Atm(1)%Pressure(k), &
                       RTSolution(1, 1)%Reflectivity(k), &                  
                       RTSolution(1, 1)%Reflectivity_Attenuated(k), &
             Atm(1)%Cloud(:)%Water_Content(k)
    END DO 
    
      write(fid) n_Channels, n_Profiles

      DO l = 1, n_Channels
         write(fid) ChannelInfo%Sensor_Channel(l)
      END DO                                                                                                  

      DO m = 1, n_Profiles                                                                                    
         n_Layers = Atm(m)%n_Layers
         write(fid) n_Layers

         DO l = 1, n_Channels                                                                                  
         DO k = 1, n_Layers
            write(fid) Atm(m)%Pressure(k), &
                       Atm(m)%Height(k), &
                       RTSolution(l, m)%Reflectivity(k), &                  
                       RTSolution(l, m)%Reflectivity_Attenuated(k)                                
         END DO                                                                                                
         END DO                                                                                                
      END DO                                                                                                  

   END SUBROUTINE Dump_Result
 
 
END PROGRAM CRTM_Forward_CloudSAT
!   
