

PROGRAM mainProgram

   USE CRTM_Module
   USE CRTM_SpcCoeff, ONLY: SC
   IMPLICIT NONE
      
   CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'crtm_Z'
   CHARACTER(*), PARAMETER :: PROGRAM_RCS_ID = ''

   INTEGER, PARAMETER  :: nStrLen=256,mxStrLen=5000
   INTEGER, PARAMETER  :: FOUTtxt=121, FOUTbin=131
   

   !--- DEFINE THE CRTM INTERFACE STRUCTURES
   TYPE(CRTM_ChannelInfo_type) ,ALLOCATABLE :: ChannelInfo(:)      ![nSensors]
   TYPE(CRTM_Atmosphere_type)  ,ALLOCATABLE :: Atm(:)              ![nProfiles]
   TYPE(CRTM_Surface_type)     ,ALLOCATABLE :: Sfc(:)              ![nProfiles]
   TYPE(CRTM_Geometry_type)    ,ALLOCATABLE :: Geo(:)              ![nProfiles]
   TYPE(CRTM_RTSolution_type)  ,ALLOCATABLE :: RTSolution(:,:)     ![nChannels,nProfiles]
         
      
   !--- Other local variables
   CHARACTER(nStrLen) :: Message
   CHARACTER(64)      :: Sensor_Ids(10) 
   INTEGER            :: Error_Status, fu
   INTEGER            :: Allocate_Status, Allocate_Status_arr(3)
   INTEGER            :: i, k, l, m, n, s, ipf                      
   INTEGER            :: n_Channels
   INTEGER            :: n_Profiles
   INTEGER            :: n_Layers, n_Levels
   INTEGER            :: n_Absorbers 
   INTEGER            :: n_Clouds    
   INTEGER            :: n_Aerosols  
   INTEGER            :: n_Sensors
   
   !--- Namelist "CRTM_Config_Type"
   TYPE CRTM_Config_Type  
      character(nStrLen) :: SensorID           !'hirs4_n18, amsua_n18'. The last work must be a valid sensor name, not a ",".
      logical            :: RADAR              !if True, the sensor is an Radar instrument
      character(nStrLen) :: input_FileName
      character(nStrLen) :: AtmosProfile_ID    !'CLOUDSAT','2BCLDGPM','UMBC48'
      character(nStrLen) :: Surface_ID         !'MIXED_SURFACES','OCEAN','LAND'
      character(nStrLen) :: Model_Type         !'FORWARD','TANGENT_LINEAR','ADJOINT','K_MATRIX'
      CHARACTER(nStrLEN) :: Coeff_dir          !CRTM-based Coeff files path 
      real               :: Sensor_zenith_angle
      integer            :: numProfToPrint     !number of profiles to write to text format output file
   ENDTYPE CRTM_Config_Type
   TYPE(CRTM_Config_Type) :: config
   NAMELIST/config_nml/config

  
   
   !--- Program header
   CALL Program_Message( PROGRAM_NAME, 'Program to compute radar relectivity', PROGRAM_RCS_ID )
   
   open( file='cfg.nml', newunit=fu, action='read')
   READ( fu, nml=config_nml )
     
   
   !--- Parse sensor IDs and put them into an array   
   n = LEN( TRIM( config%SensorID ))
   l = 1
   n_Sensors = 0
   DO i = 1, n
     IF( config%SensorID(i:i) == ",")THEN
       n_Sensors = n_Sensors + 1
       Sensor_Ids(n_Sensors) = ADJUSTL( config%SensorID(l:i-1) )
       l = i+1
     END IF
   END DO 
   n_Sensors = n_Sensors + 1
   Sensor_Ids( n_Sensors) = ADJUSTL( config%SensorID(l:n))
   
   
   
   ! --- INITIALIZE CRTM
   ! If the optional File_Path argument does not present, 
   ! the current directory is assumed for coefficient files. 
   !
   ALLOCATE( ChannelInfo( n_Sensors),  STAT=Allocate_Status )

   WRITE( *,'(/5x,"Initializing the CRTM...")' )
   Error_Status = CRTM_Init( Sensor_Ids(1:n_Sensors), &
                             ChannelInfo, &
                             Aerosol_Model       = 'CRTM', &
                             AerosolCoeff_Format = 'netCDF', &
                             AerosolCoeff_File   = trim(config%Coeff_dir)//"/AerosolCoeff.nc4", &
                             Cloud_Model         = 'CRTM', & 
                             CloudCoeff_Format   = "netCDF", &
                             CloudCoeff_File     = trim(config%Coeff_dir)//"/CloudCoeff.nc4", &
                             SpcCoeff_Format     = 'netCDF', &
                             TauCoeff_Format     = 'netCDF', &
                             File_Path           = TRIM( config%Coeff_dir)  )
   
   ! We input one sensor only.
   if (config%RADAR) SC(1)%Is_Active_Sensor  = .TRUE.



   !--- (1) Allocate memory for the interface strucutres Atm and Sfc (contain state variables) 
   !    (2) Assign data for their member variables, given a user selected case 
   CALL input_atm_sfc_geo( config%input_Filename, atm, sfc, geo, &
                           AtmosProfile_ID = config%AtmosProfile_ID, &
                           Surface_ID      = config%Surface_ID )


   !--- Get dimensions. Assuming all profiles have the same dimensions.
   ! Need a check: if (any( n_Layers(:) /= n_Layers(1) ) ) STOP
   !
   n_Profiles  = size( atm )
   !n_Layers    = atm(1)%n_Layers 
   n_Absorbers = atm(1)%n_Absorbers 
   n_Clouds    = atm(1)%n_Clouds   
   n_Aerosols  = atm(1)%n_Aerosols 


   !--- Sensor loop start 
   Sensor_Loop: DO s = 1, n_Sensors
   
      !--- print a message on screen
      WRITE( *,'(/5x,"*** Sensor ID: ",a, "; Atmos. Profiles: ",a, "; Surface Type: ",a,  "; Model Type: ",a, " ***")' )&
               TRIM( Sensor_Ids(s) ), &
               TRIM( config%AtmosProfile_ID ), &
               TRIM( config%Surface_ID ), &
               TRIM( config%Model_Type )
      
      
      !--- Allocate memory for the RTSolution arrays (needed for all models)    
      n_Channels = ChannelInfo(s)%n_Channels       
      ALLOCATE( RTSolution( n_Channels, n_Profiles ), STAT=Allocate_Status )

      do m=1,n_Profiles
         CALL CRTM_RTSolution_Create( RTSolution(:,m), atm(m)%n_Layers )
      enddo                  
      
               
      !--- Open output file
      OPEN(FOUTtxt, FILE = TRIM( Sensor_Ids(s) )//"."// &
                           TRIM( config%AtmosProfile_ID )//"."// &
                           TRIM( config%Surface_ID )//"."//&
                           TRIM( config%Model_Type )//".result.txt", STATUS = 'REPLACE')
      OPEN(FOUTbin, FILE = TRIM( Sensor_Ids(s) )//"."// &
                           TRIM( config%AtmosProfile_ID )//"."// &
                           TRIM( config%Surface_ID )//"."//&
                           TRIM( config%Model_Type )//".result.dat", STATUS = 'REPLACE', ACCESS="STREAM")
                        

      SELECT CASE (config%Model_Type)
      
         CASE( 'FORWARD' )

            Error_Status = CRTM_Forward( Atm              , &                               
                                         Sfc              , &                               
                                         Geo              , &                               
                                         ChannelInfo(s:s) , &                               
                                         RTSolution )                                       
            IF ( Error_Status /= SUCCESS ) THEN                                               
               CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Forward Model', Error_Status)                                             
               STOP                                                                            
            END IF                                                                            
            
            !--- Print out FW results 
            CALL Print_Result(FOUTtxt, ChannelInfo(s), Atm, RTSolution, config%numProfToPrint)                            
            CALL Dump_Result(FOUTbin, ChannelInfo(s), Atm, RTSolution)                
                                                                          
         CASE DEFAULT                                                                        
            Error_Status = FAILURE                                                               
            WRITE( Message,'("Error can not find a match for Model ID #: ",i0)' ) config%Model_Type   
            CALL Display_Message( PROGRAM_NAME, TRIM(Message), Error_Status)                                                    
            STOP                                                                              
         
      END SELECT !SELECT CASE (Model_Type)
            
      CLOSE( FOUTtxt )
      CLOSE( FOUTbin )
      
      
      !--- Deallocate memory for sensor dependent interface structures
      !
      Allocate_Status_arr  = 0
      CALL CRTM_RTSolution_Destroy( RTSolution )
      DEALLOCATE(RTSolution, Geo, STAT = Allocate_Status_arr(2))                                     
      
   
   END DO Sensor_Loop
   
   
   !--- Deallocate memory for the sensor/model independent interface structures   
   CALL CRTM_Atmosphere_Destroy( Atm )
   DEALLOCATE(Atm, Sfc, STAT = Allocate_Status)
      
   !  **** DESTROY THE CRTM ****
   WRITE( *, '( /5x, "Destroying the CRTM..." )' )
   Error_Status = CRTM_Destroy( ChannelInfo )

   
CONTAINS !=================== Internal subroutines =============================

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   subroutine input_atm_sfc_geo( input_Filename, atm, sfc, geo, &
                                 AtmosProfile_ID, Surface_ID )
   !----------------------------------------------------------------------------
      IMPLICIT NONE

      character(*)                             ,intent(in)    :: input_Filename
      TYPE(CRTM_Atmosphere_type)  ,ALLOCATABLE ,intent(out)   :: atm(:)       ![nProfiles]  
      TYPE(CRTM_Surface_type)     ,ALLOCATABLE ,intent(out)   :: sfc(:)       ![nProfiles]
      TYPE(CRTM_Geometry_type)    ,ALLOCATABLE ,intent(out)   :: geo(:)       ![nProfiles]
      character(*)                ,optional    ,intent(in)    :: AtmosProfile_ID
      character(*)                ,optional    ,intent(in)    :: Surface_ID
        
      integer        :: i,j,k,m
      CHARACTER(256) :: Message
      INTEGER        :: n_Profiles, n_Layers
      INTEGER        :: Error_Status, Allocate_Status
      character(64)  :: ProfID
      character(64)  :: SfcID



      ProfID = '2BCLDGPM'
      if (present( AtmosProfile_ID )) ProfID = AtmosProfile_ID
      
      SfcID = ''
      if (present( Surface_ID )) SfcID = Surface_ID


      SELECT CASE( TRIM( ProfID ) )
      case( "2BCLDGPM", "CLOUDSAT" ) 

         CALL load_atm_sfc_geo_2BCLDGPM( atm, sfc, geo, input_Filename )
   
      CASE DEFAULT
         Error_Status = FAILURE 
         WRITE( Message,'("Error can not find a match for the input case #: ",A)' ) AtmosProfile_ID
         CALL Display_Message( PROGRAM_NAME, TRIM(Message), Error_Status) 
         STOP
      END SELECT !SELECT CASE( Case_ID )
   
   end subroutine !subroutine input_atm_sfc_geo



   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   subroutine load_atm_sfc_geo_2BCLDGPM( atm, sfc, geo, input_Filename )
   !----------------------------------------------------------------------------
      USE CRTM_Module    ,ONLY: CRTM_Atmosphere_type, &
                                CRTM_Surface_type, &
                                CRTM_Geometry_type, &
                                CRTM_Atmosphere_Create, &
                                CRTM_Atmosphere_Associated, &
                                CRTM_Surface_Create, &
                                CRTM_Geometry_Create, &
                                Display_Message, &
                                H2O_ID, CO2_ID, O3_ID, N2O_ID, CO_ID, CH4_ID, O2_ID, &
                                VOLUME_MIXING_RATIO_UNITS, &
                                  MASS_MIXING_RATIO_UNITS, &
                                INVALID_CLOUD,  &
                                WATER_CLOUD, &
                                ICE_CLOUD, &
                                RAIN_CLOUD, &
                                SNOW_CLOUD, &
                                GRAUPEL_CLOUD, &
                                HAIL_CLOUD                                
      IMPLICIT NONE

      TYPE(CRTM_Atmosphere_type)  ,ALLOCATABLE ,intent(out)   :: atm(:)       ![nProfiles]  
      TYPE(CRTM_Surface_type)     ,ALLOCATABLE ,intent(out)   :: sfc(:)       ![nProfiles]
      TYPE(CRTM_Geometry_type)    ,ALLOCATABLE ,intent(out)   :: geo(:)       ![nProfiles]
      character(*)                             ,intent(in)    :: input_Filename

      INTEGER ,PARAMETER :: OC_TYP       = 0   !Ocean surface type  
      INTEGER ,PARAMETER :: SEAICE_TYP   = 1   !Sea-Ice surface type
      INTEGER ,PARAMETER :: LD_TYP       = 2   !Land surface type   
      INTEGER ,PARAMETER :: SNOW_TYP     = 3   !Snow surface type   
      INTEGER ,PARAMETER :: Desert_TYP   = 4   !Desert surface type 
      INTEGER ,PARAMETER :: COAST_TYP    = 6   !Coast surface type  
      
      integer              :: i,j,k,m, ic,ia
      INTEGER              :: Error_Stus, Allocate_Status
      integer              :: inFileLun
      integer(4)           :: nProf, nLay, nAbsorb, nCloud, nAerosol, nSensor
      integer(4)           :: cldType, iTypSfc
      real(8)              :: WindSpeed, WinDir
      real(8)              :: Tskin, lat, lon, altitude, satZenAng, solZenAng=100      
      real(8), allocatable :: Plev(:), Play(:), Tlay(:), GasLay(:,:), cldN(:), Reff(:), WC(:)



      !--- Allocate AtmProfile 
      !
      inFileLun = 101
      OPEN( inFileLun, file=trim(input_FileName), STATUS="OLD", ACCESS="STREAM", CONVERT="big_endian" )
      read( inFileLun ) nProf             ,&
                        nLay              ,&
                        nAbsorb           ,&
                        nCloud            ,&
                        nAerosol          ,&
                        nSensor

      !--- Allocate the Atm, sfc and geom structure arrays
      ALLOCATE( atm( nProf ), STAT = Allocate_Status )
      ALLOCATE( sfc( nProf ), STAT = Allocate_Status )            
      ALLOCATE( geo( nProf ), STAT = Allocate_Status )            


      !--- Loop over profiles
      Profile_Loop: DO m = 1, nProf

         call CRTM_Atmosphere_Create( Atm(m), &
                                      nLay, &
                                      nAbsorb, &
                                      nCloud, &
                                      nAerosol )
         
         !#--- allocate temporary valueable for reading float64 input data.
         allocate( Plev(0:nLay), Play(nLay), Tlay(nLay), GasLay(nLay,nAbsorb), &
                   cldN(nLay), Reff(nLay), WC(nLay))

         !--- Copy over the data
         ! In Atm structrue to be input to CRTM, the order should be from TOA -> surface
         !
         read( inFileLun ) Plev(0:nLay), Play(1:nLay), Tlay(1:nLay), GasLay(1:nLay,1:nAbsorb)
         !Atm(m)%Climatology    = Climatology_Model    ,& 
         Atm(m)%Level_Pressure( 0:nLay )            = Plev(  0:nLay)
         Atm(m)%Pressure(       1:nLay )            = Play(  1:nLay)
         Atm(m)%Temperature(    1:nLay )            = Tlay(  1:nLay)
         Atm(m)%Absorber(       1:nLay, 1:nAbsorb ) = GasLay(1:nLay,1:nAbsorb) ![H2O(mixing ratio, g/kg),O3(ppmv),CO,CO2,CH4,N2O]
         
         Atm(m)%Absorber_ID(1) = H2O_ID
         Atm(m)%Absorber_ID(2) =  O3_ID
         Atm(m)%Absorber_Units(1) =   MASS_MIXING_RATIO_UNITS   !H2O:MASS_MIXING_RATIO_UNITS(g/kg), O3:VOLUME_MIXING_RATIO_UNITS(ppmv)
         Atm(m)%Absorber_Units(2) = VOLUME_MIXING_RATIO_UNITS   !H2O:MASS_MIXING_RATIO_UNITS(g/kg), O3:VOLUME_MIXING_RATIO_UNITS(ppmv)

         read( inFileLun ) cldN(1:nLay)
         Atm(m)%Cloud_Fraction( 1:nLay ) = cldN(1:nLay)

         do ic = 1,nCloud
            read( inFileLun ) cldType, Reff(1:nLay), WC(1:nLay)
            Atm(m)%Cloud(ic)%Type = cldType
            Atm(m)%Cloud(ic)%Water_Content( 1:nLay )    = WC(1:nLay)   !kg/m^2
            Atm(m)%Cloud(ic)%Effective_Radius( 1:nLay ) = Reff(1:nLay) !micron
            !Atm(m)%Cloud(ic)%Effective_Variance( 1:nLay )=  ,&  !micron^2 (not used)
         enddo
         

         !---Fill Surface and Options Structure
         !
         read( inFileLun ) iTypSfc, Tskin  !,Wind_Speed, WinDir

         IF ( iTypSfc .eq. OC_TYP) THEN
            sfc(m)%Water_Coverage      = 1
            sfc(m)%Ice_Coverage        = 0
            sfc(m)%Land_Coverage       = 0
            sfc(m)%Snow_Coverage       = 0
            sfc(m)%Water_Temperature   = Tskin
            sfc(m)%Salinity            = 33.0
            Sfc%Wind_Speed             = 7.1 !Scene%WindDir
            Sfc%Wind_Direction         = 120.0
         ENDIF
         IF ( iTypSfc .eq. SEAICE_TYP) THEN
            sfc(m)%Water_Coverage      = 0
            sfc(m)%Ice_Coverage        = 1
            sfc(m)%Land_Coverage       = 0
            sfc(m)%Snow_Coverage       = 0
            sfc(m)%Ice_Temperature     = Tskin
         ENDIF
         IF ( iTypSfc .eq. LD_TYP) THEN
            sfc(m)%Water_Coverage      = 0
            sfc(m)%Ice_Coverage        = 0
            sfc(m)%Land_Coverage       = 1
            sfc(m)%Snow_Coverage       = 0
            sfc(m)%Land_Temperature    = Tskin
            sfc(m)%Soil_Temperature    = Tskin
         ENDIF
         IF ( iTypSfc .eq. SNOW_TYP) THEN
            sfc(m)%Water_Coverage      = 0
            sfc(m)%Ice_Coverage        = 0
            sfc(m)%Land_Coverage       = 0
            sfc(m)%Snow_Coverage       = 1
            sfc(m)%Snow_Temperature    = Tskin
         ENDIF
         IF ( iTypSfc .eq. COAST_TYP) THEN
            !---Half Ocean/Half Land
            sfc(m)%Water_Coverage      = 0.5
            sfc(m)%Ice_Coverage        = 0
            sfc(m)%Land_Coverage       = 0.5
            sfc(m)%Snow_Coverage       = 0
            sfc(m)%Water_Temperature   = Tskin
            sfc(m)%Land_Temperature    = Tskin
            sfc(m)%Soil_Temperature    = Tskin
         ENDIF



         !---Fill GeometryInfo Structure
         !
         call CRTM_Geometry_Create( geo(m) )
         
         read( inFileLun ) lat, lon, altitude, satZenAng!, solZenAng

         IF (lon .lt. 0) lon = abs(lon)+180
         IF (ABS(solZenAng) >= 180.) solZenAng = 100.  ! turn it off
         
         geo(m)%Latitude             = lat
         geo(m)%Longitude            = lon
         geo(m)%Surface_Altitude     = altitude
         geo(m)%Sensor_Zenith_Angle  = satZenAng

         deallocate( Plev, Play, Tlay, GasLay, cldN, Reff, WC)

      END DO Profile_Loop

   END subroutine


   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------
   SUBROUTINE Print_Result(fid, ChannelInfo, Atm, RTSolution, numProfToPrint)
   !----------------------------------------------------------------------------
      INTEGER,                     INTENT( IN )  :: fid
      TYPE(CRTM_ChannelInfo_type), INTENT( IN )  :: ChannelInfo
      TYPE(CRTM_Atmosphere_type),  INTENT( IN )  :: Atm(:)
      TYPE(CRTM_RTSolution_type),  INTENT( IN )  :: RTSolution(:,:)
      integer,                     intent( in )  :: numProfToPrint

      ! Local
      INTEGER  :: l, k, m, n_Profiles
      REAL(fp) :: data_out(100)
      
      n_Channels = SIZE(RTSolution, DIM=1)
      n_Profiles = SIZE(RTSolution, DIM=2)
      
      DO m = 1, min( n_Profiles, numProfToPrint)   
      DO l = 1, n_Channels                                                                                  
                                                                                                               
            WRITE( fid, '(/7x,"Profile ",i0," Reflectivity for ",a,&
                        &" channel ",i0, ", n_Layers = ", i0)') m, TRIM(ChannelInfo%Sensor_ID),& 
                         ChannelInfo%Sensor_Channel(l), Atm(m)%n_Layers
            
            WRITE( fid,'(/6x, a, i0)')" Pressure(mb)  Height(km)  Reflectivity  Reflectivity_Attenuated"
            
            DO k = 1, Atm(m)%n_Layers                                                                           
              data_out(1:4) = [ Atm(m)%Pressure(k), &
                                Atm(m)%Height(k), &
                                RTSolution(l, m)%Reflectivity(k), &                  
                                RTSolution(l, m)%Reflectivity_Attenuated(k)]                                   
              WRITE( fid,'(7x,f8.3,2x,f8.2,2x,es13.6,2x,es13.6)' ) data_out(1:4)                                        
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
      
      n_Channels = SIZE(RTSolution, DIM=1)
      n_Profiles = SIZE(RTSolution, DIM=2)
      
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


END PROGRAM
