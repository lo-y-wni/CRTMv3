!
Module read_cloudsat_module

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds              , ONLY: fp
  USE Message_Handler         , ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters         , ONLY: ZERO, ONE, TWO, FOUR, FIVE, TEN, PI, MISSING_REFL

  USE CRTM_Atmosphere_Define  , ONLY: CRTM_Atmosphere_type
  USE CRTM_Surface_Define     , ONLY: CRTM_Surface_type
  USE CRTM_Geometry_Define    , ONLY: CRTM_Geometry_type, CRTM_Geometry_IsValid  
  
  ! Disable all implicit typing
  IMPLICIT NONE
  ! --------------------
  ! Default visibilities
  ! --------------------
  ! Everything private by default
  PRIVATE
  ! RTSolution structure entities
  ! ...Datatypes
  PUBLIC :: load_atm_sfc_geo_2BCLDGPM

  ! -----------------
  ! Module parameters
  ! -----------------

CONTAINS

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
      CHARACTER(256)       :: Message
      INTEGER              :: Error_Status, Allocate_Status
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

         !--- Allocate the current profile Atmosphere structure
         !
         ! TYPE :: CRTM_Atmosphere_type
         !    LOGICAL :: Is_Allocated = .FALSE.
         ! 
         !    INTEGER :: Max_Layers   = 0  ! K dimension
         !    INTEGER :: n_Layers     = 0  ! Kuse dimension
         !    INTEGER :: n_Absorbers  = 0  ! J dimension
         !    INTEGER :: Max_Clouds   = 0  ! Nc dimension
         !    INTEGER :: n_Clouds     = 0  ! NcUse dimension
         !    INTEGER :: Max_Aerosols = 0  ! Na dimension
         !    INTEGER :: n_Aerosols   = 0  ! NaUse dimension
         !    
         !    INTEGER :: n_Added_Layers = 0
         !    
         !    INTEGER :: Climatology = US_STANDARD_ATMOSPHERE !! Climatology model associated with the profile
         ! 
         !    INTEGER, ALLOCATABLE :: Absorber_ID(:)    ! J
         !    INTEGER, ALLOCATABLE :: Absorber_Units(:) ! J
         ! 
         !    REAL(fp), ALLOCATABLE :: Level_Pressure(:)  ! 0:K
         !    REAL(fp), ALLOCATABLE :: Pressure(:)        ! K
         !    REAL(fp), ALLOCATABLE :: Temperature(:)     ! K
         !    REAL(fp), ALLOCATABLE :: Absorber(:,:)      ! K x J
         !    
         !    TYPE(CRTM_Cloud_type),   ALLOCATABLE :: Cloud(:)    ! Nc !! Clouds associated with each profile
         !    
         !    TYPE(CRTM_Aerosol_type), ALLOCATABLE :: Aerosol(:)  ! Na !! Aerosols associated with each profile
         ! END TYPE CRTM_Atmosphere_type
         !
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

            !if ( Atm(m)%Cloud(ic)%Type == WATER_CLOUD ) then
            !   Atm(m)%Cloud(1)%Type                  = WATER_CLOUD
            !   Atm(m)%Cloud(1)%Effective_Radius(:)   = 30          !micron
            !   Atm(m)%Cloud(1)%Effective_Variance(:) = 0.0         !micron^2 (not used)
            !   Atm(m)%Cloud(1)%Water_Content(:)      = lwc(:)      !kg/m^2
            !endif            
            !if ( Atm(m)%Cloud(ic)%Type == ICE_CLOUD ) then
            !   Atm(m)%Cloud(5)%Type                  = ICE_CLOUD
            !   Atm(m)%Cloud(5)%Effective_Radius(:)   = 30.0         !micron
            !   Atm(m)%Cloud(5)%Effective_Variance(:) = 0.0          !micron^2 (not used)
            !   Atm(m)%Cloud(5)%Water_Content(:)      = iwc(:)       !kg/m^2
            !endif
            !if ( Atm(m)%Cloud(ic)%Type == RAIN_CLOUD ) then
            !   Atm(m)%Cloud(2)%Type                  = RAIN_CLOUD
            !   Atm(m)%Cloud(2)%Effective_Radius(:)   = 500          !micron
            !   Atm(m)%Cloud(2)%Effective_Variance(:) = 0.0          !micron^2 (not used)
            !   Atm(m)%Cloud(2)%Water_Content(:)      = Rain(:)      !kg/m^2
            !endif
            !if ( Atm(m)%Cloud(ic)%Type == SNOW_CLOUD ) then
            !   Atm(m)%Cloud(4)%Type                  = SNOW_CLOUD
            !   Atm(m)%Cloud(4)%Effective_Radius(:)   = 300          !micron
            !   Atm(m)%Cloud(4)%Effective_Variance(:) = 0.0          !micron^2 (not used)
            !   Atm(m)%Cloud(4)%Water_Content(:)      = Snow(:)      !kg/m^2
            !endif
            !if ( Atm(m)%Cloud(ic)%Type == GRAUPEL_CLOUD ) then
            !   Atm(m)%Cloud(3)%Type                  = GRAUPEL_CLOUD
            !   Atm(m)%Cloud(3)%Effective_Radius(:)   = 1000         !micron
            !   Atm(m)%Cloud(3)%Effective_Variance(:) = 0.0          !micron^2 (not used)
            !   Atm(m)%Cloud(3)%Water_Content(:)      = Graupel(:)   !kg/m^2
            !endif
         enddo
         
         !do ia = 1,nAerosol
         !enddo 

         !---Fill Surface and Options Structure
         !   Don't call "CRTM_Surface_Create" if not intent to retrieve the surface emissivity. Once
         !   the sfc structure allocated array for "sfc%Sensor_data", the surface structure need to 
         !   have the all "sfc%Sensor_data" components to be assgined.
         !   
         !call CRTM_Surface_Create( sfc(m), n_Channels=1 )

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
            !Sfc%Ice_Type            = DEFAULT_ICE_TYPE
            !Sfc%Ice_Thickness       = DEFAULT_ICE_THICKNESS
            !Sfc%Ice_Density         = DEFAULT_ICE_DENSITY
            !Sfc%Ice_Roughness       = DEFAULT_ICE_ROUGHNESS
         ENDIF
         IF ( iTypSfc .eq. LD_TYP) THEN
            sfc(m)%Water_Coverage      = 0
            sfc(m)%Ice_Coverage        = 0
            sfc(m)%Land_Coverage       = 1
            sfc(m)%Snow_Coverage       = 0
            sfc(m)%Land_Temperature    = Tskin
            sfc(m)%Soil_Temperature    = Tskin
            !---Randomize or fix surface properties
            !rn=1
            !call random_number(rn)
            !Sfc%Soil_Moisture_Content = 0.2*rn
            !call random_number(rn)
            !Sfc%Vegetation_Fraction   = 0.5*rn
            !call random_number(rn)
            !Sfc%Canopy_Water_Content  = 0.2*rn
         ENDIF
         IF ( iTypSfc .eq. SNOW_TYP) THEN
            sfc(m)%Water_Coverage      = 0
            sfc(m)%Ice_Coverage        = 0
            sfc(m)%Land_Coverage       = 0
            sfc(m)%Snow_Coverage       = 1
            sfc(m)%Snow_Temperature    = Tskin
            !---Vary Snow parameters
            !call random_number(rn)
            !Sfc%Snow_Depth          = 10+(rn*490)
            !call random_number(rn)
            !Sfc%Snow_Grain_Size     = 0.1+(rn*0.7)
            !---Fix Snow Parameters
            !Sfc%Snow_Depth          = 500
            !Sfc%Snow_Grain_Size     = 0.5
            !Sfc%Snow_Density        = 0.25
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
            !sfc(m)%Wind_Direction      = WindDir
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
         !!GeometryInfo%Year                 = Scene%scanYear
         !!---Calculate day and month
         !CALL day_month(Scene%scanYear,GeometryInfo(1)%Month,GeometryInfo(1)%Day,Scene%scanDay)
         !GeometryInfo%iFOV                 = Scene%iScanPos
         geo(m)%Sensor_Zenith_Angle  = satZenAng
         !!GeometryInfo%Sensor_Azimuth_Angle = 0
         !!GeometryInfo%Source_Zenith_Angle  = 0.
         !!GeometryInfo%Source_Azimuth_Angle = 0.
         !!---esm 14-10-28 :: trap bad angles from BUFR conversion
         !IF (ABS(Scene%SolZenAngle) >= 180.) Scene%SolZenAngle = 90.  ! turn it off
         !!---esm 14-10-20 :: turn on source for IR/NLTE-calculations 
         !geo(m)%Source_Zenith_Angle  = solZenAng

         deallocate( Plev, Play, Tlay, GasLay, cldN, Reff, WC)

      END DO Profile_Loop

   END subroutine load_atm_sfc_geo_2BCLDGPM
!
END Module read_cloudsat_module
!
