!
! Test_Options
!
! Program to test the CRTM Options structure manipulation and
! I/O functions
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 28-Jan-2009
!                       paul.vandelst@noaa.gov
!

PROGRAM Test_Options

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE Type_Kinds
  USE Message_Handler
  USE UnitTest_Define
  USE CRTM_Options_Define
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*),  PARAMETER :: PROGRAM_NAME   = 'Test_Options'
  CHARACTER(*),  PARAMETER :: PROGRAM_VERSION_ID = &
    '$Id: Test_Options.f90 18969 2012-04-29 20:21:03Z paul.vandelst@noaa.gov $'
  ! Filenames
  CHARACTER(*), PARAMETER :: TEST_FILENAME = 'Test.CRTM_Options.bin'
  ! Dimensions
  INTEGER, PARAMETER :: N_CHANNELS = 5
  INTEGER, PARAMETER :: N_PROFILES = 2
  ! Loop number
  INTEGER, PARAMETER :: N_LOOPS = 10



  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: id
  INTEGER :: err_stat
  INTEGER :: m, n
  TYPE(UnitTest_type) :: utest
  TYPE(CRTM_Options_type) :: sopt     , r1opt(N_PROFILES)
  TYPE(CRTM_Options_type) :: sopt_copy, r1opt_copy(N_PROFILES)


  ! Output header
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to test the CRTM Options definition procedures.', &
                        '$Revision: 18969 $' )

  ! Test initialisation
  CALL UnitTest_Init(UTest)


  ! Creation test
  CALL UnitTest_Setup(utest, 'Creation Test', Caller=PROGRAM_NAME)
  DO n = 1, N_LOOPS
    ! ...Scalar
    CALL CRTM_Options_Create( sopt, N_CHANNELS )
    CALL UnitTest_Assert(utest,CRTM_Options_Associated(sopt))
    ! ...Rank1
    CALL CRTM_Options_Create( r1opt, N_CHANNELS )
    CALL UnitTest_Assert(utest,ANY(CRTM_Options_Associated(r1opt)))
  END DO
  CALL UnitTest_Report(utest)


  ! Assignment and comparison test
  CALL UnitTest_Setup(utest, 'Assignment and comparison Test', Caller=PROGRAM_NAME)
  DO n = 1, N_LOOPS
    ! ...Scalar
    sopt%Check_Input             = .FALSE.
    sopt%Use_Old_MWSSEM          = .TRUE.
    sopt%Use_Antenna_Correction  = .TRUE.
    sopt%Apply_NLTE_Correction   = .FALSE.
    sopt%RT_Algorithm_ID         = 25
    sopt%Aircraft_Pressure       = 100.0_fp
    sopt%Use_n_Streams           = .TRUE.
    sopt%n_Streams               = 2
    sopt%Include_Scattering      = .FALSE.
    sopt%Use_Emissivity          = .TRUE.
    sopt%Emissivity              = 0.99_fp
    sopt%Use_Direct_Reflectivity = .TRUE.
    sopt%Direct_Reflectivity     = 0.015_fp
    ! ...SSU component
    CALL SSU_Input_SetValue(sopt%SSU, &
                            Time=0.1234_fp, &
                            Cell_Pressure = 0.0123_fp)
    ! ...Zeeman component
    CALL Zeeman_Input_SetValue( sopt%Zeeman, &
                                Field_Strength = 1.2345_fp, &
                                Cos_ThetaB     = 2.3456_fp, &
                                Cos_PhiB       = 3.4567_fp, &
                                Doppler_Shift  = 4.5678_fp  )
    sopt_copy = sopt
    CALL UnitTest_Assert(utest,sopt_copy == sopt)
    ! ...Rank1
    DO m = 1, N_PROFILES
      r1opt(m) = sopt
    END DO
    r1opt_copy = r1opt
    CALL UnitTest_Assert(utest,ALL(r1opt_copy == r1opt))
  END DO
  CALL CRTM_Options_Inspect( sopt )
  CALL UnitTest_Report(utest)


  ! Write and Read test
  CALL UnitTest_Setup(utest, 'WriteFile and ReadFile Test', Caller=PROGRAM_NAME)
  ! ...write
  err_stat = CRTM_Options_WriteFile( TEST_FILENAME, r1opt, Quiet = .TRUE. )
  CALL UnitTest_Assert(utest, err_stat==SUCCESS)
  ! ...read
  err_stat = CRTM_Options_ReadFile( TEST_FILENAME, r1opt_copy, Quiet   = .TRUE. )
  CALL UnitTest_Assert(utest, err_stat==SUCCESS)
  ! ...check results
  CALL UnitTest_Assert(utest, ALL(r1opt_copy == r1opt))
  DO n = 1, N_PROFILES
    CALL CRTM_Options_Inspect(r1opt(n))
  END DO
  DO n = 1, N_PROFILES
    CALL CRTM_Options_Inspect(r1opt_copy(n))
  END DO
  CALL UnitTest_Report(utest)


  ! Destruction test
  CALL UnitTest_Setup(utest, 'Destruction Test', Caller=PROGRAM_NAME)
  ! ...Scalar
  CALL CRTM_Options_Destroy( sopt )
  CALL UnitTest_Assert(utest,.NOT. CRTM_Options_Associated(sopt))
  CALL UnitTest_Assert(utest,sopt%n_Channels == 0)
  CALL UnitTest_Assert(utest,.NOT. sopt%Use_Emissivity)
  CALL UnitTest_Assert(utest,.NOT. sopt%Use_Direct_Reflectivity)
  CALL UnitTest_Assert(utest,.NOT. sopt%Use_Antenna_Correction)
  CALL CRTM_Options_Inspect( sopt )
  CALL CRTM_Options_Destroy( sopt_copy )
  CALL UnitTest_Assert(utest,.NOT. CRTM_Options_Associated(sopt_copy))
  CALL UnitTest_Assert(utest,sopt_copy%n_Channels == 0)
  CALL UnitTest_Assert(utest,.NOT. sopt_copy%Use_Emissivity)
  CALL UnitTest_Assert(utest,.NOT. sopt_copy%Use_Direct_Reflectivity)
  CALL UnitTest_Assert(utest,.NOT. sopt_copy%Use_Antenna_Correction)
  ! ...Rank1
  CALL CRTM_Options_Destroy( r1opt )
  CALL UnitTest_Assert(utest,ANY(.NOT. CRTM_Options_Associated(r1opt)))
  CALL UnitTest_Assert(utest,ALL(r1opt%n_Channels == 0))
  CALL UnitTest_Assert(utest,ALL(.NOT. r1opt%Use_Emissivity))
  CALL UnitTest_Assert(utest,ALL(.NOT. r1opt%Use_Direct_Reflectivity))
  CALL UnitTest_Assert(utest,ALL(.NOT. r1opt%Use_Antenna_Correction))
  CALL CRTM_Options_Destroy( r1opt_copy )
  CALL UnitTest_Assert(utest,ANY(.NOT. CRTM_Options_Associated(r1opt_copy)))
  CALL UnitTest_Assert(utest,ALL(r1opt_copy%n_Channels == 0))
  CALL UnitTest_Assert(utest,ALL(.NOT. r1opt_copy%Use_Emissivity))
  CALL UnitTest_Assert(utest,ALL(.NOT. r1opt_copy%Use_Direct_Reflectivity))
  CALL UnitTest_Assert(utest,ALL(.NOT. r1opt_copy%Use_Antenna_Correction))
  CALL UnitTest_Report(utest)

  CALL UnitTest_Summary(utest)

END PROGRAM Test_Options
