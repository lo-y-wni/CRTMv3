!
! CloudCoeff_BIN2NC
!
! Program to convert Binary format CloudCoeff files to the netCDF format.
!
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson
!                       benjamin.t.johnson@noaa.gov / bjohns@ucar.edu
!       Updated:        Benjamin Johnson 8-9-2024 (bjohns@ucar.edu)
!                       Converted from CloudCoeff_BIN2NC.f90
!

PROGRAM CloudCoeff_BIN2NC

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage
  USE CRTM_Module
  USE CloudCoeff_Define   , ONLY: CloudCoeff_type
  USE CloudCoeff_IO       , ONLY: CloudCoeff_Binary_to_netCDF
  USE Endian_Utility    , ONLY: Big_Endian
  ! Disable implicit typing
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'CloudCoeff_BIN2NC'

  ! ---------
  ! Variables
  ! ---------
  INTEGER :: err_stat
  CHARACTER(256) :: msg
  CHARACTER(256) :: BIN_filename
  CHARACTER(256) :: NC_filename
  CHARACTER(256) :: answer
  CHARACTER(256) :: version_str
  
  ! Program header
  CALL CRTM_Version(version_str)
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to convert a CRTM CloudCoeff data file from Binary to netCDF format.', &
                        'CRTM Version: '//TRIM(version_str) )

  ! Get the filenames
  IF ( Big_Endian() ) THEN
     WRITE(*,FMT='(/5x,"Enter the INPUT Binary [Big Endian] CloudCoeff filename : ")', ADVANCE='NO')
  ELSE
     WRITE(*,FMT='(/5x,"Enter the INPUT Binary [Little Endian] CloudCoeff filename : ")', ADVANCE='NO')
  END IF
   READ(*,'(a)') BIN_filename
  BIN_filename = ADJUSTL(BIN_filename)
  WRITE(*,FMT='(/5x,"Enter the OUTPUT netCDF CloudCoeff filename: ")', ADVANCE='NO')
  READ(*,'(a)') NC_filename
  NC_filename = ADJUSTL(NC_filename)

  ! ...Sanity check that they're not the same
  IF ( bin_filename == nc_filename ) THEN
    msg = 'CloudCoeff netCDF and Binary filenames are the same!'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  END IF  

  ! Perform the conversion
  err_stat = CloudCoeff_Binary_to_netCDF( BIN_filename, NC_filename, Quiet=.FALSE.)
  IF ( err_stat /= SUCCESS ) THEN
    msg = 'CloudCoeff Binary -> netCDF conversion failed!'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  ELSE
    msg = 'CloudCoeff Binary -> netCDF conversion successful!'
    CALL Display_Message( PROGRAM_NAME, msg, err_stat )
  END IF
  
END PROGRAM CloudCoeff_BIN2NC
