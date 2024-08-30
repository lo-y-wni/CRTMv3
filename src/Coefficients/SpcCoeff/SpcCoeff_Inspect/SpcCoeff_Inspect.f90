!
! SpcCoeff_Inspect
!
! Program to inspect the contents of a CRTM SpcCoeff file (binary or netCDF)
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 03-Feb-2011
!                       paul.vandelst@noaa.gov
!

PROGRAM SpcCoeff_Inspect

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage
  USE File_Utility      , ONLY: File_Exists
  USE Message_Handler   , ONLY: SUCCESS, FAILURE, Program_Message, Display_Message
  USE SpcCoeff_Define   , ONLY: SpcCoeff_type, SpcCoeff_Destroy, &
                                Inspect => SpcCoeff_Inspect
  USE SpcCoeff_IO, ONLY: SpcCoeff_ReadFile
  ! Disable implicit typing
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'SpcCoeff_Inspect'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = ''

  ! ---------
  ! Variables
  ! ---------
  INTEGER :: err_stat
  CHARACTER(256) :: filename, msg
  INTEGER :: n_args
  TYPE(SpcCoeff_type) :: sc
  LOGICAL :: is_nc, is_bin

  ! Generate a string containing the SpcCoeff release for info
  WRITE(msg,'(i10)') sc%Release
  
  
  ! Output program header
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to display the contents of a CRTM '//&
                        'Binary/netCDF format R'//TRIM(ADJUSTL(msg))//' SpcCoeff '//&
                        'file to stdout.', &
                        '$Revision$' )

  ! Get the filename
  n_args = COMMAND_ARGUMENT_COUNT() 
  IF ( n_args > 0 ) THEN 
    CALL GET_COMMAND_ARGUMENT(1, filename) 
  ELSE 
    WRITE( *,FMT='(/5x,"Enter the SpcCoeff filename: ")',ADVANCE='NO' ) 
    READ( *,'(a)' ) filename 
  END IF
  filename = ADJUSTL(filename)
  IF ( .NOT. File_Exists( TRIM(filename) ) ) THEN
    msg = 'File '//TRIM(filename)//' not found.'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  END IF

  ! Check if filename ends in ".nc"
  is_nc = (INDEX(TRIM(filename), '.nc') == LEN_TRIM(filename) - 2)

  ! Check if filename ends in ".bin"
  is_bin = (INDEX(TRIM(filename), '.bin') == LEN_TRIM(filename) - 3)
  
  err_stat = FAILURE
  ! Read the data file
  If (is_bin) err_stat = SpcCoeff_ReadFile( filename, sc )
  IF (is_nc)  err_stat = SpcCoeff_ReadFile( filename, sc, netCDF=.TRUE. )

  IF ( err_stat /= SUCCESS ) THEN
    msg = 'Error reading SpcCoeff file '//TRIM(filename)
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  END IF

  ! Display the contents
  CALL Inspect( sc )

  ! Clean up
  CALL SpcCoeff_Destroy( sc )

END PROGRAM SpcCoeff_Inspect
