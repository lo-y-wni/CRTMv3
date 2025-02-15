!
! ODAS_BIN2NC
!
! Program to convert netCDF format ODAS files to the CRTM Binary
! format.
!
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson (benjamin.t.johnson@noaa.gov) August 9, 2024
!                       Based on ODPS_BIN2NC.f90
!

PROGRAM ODAS_BIN2NC

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage
  USE File_Utility   , ONLY: File_Exists
  USE Message_Handler, ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, &
                             Program_Message, Display_Message
  USE ODAS_Define    , ONLY: ODAS_type, &
                             Destroy_ODAS, Equal_ODAS
  USE ODAS_Binary_IO , ONLY: Read_ODAS_Binary, Write_ODAS_Binary
  USE ODAS_netCDF_IO , ONLY: Read_ODAS_netCDF, Write_ODAS_netCDF
  ! Disable implicit typing
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'ODAS_BIN2NC'

  ! ---------
  ! Variables
  ! ---------
  INTEGER :: Error_Status, n_args
  CHARACTER(256) :: NC_Filename, tmp_filename
  CHARACTER(256) :: BIN_Filename
  TYPE(ODAS_type) :: ODAS
  TYPE(ODAS_type) :: ODAS_Test


  ! Get the input and output filenames
  ! Get the filename
  n_args = COMMAND_ARGUMENT_COUNT()
  IF ( n_args > 0 ) THEN
     CALL GET_COMMAND_ARGUMENT(1, BIN_filename)
     !** automatically generate the binary filename based on the command line netcdf filename
     tmp_filename = BIN_filename(1:LEN_TRIM(BIN_filename) - 4)//".nc"
     NC_filename = TRIM(ADJUSTL(tmp_filename))
     PRINT *, "Output filename:", NC_filename
  ELSE     
     ! Output prgram header
     CALL Program_Message( PROGRAM_NAME, &
          'Program to convert netCDF format ODAS files to '//&
          'their CRTM Binary format.', &
          '$Revision$' )
     
     ! Get the input and output filenames
     WRITE( *, FMT     = '( /5x, "Enter the INPUT Binary ODAS file: " )', &
          ADVANCE = 'NO' )
     READ( *, '( a )' ) BIN_Filename
     BIN_Filename = ADJUSTL( BIN_FileNAME )
     IF ( .NOT. File_Exists( TRIM( BIN_Filename ) ) ) THEN
        CALL Display_Message( PROGRAM_NAME, &
             'File '//TRIM( BIN_Filename )//' not found.', &
             FAILURE )
        STOP
     END IF

     WRITE( *, FMT     = '( /5x, "Enter the OUTPUT netCDF ODAS file: " )', &
          ADVANCE = 'NO' )
     READ( *, '( a )' ) NC_Filename
     NC_Filename = ADJUSTL( NC_Filename )
  END IF

  ! Check that the BIN file isn't accidentally overwritten
  IF ( TRIM( BIN_Filename ) == TRIM( NC_Filename ) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Output filename is the same as the input filename!', &
                          FAILURE )
    STOP
  END IF

  ! Read the input binary file
  WRITE( *, '( /5x, "Reading Binary ODAS data ..." )' )
  Error_Status = Read_ODAS_Binary( BIN_Filename, ODAS )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error reading Binary ODAS file '//&
                          TRIM( BIN_Filename ), &
                          Error_Status )
    STOP
  END IF

  ! Write the netCDF file
  WRITE( *, '( /5x, "Writeing netCDF ODAS data ..." )' )
  Error_Status = Write_ODAS_netCDF( NC_Filename, ODAS )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error writing netCDF ODAS file '//&
                          TRIM( NC_Filename ), &
                          Error_Status )
    STOP
  END IF


  ! Test read the netCDF data file
  WRITE( *, '( /5x, "Test reading the netCDF ODAS data file ..." )' )
  Error_Status = Read_ODAS_netCDF( NC_Filename, ODAS_Test )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error reading netCDF ODAS file '//&
                          TRIM( NC_Filename ), &
                          Error_Status )
    STOP
  END IF

  ! Compare the two structures
  WRITE( *, '( /5x, "Comparing the Binary and netCDF ODAS structures ..." )' )
  Error_Status = Equal_ODAS( ODAS_Test, ODAS, Check_All=1  )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Differences found in Binary and netCDF'//&
                          'file ODAS structure comparison.', &
                          Error_Status )
  ELSE
    CALL Display_Message( PROGRAM_NAME, &
                          'Binary and netCDF file ODAS structures are equal.', &
                          INFORMATION )
  END IF

  ! Destroy the structures
  Error_Status = Destroy_ODAS( ODAS )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error destroying ODAS structure.', &
                          WARNING )
  END IF

  Error_Status = Destroy_ODAS( ODAS_Test )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'Error destroying ODAS_Test structure.', &
                          WARNING )
  END IF

END PROGRAM ODAS_BIN2NC
