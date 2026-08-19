!> Required model hook for model-specific live-data metadata.
subroutine init_live_data_model(file_handle)
  implicit none

  integer, intent(in) :: file_handle

  ! Model011 currently defines no analytic input profiles to append.  Its
  ! physical background is imported from stored JOREK fields instead.
end subroutine init_live_data_model
