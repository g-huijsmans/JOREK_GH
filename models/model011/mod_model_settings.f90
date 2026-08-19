!> Hard-coded settings for the electrostatic gyrokinetic Poisson model.
module mod_model_settings
  implicit none

  integer, parameter :: jorek_model = 011

  logical, parameter :: hydrodynamics = .false.
  logical, parameter :: reduced_MHD   = .false.
  logical, parameter :: full_MHD      = .false.

  logical, parameter :: with_rho        = .true.
  logical, parameter :: with_TiTe       = .false.
  logical, parameter :: with_neutrals   = .false.
  logical, parameter :: with_impurities = .false.
  logical, parameter :: with_Vpar       = .true.
  logical, parameter :: with_refluid    = .false.

  integer, parameter :: n_mod_ext = 0

  ! One equation is solved, while imported equilibrium nodes retain the
  ! seven physical fields of the reduced-MHD equilibrium.
  integer, parameter :: n_var    = 1
  integer, parameter :: n_eq_var = 7

  integer, parameter :: var_psi  = 1
  integer, parameter :: var_u    = 2
  integer, parameter :: var_zj   = 3
  integer, parameter :: var_w    = 4
  integer, parameter :: var_rho  = 5
  integer, parameter :: var_T    = 6
  integer, parameter :: var_Vpar = 7

  integer, parameter :: var_A3     = 0
  integer, parameter :: var_AR     = 0
  integer, parameter :: var_AZ     = 0
  integer, parameter :: var_uR     = 0
  integer, parameter :: var_uZ     = 0
  integer, parameter :: var_up     = 0
  integer, parameter :: var_rhon   = 0
  integer, parameter :: var_Ti     = 0
  integer, parameter :: var_Te     = 0
  integer, parameter :: var_jec    = 0
  integer, parameter :: var_jec1   = 0
  integer, parameter :: var_jec2   = 0
  integer, parameter :: var_nre    = 0
  integer, parameter :: var_rhoimp = 0

  ! Equation variable 1 is stored in physical node field var_u.
  integer, dimension(n_var), parameter :: var_index = (/ var_u /)

  ! The model supplies an analytic harmonic routine.  The compatibility
  ! element_matrix_fft entry point is only a wrapper and performs no FFT.
  logical, parameter :: unified_element_matrix = .false.
end module mod_model_settings
