!> Stage A: passive full-orbit hydrogenic neutrals; no fields or fluid feedback.
module mod_gk_neutrals
use mpi
use mod_particle_sim, only: particle_sim
use mod_particle_types, only: particle_kinetic_leapfrog
use mod_find_rz_nearby, only: find_RZ_nearby
use constants, only: ATOMIC_MASS_UNIT
implicit none
private
public :: gk_neutral_state, initialize_gk_neutrals, advance_gk_neutrals
public :: report_gk_neutral_diagnostics, push_gk_neutral_ballistic
public :: get_gk_neutral_events, consume_gk_neutral_events

type :: gk_neutral_state
  integer, private :: group_num = 0, lost_count = 0
  logical, private :: enabled = .false.
  real*8, private :: lost_weight = 0.d0, lost_energy = 0.d0
  ! Transient rank-local events, never serialized or placed in a free-slot pool.
  integer, allocatable :: newly_lost_indices(:)
  integer :: n_newly_lost_indices = 0
end type gk_neutral_state

contains

!> The dedicated main supplies the single hydrogenic neutral group explicitly.
!> Fresh neutral slots remain free; populated neutrals use the common particle restart.
subroutine initialize_gk_neutrals(sim, neutral_group, neutrals)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: neutral_group
  type(gk_neutral_state), intent(out) :: neutrals
  integer :: g, ierr

  neutrals%group_num = neutral_group
  if (neutral_group > size(sim%groups)) return
  neutrals%enabled = .true.
  g = neutrals%group_num
  if (.not. allocated(sim%groups(g)%particles)) then
    write(*,*) 'ERROR: GK neutral group ', g, ' ID ', sim%groups(g)%id, ': particles not allocated'
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  endif
  ! Neutrals have no guiding centre or gyromotion: retain all three velocities.
  select type (p => sim%groups(g)%particles)
  type is (particle_kinetic_leapfrog)
    call validate_neutral_charges(sim, g)
  class default
    write(*,*) 'ERROR: GK neutral group ', g, ' ID ', sim%groups(g)%id, &
               ': expected particle_kinetic_leapfrog'
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end select
  allocate(neutrals%newly_lost_indices(size(sim%groups(g)%particles)),source=0)
  if (sim%my_id == 0) write(*,*) 'Stage A neutral group ', g, ' ID ', sim%groups(g)%id, &
                                ': particle_kinetic_leapfrog, active charges validated q=0'
end subroutine initialize_gk_neutrals

subroutine validate_neutral_charges(sim, g)
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: g
  integer :: j, ierr
  select type (p => sim%groups(g)%particles)
  type is (particle_kinetic_leapfrog)
    do j=1,size(p)
      if (p(j)%i_elm <= 0) cycle
      if (p(j)%q == 0) cycle
      write(*,*) 'ERROR: GK neutral group ', g, ' ID ', sim%groups(g)%id, &
                 ' particle ', j, ' detected charge ', p(j)%q, ': require q=0'
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    enddo
  end select
end subroutine validate_neutral_charges

!> Precondition: q=0, checked for every active marker by the group routines.
!> With zero charge the Cartesian velocity is constant, independent of E and B.
!> Leave phi continuous, as in the existing cylindrical Boris/GC pushers;
!> the field/geometry interpolation owns the toroidal-period convention.
pure subroutine push_gk_neutral_ballistic(p, dt)
  type(particle_kinetic_leapfrog), intent(inout) :: p
  real*8, intent(in) :: dt
  real*8 :: r_local, phi_local, radius, angle, c, s, vr, vphi

  r_local = p%x(1) + dt*p%v(1)
  phi_local = dt*p%v(3)
  radius = sqrt(r_local**2 + phi_local**2)
  angle = 0.d0
  ! At R=0 the basis is arbitrary; retaining the old basis avoids division by zero.
  if (radius > 0.d0) angle = atan2(phi_local, r_local)
  c = cos(angle)
  s = sin(angle)
  vr = p%v(1)
  vphi = p%v(3)
  p%x = [radius, p%x(2) + dt*p%v(2), p%x(3) + angle]
  p%v(1) = c*vr + s*vphi
  p%v(3) = -s*vr + c*vphi
end subroutine push_gk_neutral_ballistic

!> One call per outer step, using the same dt and substep count as the ions.
!> Loss counters refer only to this outer step, including steps not printed.
subroutine advance_gk_neutrals(sim, neutrals, dt, nsubsteps)
  class(particle_sim), intent(inout) :: sim
  type(gk_neutral_state), intent(inout) :: neutrals
  real*8, intent(in) :: dt
  integer, intent(in) :: nsubsteps
  integer :: g, j, step, old_element, ifail
  real*8 :: old_rz(2), old_st(2)

  if (.not. neutrals%enabled) return
  neutrals%lost_count = 0
  neutrals%lost_weight = 0.d0
  neutrals%lost_energy = 0.d0
  g = neutrals%group_num
  neutrals%n_newly_lost_indices=0
  if (allocated(neutrals%newly_lost_indices)) then
    if (size(neutrals%newly_lost_indices)/=size(sim%groups(g)%particles)) deallocate(neutrals%newly_lost_indices)
  endif
  if (.not.allocated(neutrals%newly_lost_indices)) &
    allocate(neutrals%newly_lost_indices(size(sim%groups(g)%particles)),source=0)
  ! Also catches a previously free slot activated with an invalid charge.
  call validate_neutral_charges(sim, g)
  select type (p => sim%groups(g)%particles)
  type is (particle_kinetic_leapfrog)
    do j=1,size(p)
      if (p(j)%i_elm <= 0) cycle
      do step=1,nsubsteps
        old_rz = p(j)%x(1:2)
        old_st = p(j)%st
        old_element = p(j)%i_elm
        call push_gk_neutral_ballistic(p(j), dt)
        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, &
             old_rz(1), old_rz(2), old_st(1), old_st(2), old_element, &
             p(j)%x(1), p(j)%x(2), p(j)%st(1), p(j)%st(2), p(j)%i_elm, ifail, phi=p(j)%x(3))
        ! Keep the locator's negative boundary-element index and boundary st.
        if (p(j)%i_elm <= 0) then
          if (p(j)%i_elm<0) then
            neutrals%n_newly_lost_indices=neutrals%n_newly_lost_indices+1
            neutrals%newly_lost_indices(neutrals%n_newly_lost_indices)=j
          endif
          neutrals%lost_count = neutrals%lost_count + 1
          neutrals%lost_weight = neutrals%lost_weight + p(j)%weight
          neutrals%lost_energy = neutrals%lost_energy + &
               0.5d0*p(j)%weight*sim%groups(g)%mass*ATOMIC_MASS_UNIT*dot_product(p(j)%v,p(j)%v)
          exit
        endif
      enddo
    enddo
  end select
end subroutine advance_gk_neutrals

!> Read-only event snapshot: allocation is proportional to events, not all slots.
!> Only event buffer/count are public for inspection; other passive state stays private.
subroutine get_gk_neutral_events(neutrals,indices,group,ierr)
  type(gk_neutral_state), intent(in) :: neutrals
  integer, allocatable, intent(out) :: indices(:)
  integer, intent(out) :: group,ierr
  integer :: n
  group=neutrals%group_num; ierr=1; n=neutrals%n_newly_lost_indices
  if(n<0) return
  if(.not.allocated(neutrals%newly_lost_indices)) then
    if(n/=0) return
  else
    if(n>size(neutrals%newly_lost_indices)) return
  endif
  allocate(indices(n))
  if(n>0) indices=neutrals%newly_lost_indices(:n)
  ierr=0
end subroutine

!> Consume only after global success. Scalar loss diagnostics remain wall crossings
!> before reflection, even when those markers are now active again.
subroutine consume_gk_neutral_events(neutrals)
  type(gk_neutral_state), intent(inout) :: neutrals
  neutrals%n_newly_lost_indices=0
end subroutine

!> Global active totals and losses during the most recent outer step; energy in J.
subroutine report_gk_neutral_diagnostics(sim, neutrals)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_state), intent(in) :: neutrals
  integer :: g, j, counts(2), global_counts(2), ierr
  real*8 :: totals(4), global_totals(4)

  if (.not. neutrals%enabled) return
  g = neutrals%group_num
  counts = [0, neutrals%lost_count]
  totals = [0.d0, 0.d0, neutrals%lost_weight, neutrals%lost_energy]
  select type (p => sim%groups(g)%particles)
  type is (particle_kinetic_leapfrog)
    do j=1,size(p)
      if (p(j)%i_elm <= 0) cycle
      counts(1) = counts(1) + 1
      totals(1) = totals(1) + p(j)%weight
      totals(2) = totals(2) + &
           0.5d0*p(j)%weight*sim%groups(g)%mass*ATOMIC_MASS_UNIT*dot_product(p(j)%v,p(j)%v)
    enddo
  end select
  call MPI_REDUCE(counts, global_counts, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(totals, global_totals, 4, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (sim%my_id == 0) then
    write(*,'(A,I0,2A,A,ES16.8,A,I0,A,ES24.16,A,ES24.16,A,I0,A,ES24.16,A,ES24.16)') &
         'GK neutrals group=', g, ' id=', sim%groups(g)%id, ' time=', sim%time, &
         ' active=', global_counts(1), ' weight=', global_totals(1), ' energy_J=', global_totals(2), &
         ' newly_lost=', global_counts(2), ' lost_weight=', global_totals(3), ' lost_energy_J=', global_totals(4)
  endif
end subroutine report_gk_neutral_diagnostics

end module mod_gk_neutrals
