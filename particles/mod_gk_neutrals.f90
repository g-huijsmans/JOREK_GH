!> Stage A: passive full-orbit hydrogenic neutrals; no fields or fluid feedback.
module mod_gk_neutrals
use mpi
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use mod_gk_neutral_trajectory, only: verify_neutral_trajectory, side_correction, &
  WALL_CROSSING, INTERIOR_ENDPOINT, UNRESOLVED, max_trace_elements
use mod_particle_sim, only: particle_sim
use mod_particle_types, only: particle_kinetic_leapfrog
use mod_find_rz_nearby, only: find_RZ_nearby
use constants, only: ATOMIC_MASS_UNIT
implicit none
private
public :: gk_neutral_state, initialize_gk_neutrals, advance_gk_neutrals
public :: report_gk_neutral_diagnostics, push_gk_neutral_ballistic
public :: get_gk_neutral_loss_totals
public :: get_gk_neutral_events, consume_gk_neutral_events, gk_neutral_event

! Fixed-size diagnostics for verifier calls only; the positive locator fast path
! is not sampled. Columns distinguish successful classifications from unresolved ones.
integer, parameter, public :: neutral_visit_bins=9
integer, parameter, public :: neutral_visit_upper(neutral_visit_bins)=[0,1,4,8,16,32,64,96,max_trace_elements]

! Pre-substep provenance for a current-step loss, not restart/particle storage.
type :: gk_neutral_event
  integer :: marker=0, element=0, substep=0
  real*8 :: x_old(3)=0.d0, v_old(3)=0.d0, st_old(2)=0.d0, dt=0.d0
end type

type :: gk_neutral_state
  integer, private :: group_num = 0, lost_count = 0
  logical, private :: enabled = .false.
  real*8, private :: lost_weight = 0.d0, lost_energy = 0.d0
  ! Transient rank-local events, never serialized or placed in a free-slot pool.
  integer, allocatable :: newly_lost_indices(:)
  type(gk_neutral_event), allocatable :: newly_lost_events(:)
  integer :: n_newly_lost_indices = 0
  ! Cumulative verification counters, separate from per-step physical losses.
  integer :: recovered_count=0, verified_wall_count=0, unresolved_count=0, max_traversed=0
  integer :: max_successful_traversed=0
  integer(int64) :: refinement_counts(3)=0_int64
  integer :: max_refinement_iterations=0
  real*8 :: refinement_residuals(2)=0.d0
  integer(int64) :: visit_histogram(neutral_visit_bins,2)=0_int64
  integer :: zero_count=0, zero_wall_count=0, zero_interior_count=0, zero_unresolved_count=0, roundoff_count=0
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
  allocate(neutrals%newly_lost_events(size(sim%groups(g)%particles)))
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
  integer :: g, j, step, old_element, ifail, locator_element, visit_bin, outcome
  type(side_correction) :: verification
  real*8 :: old_rz(2), old_st(2), old_x(3), old_v(3)

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
  if (allocated(neutrals%newly_lost_events)) then
    if (size(neutrals%newly_lost_events)/=size(sim%groups(g)%particles)) deallocate(neutrals%newly_lost_events)
  endif
  if (.not.allocated(neutrals%newly_lost_events)) allocate(neutrals%newly_lost_events(size(sim%groups(g)%particles)))
  ! Also catches a previously free slot activated with an invalid charge.
  call validate_neutral_charges(sim, g)
  select type (p => sim%groups(g)%particles)
  type is (particle_kinetic_leapfrog)
    do j=1,size(p)
      if (p(j)%i_elm <= 0) cycle
      do step=1,nsubsteps
        old_x = p(j)%x
        old_v = p(j)%v
        old_rz = p(j)%x(1:2)
        old_st = p(j)%st
        old_element = p(j)%i_elm
        if(.not.all(ieee_is_finite([old_x,old_v,old_st,dt]))) &
          call unresolved_loss(sim,neutrals,j,old_element,old_element,'nonfinite pre-substep state',0)
        call push_gk_neutral_ballistic(p(j), dt)
        if(.not.all(ieee_is_finite([p(j)%x,p(j)%v]))) &
          call unresolved_loss(sim,neutrals,j,old_element,old_element,'nonfinite ballistic endpoint',0)
        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, &
             old_rz(1), old_rz(2), old_st(1), old_st(2), old_element, &
             p(j)%x(1), p(j)%x(2), p(j)%st(1), p(j)%st(2), p(j)%i_elm, ifail, phi=p(j)%x(3))
        ! A clipped locator search is not a physical trajectory. Verify only
        ! suspicious returns, before any loss counters or event records change.
        locator_element=p(j)%i_elm
        if(locator_element<=0) then
          if(any([neutrals%zero_count,neutrals%zero_wall_count,neutrals%zero_interior_count, &
                  neutrals%zero_unresolved_count,neutrals%roundoff_count]==huge(0))) &
            call unresolved_loss(sim,neutrals,j,locator_element,old_element,'verification counter overflow',0)
          if(locator_element==0) neutrals%zero_count=neutrals%zero_count+1
          call verify_neutral_trajectory(sim%fields%node_list,sim%fields%element_list,old_element, &
            old_st,old_x,old_v,dt,verification,p(j)%x,p(j)%v)
          if(any(neutrals%refinement_counts>huge(0_int64)-int(verification%refinement_counts(:3),int64))) &
            call unresolved_loss(sim,neutrals,j,locator_element,old_element,'refinement counter overflow',0)
          neutrals%refinement_counts=neutrals%refinement_counts+int(verification%refinement_counts(:3),int64)
          neutrals%max_refinement_iterations=max(neutrals%max_refinement_iterations,verification%refinement_counts(4))
          neutrals%refinement_residuals=max(neutrals%refinement_residuals,verification%refinement_residuals)
          neutrals%max_traversed=max(neutrals%max_traversed,verification%n_visits)
          ! Nine bins and two outcomes: constant memory, no allocation or MPI here.
          visit_bin=findloc(verification%n_visits<=neutral_visit_upper,.true.,dim=1)
          if(verification%n_visits<0.or.visit_bin==0) &
            call unresolved_loss(sim,neutrals,j,locator_element,old_element,'invalid verifier visit count', &
                                 verification%n_visits)
          outcome=2
          if(verification%status==WALL_CROSSING.or.verification%status==INTERIOR_ENDPOINT) outcome=1
          if(neutrals%visit_histogram(visit_bin,outcome)==huge(0_int64)) &
            call unresolved_loss(sim,neutrals,j,locator_element,old_element,'visit histogram overflow', &
                                 verification%n_visits)
          neutrals%visit_histogram(visit_bin,outcome)=neutrals%visit_histogram(visit_bin,outcome)+1_int64
          if(outcome==1) neutrals%max_successful_traversed=max(neutrals%max_successful_traversed,verification%n_visits)
          if(verification%start_clamped.and.verification%status/=UNRESOLVED) &
            neutrals%roundoff_count=neutrals%roundoff_count+1
          select case(verification%status)
          case(INTERIOR_ENDPOINT)
            if(neutrals%recovered_count==huge(0)) &
              call unresolved_loss(sim,neutrals,j,locator_element,old_element,'recovered count overflow', &
                                   verification%n_visits)
            write(*,*) 'GK neutral false locator loss recovered: rank/group/marker=',sim%my_id,g,j, &
              ' locator element=',locator_element,' verified element=',verification%element
            write(*,*) ' locator st=',p(j)%st,' corrected st=',verification%st, &
              ' interior crossings=',verification%n_trace,' physical wall crossings=0', &
              ' elements visited=',verification%n_visits,' max geometric residual=',verification%max_residual
            p(j)%i_elm=verification%element
            p(j)%st=verification%st
            neutrals%recovered_count=neutrals%recovered_count+1
            if(locator_element==0) neutrals%zero_interior_count=neutrals%zero_interior_count+1
            cycle ! Continue any remaining substeps, without a loss event or RNG.
          case(WALL_CROSSING)
            if(neutrals%verified_wall_count==huge(0)) &
              call unresolved_loss(sim,neutrals,j,locator_element,old_element,'verified wall count overflow', &
                                   verification%n_visits)
            ! The verifier's FIRST physical wall is authoritative, even when the
            ! locator returned zero, another element, or a different wall coordinate.
            ! Retain the ballistic endpoint/velocity and original event provenance.
            p(j)%i_elm=-verification%element
            p(j)%st=verification%st
            neutrals%verified_wall_count=neutrals%verified_wall_count+1
            if(locator_element==0) neutrals%zero_wall_count=neutrals%zero_wall_count+1
          case default
            if(locator_element==0) neutrals%zero_unresolved_count=neutrals%zero_unresolved_count+1
            call capture_unresolved_trajectory(sim,g,j,step,old_element,old_st,old_x,old_v,dt, &
              locator_element,ifail,p(j),verification)
            write(*,*) 'GK neutral unresolved provenance: rank/group/marker=',sim%my_id,g,j, &
              ' substep/dt=',step,dt,' old element/st=',old_element,old_st
            write(*,*) ' old x/v=',old_x,old_v,' endpoint x/v/st=',p(j)%x,p(j)%v,p(j)%st
            call unresolved_loss(sim,neutrals,j,locator_element,old_element,verification%reason,verification%n_visits)
          end select
        endif
        ! Only a geometrically verified wall reaches accounting. No pool operation
        ! occurs here; reflection stages these events and commits after its MPI vote.
        if (p(j)%i_elm <= 0) then
          if (p(j)%i_elm<0) then
            neutrals%n_newly_lost_indices=neutrals%n_newly_lost_indices+1
            neutrals%newly_lost_indices(neutrals%n_newly_lost_indices)=j
            neutrals%newly_lost_events(neutrals%n_newly_lost_indices)= &
              gk_neutral_event(j,old_element,step,old_x,old_v,old_st,dt)
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

! Fatal geometric failures preserve a complete replay input before MPI_ABORT.
! Separate rank/marker files avoid interleaved output; never overwrite an earlier capture.
subroutine capture_unresolved_trajectory(sim,g,j,step,e,st,x,v,dt,locator,ifail,endpoint,verification)
  use phys_module, only: index_now,find_RZ_nearby_iter,find_RZ_nearby_tol
  class(particle_sim), intent(in) :: sim
  integer, intent(in) :: g,j,step,e,locator,ifail
  real*8, intent(in) :: st(2),x(3),v(3),dt
  type(particle_kinetic_leapfrog), intent(in) :: endpoint
  type(side_correction), intent(in) :: verification
  integer :: u,ios,k
  character(len=120) :: filename,stem
  write(stem,'(A,I0,A,I0,A,I0)') 'neutral_unresolved_rank',sim%my_id,'_marker',j,'_step',step
  filename=trim(stem)//'.txt'
  open(newunit=u,file=trim(filename),status='new',action='write',iostat=ios)
  if(ios/=0) then
    write(*,*) 'ERROR: cannot save unresolved trajectory state: ',trim(filename),' iostat=',ios
    return
  endif
  write(u,'(A,6I12)') 'identity ',sim%my_id,0,g,j,index_now,step
  write(u,'(A,3ES26.17)') 'mass_weight_q ',sim%groups(g)%mass,endpoint%weight,real(endpoint%q,8)
  write(u,'(A,I12)') 'locator_iterations ',find_RZ_nearby_iter
  write(u,'(A,ES26.17)') 'locator_tolerance ',find_RZ_nearby_tol
  write(u,'(A,2ES26.17)') 'time_dt ',sim%time,dt
  write(u,'(A,I12)') 'old_element ',e
  write(u,'(A,2ES26.17)') 'old_st ',st
  write(u,'(A,3ES26.17)') 'old_x ',x
  write(u,'(A,3ES26.17)') 'old_v ',v
  write(u,'(A,2I12)') 'locator ',locator,ifail
  write(u,'(A,2ES26.17)') 'locator_st ',endpoint%st
  write(u,'(A,3ES26.17)') 'ballistic_x ',endpoint%x
  write(u,'(A,3ES26.17)') 'ballistic_v ',endpoint%v
  write(u,'(A,ES26.17,2I12,ES26.17)') 'weight_charge_life_birth ', &
    endpoint%weight,endpoint%q,endpoint%i_life,real(endpoint%t_birth,8)
  write(u,'(A,A)') 'reason ',trim(verification%reason)
  write(u,'(A,4I16,2ES26.17)') 'refinement ',verification%refinement_counts,verification%refinement_residuals
  write(u,'(A,2I12,ES26.17)') 'visits_trace_residual ', &
    verification%n_visits,verification%n_trace,verification%max_residual
  do k=1,verification%n_visits
    write(u,'(A,3I12,ES26.17)') 'trace ',k,verification%trace_elements(k), &
      verification%trace_faces(k),verification%trace_fractions(k)
  enddo
  close(u)
  filename=trim(stem)//'.bin'
  open(newunit=u,file=trim(filename),status='new',access='stream',form='unformatted',iostat=ios)
  if(ios/=0) then
    write(*,*) 'ERROR: cannot save unresolved trajectory geometry: ',trim(filename),' iostat=',ios
    return
  endif
  write(u) sim%fields%node_list%n_nodes,sim%fields%element_list%n_elements
  write(u) shape(sim%fields%node_list%node(1)%x),shape(sim%fields%element_list%element(1)%size)
  do k=1,sim%fields%node_list%n_nodes
    write(u) sim%fields%node_list%node(k)%x,sim%fields%node_list%node(k)%boundary
  enddo
  do k=1,sim%fields%element_list%n_elements
    write(u) sim%fields%element_list%element(k)%size, &
      sim%fields%element_list%element(k)%vertex,sim%fields%element_list%element(k)%neighbours
  enddo
  close(u)
end subroutine

!> Passive advance retains its existing fatal-error contract. MPI_ABORT stops
!> every rank before the reflection transaction; no unresolved event is published.
subroutine unresolved_loss(sim,neutrals,j,locator_element,start_element,reason,visits)
  class(particle_sim), intent(in) :: sim
  type(gk_neutral_state), intent(inout) :: neutrals
  integer, intent(in) :: j,locator_element,start_element,visits
  character(len=*), intent(in) :: reason
  integer :: error
  if(neutrals%unresolved_count<huge(0)) neutrals%unresolved_count=neutrals%unresolved_count+1
  write(*,*) 'ERROR: GK neutral unresolved trajectory rank/group/marker=',sim%my_id,neutrals%group_num,j
  write(*,*) ' original locator element=',locator_element,' starting element=',start_element, &
             ' visits=',visits,' unresolved count=',neutrals%unresolved_count,': ',trim(reason)
  write(*,'(A,9I20)') 'GK neutral local visit histogram successful=',neutrals%visit_histogram(:,1)
  write(*,'(A,9I20)') 'GK neutral local visit histogram unresolved=',neutrals%visit_histogram(:,2)
  write(*,'(A,I0)') 'GK neutral local max successful visits=',neutrals%max_successful_traversed
  write(*,'(A,3I20,I8,2ES26.17)') 'GK neutral local refinement=',neutrals%refinement_counts, &
    neutrals%max_refinement_iterations,neutrals%refinement_residuals
  call MPI_ABORT(MPI_COMM_WORLD,1,error)
end subroutine

!> Read-only per-outer-step physical loss diagnostics, useful to independent audits.
subroutine get_gk_neutral_loss_totals(neutrals,count,weight,energy)
  type(gk_neutral_state), intent(in) :: neutrals
  integer, intent(out) :: count
  real*8, intent(out) :: weight,energy
  count=neutrals%lost_count;weight=neutrals%lost_weight;energy=neutrals%lost_energy
end subroutine

!> Read-only event snapshot: allocation is proportional to events, not all slots.
!> Only event buffer/count are public for inspection; other passive state stays private.
subroutine get_gk_neutral_events(neutrals,indices,group,ierr,events)
  type(gk_neutral_state), intent(in) :: neutrals
  integer, allocatable, intent(out) :: indices(:)
  integer, intent(out) :: group,ierr
  type(gk_neutral_event), allocatable, optional, intent(out) :: events(:)
  integer :: n
  group=neutrals%group_num; ierr=1; n=neutrals%n_newly_lost_indices
  if(n<0) return
  if(.not.allocated(neutrals%newly_lost_indices)) then
    if(n/=0) return
  else
    if(n>size(neutrals%newly_lost_indices)) return
  endif
  ! Legacy/manual index-only events still support the unchanged normal path.
  ! Missing provenance is represented by marker=0 and cannot authorize a fallback.
  if(present(events)) then
    allocate(events(n))
    if(allocated(neutrals%newly_lost_events)) then
      if(n>size(neutrals%newly_lost_events)) return
      if(n>0) events=neutrals%newly_lost_events(:n)
    endif
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
  integer :: g, j, counts(2), global_counts(2), ierr, checks(3), global_checks(3), zeros(5), global_zeros(5), max_visits
  integer :: max_successful_visits,global_refinement_max
  integer(int64) :: global_refinement(3)
  real*8 :: global_refinement_residuals(2)
  integer(int64) :: global_visits(neutral_visit_bins,2)
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
  checks=[neutrals%recovered_count,neutrals%verified_wall_count,neutrals%unresolved_count]
  call MPI_REDUCE(checks,global_checks,3,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(neutrals%max_traversed,max_visits,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  zeros=[neutrals%zero_count,neutrals%zero_wall_count,neutrals%zero_interior_count, &
         neutrals%zero_unresolved_count,neutrals%roundoff_count]
  call MPI_REDUCE(zeros,global_zeros,5,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(neutrals%visit_histogram,global_visits,2*neutral_visit_bins,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(neutrals%max_successful_traversed,max_successful_visits,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(neutrals%refinement_counts,global_refinement,3,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(neutrals%max_refinement_iterations,global_refinement_max,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(neutrals%refinement_residuals,global_refinement_residuals,2,MPI_DOUBLE_PRECISION, &
                  MPI_MAX,0,MPI_COMM_WORLD,ierr)
  if (sim%my_id == 0) then
    write(*,'(A,3I20,I8,2ES26.17)') 'GK neutral refinement attempts/failed/iterations/max/before/after=', &
      global_refinement,global_refinement_max,global_refinement_residuals
    write(*,'(A,9I5)') 'GK neutral visit histogram bin upper bounds=',neutral_visit_upper
    write(*,'(A,9I20)') 'GK neutral visit histogram successful=',global_visits(:,1)
    write(*,'(A,9I20)') 'GK neutral visit histogram unresolved=',global_visits(:,2)
    write(*,'(A,I0)') 'GK neutral max successful visits=',max_successful_visits
    write(*,*) 'GK neutral cumulative zero/wall/interior/unresolved/roundoff=',global_zeros
    write(*,*) 'GK neutral cumulative recovered/verified_wall/unresolved/max_visits=',global_checks,max_visits
    write(*,'(A,I0,2A,A,ES16.8,A,I0,A,ES24.16,A,ES24.16,A,I0,A,ES24.16,A,ES24.16)') &
         'GK neutrals group=', g, ' id=', sim%groups(g)%id, ' time=', sim%time, &
         ' active=', global_counts(1), ' weight=', global_totals(1), ' energy_J=', global_totals(2), &
         ' newly_lost=', global_counts(2), ' lost_weight=', global_totals(3), ' lost_energy_J=', global_totals(4)
  endif
end subroutine report_gk_neutral_diagnostics

end module mod_gk_neutrals
