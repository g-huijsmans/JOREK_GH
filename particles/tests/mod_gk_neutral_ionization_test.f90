module mod_gk_neutral_ionization_test
use mpi
use fruit
use, intrinsic :: ieee_arithmetic
use, intrinsic :: iso_fortran_env, only: int64
use data_structure
use mod_particle_sim, only: particle_sim,configure_particle_groups
use mod_particle_types
use mod_gk_neutral_ionization
use mod_gk_adas_ionization, only: sample_ionization_time, ionization_probability, stage_ionization_motion, &
  initialize_gk_adas_ionization, cache_ionization_profiles, advance_gk_adas_ionization, gk_adas_ionization_state
use mod_gk_neutrals
use mod_gk_neutral_trajectory
use mod_gk_neutral_reflection
use mod_pcg32_rng, only: pcg32_rng
use mod_find_rz_nearby, only: find_RZ_nearby
use phys_module, only: use_manual_random_seed, manual_seed, use_fixed_rng_value, fixed_rng_value
use mod_particle_slot_pool
use mod_particle_slot_pool_test, only: assert_slot_registry_equal, snapshot_slot_registry
use mod_fields_linear, only: jorek_fields_interp_linear
use mod_projection_helpers_test_tools, only: default_square_grid
use mod_neighbours, only: update_neighbours
use mod_interp, only: interp_RZ
use constants, only: ATOMIC_MASS_UNIT,EL_CHG
implicit none
private
public :: run_ionization_tests,run_ionization_invalid,run_ionization_event_time_tests
integer, parameter :: ion_group=1,electron_group=2,neutral_group=3
type(particle_slot_registry) :: registry
integer :: rank,tasks,interior
character(len=512) :: event_test_data_directory
real*8, parameter :: reaction_time=0.25d0
real*8 :: maximum_energy_residual=0.d0
type(type_node_list), target :: nodes
type(type_element_list), target :: elements
type, extends(jorek_fields_interp_linear) :: analytic_fields
  real*8 :: vector(3)=[1.d0,2.d0,3.d0]
contains
  procedure :: calc_EBpsiU => known_field
end type
type :: particle_snapshot
  type(particle_gc_vpar) :: ions(5),electrons(5)
  type(particle_kinetic_leapfrog) :: neutrals(3)
end type
contains


subroutine checked_report(sim,state)
  type(particle_sim), intent(in) :: sim
  type(gk_neutral_ionization_state), intent(inout) :: state
  type(particle_slot_registry) :: before
  call snapshot_slot_registry(registry,before)
  call report_gk_neutral_ionization(sim,state)
  call assert_slot_registry_equal(registry,before,'report leaves registry unchanged')
end subroutine

subroutine checked_ionize(sim,state,indices,time,ierr)
  type(particle_sim), intent(inout) :: sim
  type(gk_neutral_ionization_state), intent(inout) :: state
  integer, intent(in) :: indices(:)
  real*8, intent(in) :: time
  integer, optional, intent(out) :: ierr
  type(particle_slot_registry) :: before
  integer :: error,g,n,old_free(3)
  call validate_particle_slot_registry(sim,registry,[1,2,3],error)
  call assert_equals(0,error,'pool agrees before ionization')
  call snapshot_slot_registry(registry,before); old_free=registry%pool%n_free
  call ionize_selected_gk_neutrals(sim,state,registry,indices,time,error)
  if(error/=0) then
    call assert_slot_registry_equal(registry,before,'ionization rollback')
  else
    n=size(indices)
    call assert_true(all(registry%pool%n_free==old_free+[-n,-n,n]),'ionization pool transfers')
  endif
  call validate_particle_slot_registry(sim,registry,[1,2,3],g)
  call assert_equals(0,g,'pool agrees after ionization')
  if(present(ierr)) then
    ierr=error
  else if(error/=0) then
    call MPI_ABORT(MPI_COMM_WORLD,1,g)
  endif
end subroutine

subroutine check_restart_pool(sim)
  type(particle_sim), intent(in) :: sim
  type(particle_slot_registry) :: rebuilt
  integer :: g,j,error,next(1)
  call initialize_particle_slot_registry(sim,[1,2,3],rebuilt)
  call validate_particle_slot_registry(sim,rebuilt,[1,2,3],error)
  call assert_equals(0,error,'restart rebuilt registry valid')
  do g=1,3
    call assert_equals(registry%pool(g)%n_slots,rebuilt%pool(g)%n_slots,'restart pool size')
    call assert_equals(registry%pool(g)%n_free,rebuilt%pool(g)%n_free,'restart free count')
    do j=1,rebuilt%pool(g)%n_slots
      call assert_true(is_particle_slot_free(registry%pool(g),j).eqv. &
                       is_particle_slot_free(rebuilt%pool(g),j),'restart membership')
    enddo
    if(rebuilt%pool(g)%n_free>0) then
      call peek_particle_slots(rebuilt%pool(g),1,next)
      call assert_equals(findloc(sim%groups(g)%particles%i_elm,0,dim=1),next(1),'restart lowest free first')
    endif
  enddo
end subroutine

subroutine known_field(fields,time,i_elm,st,phi,E,B,psi,U)
  class(analytic_fields), intent(in) :: fields
  real*8, intent(in) :: time,st(2),phi
  integer, intent(in) :: i_elm
  real*8, intent(out) :: E(3),B(3),psi,U
  real*8 :: r,z
  call interp_RZ(fields%node_list,fields%element_list,i_elm,st(1),st(2),r,z)
  B=fields%vector*(1.d0+time+0.1d0*r)
  E=0.d0; psi=0.d0; U=0.d0
end subroutine

subroutine setup_mesh()
  integer :: ierr,e
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,tasks,ierr)
  call init_node_list(nodes,n_nodes_max,nodes%n_dof,n_var)
  call default_square_grid(rank,tasks,5,5,nodes,elements,ierr,0.5d0,2.5d0,-1.d0,1.d0)
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
  interior=0
  do e=1,elements%n_elements
    if(all(nodes%node(elements%element(e)%vertex)%boundary==0)) interior=e
  enddo
  call assert_true(interior>0,'interior test element exists')
end subroutine

logical function owns(j)
  integer, intent(in) :: j
  owns=rank==mod(j-1,tasks)
end function

subroutine fixture(sim,n,indices)
  type(particle_sim), intent(out) :: sim
  integer, intent(in) :: n
  integer, allocatable, intent(out) :: indices(:)
  integer :: j,g
  real*8 :: r,z
  sim%my_id=rank; sim%n_mpi=tasks; sim%time=99.d0 ! must not be used for reaction fields/birth
  allocate(analytic_fields::sim%fields)
  sim%fields%node_list=>nodes; sim%fields%element_list=>elements
  allocate(sim%groups(3))
  sim%groups%mass=2.01410178d0; sim%groups(2)%mass=5.485799d-4
  sim%groups%coupling_scheme='non'; sim%groups%ad%suffix=''
  sim%groups(1)%id='ion'; sim%groups(2)%id='ele'; sim%groups(3)%id='neu'
  allocate(particle_gc_vpar::sim%groups(1)%particles(5),sim%groups(2)%particles(5))
  allocate(particle_kinetic_leapfrog::sim%groups(3)%particles(3))
  call interp_RZ(nodes,elements,interior,0.5d0,0.5d0,r,z)
  do g=1,2
    select type(p=>sim%groups(g)%particles)
    type is(particle_gc_vpar)
      ! Protect active and unresolved/lost markers, overwrite all stale fields of free slots.
      p(1)%i_elm=interior; p(2)%i_elm=-interior
      p%q=3-2*g; p%weight=123.d0; p%i_life=7+g; p%t_birth=-2.e0
      do j=1,5
        p(j)%x=[r,z,0.2d0]; p(j)%st=0.5d0
        p(j)%vpar=999.d0; p(j)%mu=888.d0; p(j)%B_norm=777.d0
      enddo
    end select
  enddo
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,n
      if(.not.owns(j)) cycle
      p(j)%i_elm=interior; p(j)%st=0.5d0; p(j)%x=[r,z,0.7d0*j]
      p(j)%v=[1200.d0,-340.d0,560.d0]*j
      p(j)%q=0; p(j)%weight=1.d10*j; p(j)%i_life=4; p(j)%t_birth=-1.e0
    enddo
  end select
  indices=pack([(j,j=1,n)],[(owns(j),j=1,n)])
  call initialize_particle_slot_registry(sim,[1,2,3],registry)
end subroutine

subroutine snapshot(sim,saved)
  type(particle_sim), intent(in) :: sim
  type(particle_snapshot), intent(out) :: saved
  select type(p=>sim%groups(1)%particles)
  type is(particle_gc_vpar)
    saved%ions=p
  end select
  select type(p=>sim%groups(2)%particles)
  type is(particle_gc_vpar)
    saved%electrons=p
  end select
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    saved%neutrals=p
  end select
end subroutine

! Bitwise comparisons also verify rollback of intentionally injected NaNs.
logical function same_gc(a,b) result(same)
  type(particle_gc_vpar), intent(in) :: a,b
  real*8 :: x(10),y(10)
  x=[a%x,a%st,a%weight,a%vpar,a%mu,a%B_norm,dble(a%t_birth)]
  y=[b%x,b%st,b%weight,b%vpar,b%mu,b%B_norm,dble(b%t_birth)]
  same=all(transfer(x,[0_int64],10)==transfer(y,[0_int64],10)).and. &
       a%i_elm==b%i_elm.and.a%i_life==b%i_life.and.a%q==b%q
end function

logical function same_neutral(a,b) result(same)
  type(particle_kinetic_leapfrog), intent(in) :: a,b
  real*8 :: x(10),y(10)
  x=[a%x,a%st,a%weight,a%v,dble(a%t_birth)]
  y=[b%x,b%st,b%weight,b%v,dble(b%t_birth)]
  same=all(transfer(x,[0_int64],10)==transfer(y,[0_int64],10)).and. &
       a%i_elm==b%i_elm.and.a%i_life==b%i_life.and.a%q==b%q
end function

subroutine unchanged(sim,saved,label)
  type(particle_sim), intent(in) :: sim
  type(particle_snapshot), intent(in) :: saved
  character(len=*), intent(in) :: label
  integer :: j
  select type(p=>sim%groups(1)%particles)
  type is(particle_gc_vpar)
    do j=1,5
      call assert_true(same_gc(p(j),saved%ions(j)),label//' ion')
    enddo
  end select
  select type(p=>sim%groups(2)%particles)
  type is(particle_gc_vpar)
    do j=1,5
      call assert_true(same_gc(p(j),saved%electrons(j)),label//' electron')
    enddo
  end select
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,3
      call assert_true(same_neutral(p(j),saved%neutrals(j)),label//' neutral')
    enddo
  end select
end subroutine

subroutine run_ionization_tests()
  call setup_mesh()
  call run_test_case(empty_selection,'initialization and empty collective no-op')
  call run_test_case(single_conversion,'single conversion, Cartesian GC reference and full metadata')
  call run_test_case(parallel_perpendicular,'parallel and perpendicular velocities')
  call run_test_case(multiple_restart,'unequal markers, global conservation and complete restart')
  call run_test_case(invalid_transactions,'source/capacity/field/GC failures and cross-rank rollback')
  call run_test_case(invalid_groups,'explicit indices, types, charges, masses and kinetic electrons')
  if(rank==0) write(*,*) 'Maximum focused relative energy residual=',maximum_energy_residual
end subroutine

subroutine empty_selection()
  type(particle_sim) :: sim
  type(particle_snapshot) :: saved
  type(gk_neutral_ionization_state) :: state
  integer, allocatable :: indices(:)
  call fixture(sim,3,indices)
  call snapshot(sim,saved)
  call initialize_gk_neutral_ionization(sim,ion_group,electron_group,neutral_group,state)
  call assert_true(state%initialized,'valid initialization')
  call checked_ionize(sim,state,[integer::],reaction_time)
  call checked_report(sim,state)
  call unchanged(sim,saved,'empty selection')
  call assert_true(all(state%global_counts==0).and.all(state%global_totals==0.d0),'zero no-op diagnostics')
end subroutine

subroutine single_conversion()
  type(particle_sim) :: sim
  type(particle_snapshot) :: saved,after
  type(gk_neutral_ionization_state) :: state
  type(particle_gc_vpar) :: ion,electron
  type(particle_kinetic_leapfrog) :: source
  integer, allocatable :: indices(:)
  real*8 :: b(3),E(3),psi,U,scale,d(3),cart(2),expected(3),energy,bgc
  integer :: ierr
  call fixture(sim,1,indices)
  call snapshot(sim,saved)
  call initialize_gk_neutral_ionization(sim,ion_group,electron_group,neutral_group,state)
  call checked_ionize(sim,state,indices,reaction_time)
  call snapshot(sim,after)
  call assert_true(same_gc(saved%ions(1),after%ions(1)),'active ion protected')
  call assert_true(same_gc(saved%ions(2),after%ions(2)),'lost ion protected')
  call assert_true(same_gc(saved%electrons(1),after%electrons(1)),'active electron protected')
  call assert_true(same_gc(saved%electrons(2),after%electrons(2)),'lost electron protected')
  if(owns(1)) then
    source=saved%neutrals(1); ion=after%ions(3); electron=after%electrons(3)
    call assert_equals(0,after%neutrals(1)%i_elm,'source inactive')
    source%i_elm=0
    call assert_true(same_neutral(source,after%neutrals(1)),'only source element changed')
    source=saved%neutrals(1)
    call assert_equals(2,count(after%ions%i_elm>0),'one ion activated')
    call assert_equals(2,count(after%electrons%i_elm>0),'one electron activated')
    call assert_equals(1,int(ion%q),'ion charge'); call assert_equals(-1,int(electron%q),'electron charge')
    call assert_equals(source%weight,ion%weight,0.d0,'ion weight')
    call assert_equals(source%weight,electron%weight,0.d0,'electron weight')
    call assert_equals(0.d0,electron%vpar,0.d0,'cold electron vpar')
    call assert_equals(0.d0,electron%mu,0.d0,'cold electron mu')
    call assert_true(all(electron%x==source%x).and.all(electron%st==source%st).and. &
                     electron%i_elm==source%i_elm,'electron reaction position')
    call assert_true(ion%i_elm>0.and.all(ion%st>=0.d0).and.all(ion%st<=1.d0),'valid GC location')
    call sim%fields%calc_EBpsiU(reaction_time,source%i_elm,source%st,source%x(3),E,b,psi,U)
    call assert_equals(norm2(b),electron%B_norm,1.d-14,'electron source B at explicit time')
    ! Independent Cartesian position reference; explicit component cross product,
    ! then rotate the physical displacement, never add a length to an angle.
    scale=sim%groups(1)%mass*ATOMIC_MASS_UNIT/(EL_CHG*sum(b*b))
    d=scale*[source%v(2)*b(3)-source%v(3)*b(2), &
             source%v(3)*b(1)-source%v(1)*b(3),source%v(1)*b(2)-source%v(2)*b(1)]
    cart=[source%x(1)*cos(source%x(3)),source%x(1)*sin(source%x(3))]+ &
         [d(1)*cos(source%x(3))-d(3)*sin(source%x(3)), &
          d(1)*sin(source%x(3))+d(3)*cos(source%x(3))]
    expected=[norm2(cart),source%x(2)+d(2),atan2(cart(2),cart(1))]
    call assert_equals(0.d0,maxval(abs(ion%x-expected)),2.d-14,'Cartesian GC reference')
    call assert_true(abs(ion%x(3)-(source%x(3)+d(3)))>1.d-8,'detect direct length-to-phi bug')
    bgc=sqrt(14.d0)*(1.d0+reaction_time+0.1d0*ion%x(1))
    call assert_equals(bgc,ion%B_norm,2.d-12,'B reevaluated at GC')
    energy=(ion%vpar**2+2.d0*ion%mu*ion%B_norm)/sum(source%v**2)
    call assert_equals(1.d0,energy,3.d-14,'ion energy equals source energy')
    call assert_equals(saved%ions(3)%i_life+1,ion%i_life,'ion destination lifetime')
    call assert_equals(saved%electrons(3)%i_life+1,electron%i_life,'electron destination lifetime')
    call assert_equals(reaction_time,dble(ion%t_birth),0.d0,'ion explicit birth')
    call assert_equals(reaction_time,dble(electron%t_birth),0.d0,'electron explicit birth')
  endif
  call checked_report(sim,state)
  call assert_equals(1,state%global_counts(2),'one physical conversion')
  call assert_equals(0.d0,state%global_totals(7),0.d0,'zero created charge')
  call checked_ionize(sim,state,indices,reaction_time,ierr)
  call assert_true(ierr/=0,'double conversion rejected globally')
  call unchanged(sim,after,'double conversion rollback')
end subroutine

subroutine parallel_perpendicular()
  type(particle_sim) :: sim
  type(particle_snapshot) :: before,after
  type(gk_neutral_ionization_state) :: state
  integer, allocatable :: indices(:)
  integer :: k
  real*8 :: displacement
  do k=1,2
    call fixture(sim,1,indices)
    select type(f=>sim%fields)
    type is(analytic_fields)
      f%vector=[0.d0,0.d0,2.d0]
    end select
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      p(1)%v=[0.d0,0.d0,1000.d0]
      if(k==2) p(1)%v=[0.d0,1000.d0,0.d0]
    end select
    call snapshot(sim,before)
    call initialize_gk_neutral_ionization(sim,ion_group,electron_group,neutral_group,state)
    call checked_ionize(sim,state,indices,reaction_time)
    call snapshot(sim,after)
    if(.not.owns(1)) cycle
    if(k==1) then
      call assert_equals(0.d0,maxval(abs(after%ions(3)%x-before%neutrals(1)%x)),0.d0,'parallel no displacement')
      call assert_equals(0.d0,after%ions(3)%mu,0.d0,'parallel mu zero')
      call assert_equals(1000.d0,after%ions(3)%vpar,1.d-12,'parallel energy in vpar')
    else
      displacement=sim%groups(1)%mass*ATOMIC_MASS_UNIT*1000.d0/ &
                   (EL_CHG*2.d0*(1.d0+reaction_time+0.1d0*before%neutrals(1)%x(1)))
      call assert_equals(before%neutrals(1)%x(1)+displacement,after%ions(3)%x(1),1.d-14,'perpendicular shift sign')
      call assert_equals(0.d0,after%ions(3)%vpar,0.d0,'perpendicular vpar zero')
      call assert_equals(5.d5,after%ions(3)%mu*after%ions(3)%B_norm,1.d-8,'perpendicular energy in mu B')
    endif
  enddo
end subroutine

subroutine multiple_restart()
  use hdf5, only: H5open_f
  use mod_particle_io, only: write_simulation_hdf5,read_simulation_hdf5
  use phys_module, only: n_part_groups,part_groups_in_use,part_group_configs
  use mod_particle_group_id, only: matching_part_config_indices
  type(particle_sim) :: sim,restored
  type(particle_snapshot) :: saved
  type(gk_neutral_ionization_state) :: state
  integer, allocatable :: indices(:)
  integer :: j,g,ierr,unit
  real*8 :: expected,totals(10)
  character(len=128) :: filename
  call fixture(sim,3,indices)
  call initialize_gk_neutral_ionization(sim,ion_group,electron_group,neutral_group,state)
  call checked_ionize(sim,state,indices,reaction_time)
  call checked_report(sim,state)
  expected=0.5d0*sim%groups(3)%mass*ATOMIC_MASS_UNIT*sum([1200.d0,-340.d0,560.d0]**2)*1.d10*36.d0
  call assert_equals(3,state%global_counts(1),'three unique global sources')
  call assert_equals(3,state%global_counts(2),'three global pairs')
  call assert_true(all(state%global_totals(1:3)==6.d10),'equal source/ion/electron weights')
  call assert_equals(1.d0,state%global_totals(4)/expected,3.d-14,'source energy')
  call assert_equals(1.d0,state%global_totals(5)/expected,3.d-14,'ion energy')
  call assert_equals(0.d0,state%global_totals(6),0.d0,'zero electron energy')
  call assert_true(all(state%global_totals(7:9)==0.d0),'zero charge and nuclear residuals')
  maximum_energy_residual=max(maximum_energy_residual,abs(state%global_totals(10))/expected)
  totals=state%global_totals
  call checked_report(sim,state)
  call assert_true(all(state%global_totals==totals),'report does not double reduce')
  call snapshot(sim,saved)
  call export_products(saved,indices)
  n_part_groups=3; part_groups_in_use(1:3)=sim%groups%id
  write(filename,'(A,I0,A)') 'ionization_restart_np',tasks,'.h5'
  call H5open_f(ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call write_simulation_hdf5(sim,trim(filename))
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  do g=1,3
    matching_part_config_indices(g)=g
    part_group_configs(g)%id=sim%groups(g)%id
    part_group_configs(g)%mass=sim%groups(g)%mass
    part_group_configs(g)%Z=1; part_group_configs(g)%coupling_scheme='non'
    part_group_configs(g)%n_particles=size(sim%groups(g)%particles)*tasks
    part_group_configs(g)%type='particle_gc_vpar'
  enddo
  part_group_configs(3)%type='particle_kinetic_leapfrog'
  restored%my_id=rank; restored%n_mpi=tasks
  call restored%allocate_groups(3)
  call configure_particle_groups(restored)
  call read_simulation_hdf5(restored,trim(filename),test_in=.true.)
  do g=1,3
    call assert_equals(sim%groups(g)%id,restored%groups(g)%id,'restart group ID')
    call assert_equals(codify_particle_type(sim%groups(g)%particles), &
                       codify_particle_type(restored%groups(g)%particles),'restart dynamic type')
  enddo
  call unchanged(restored,saved,'complete restart')
  call check_restart_pool(restored)
  call unchanged(sim,saved,'report/restart write does not mutate particles')
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(rank==0) then
    open(newunit=unit,file=trim(filename),status='old',iostat=ierr)
    if(ierr==0) close(unit,status='delete')
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine

! Export by global physical marker identity, independent of destination slot/rank.
! The scratch checker compares these complete products across MPI decompositions.
subroutine export_products(saved,indices)
  type(particle_snapshot), intent(in) :: saved
  integer, intent(in) :: indices(:)
  type(particle_gc_vpar) :: p
  real*8 :: local(13,2,3),global(13,2,3)
  integer :: k,j,g,ierr,unit
  character(len=80) :: filename
  local=0.d0
  do k=1,size(indices)
    j=indices(k)
    do g=1,2
      p=saved%ions(k+2)
      if(g==2) p=saved%electrons(k+2)
      local(:,g,j)=[p%x,p%st,p%weight,p%vpar,p%mu,p%B_norm,dble(p%t_birth), &
                    dble(p%i_elm),dble(p%i_life),dble(p%q)]
    enddo
  enddo
  call MPI_REDUCE(local,global,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if(rank==0) then
    write(filename,'(A,I0,A)') 'physical_products_np',tasks,'.txt'
    open(newunit=unit,file=trim(filename),status='replace')
    do j=1,3
      do g=1,2
        write(unit,'(2I4,13ES26.17)') j,g,global(:,g,j)
      enddo
    enddo
    close(unit)
  endif
end subroutine

subroutine invalid_transactions()
  type(particle_sim) :: sim
  type(particle_snapshot) :: saved
  type(gk_neutral_ionization_state) :: state
  integer, allocatable :: indices(:)
  integer :: mode,ierr,bad_rank,e,boundary_element
  real*8 :: nan,r,z,max_r
  character(len=64) :: label
  nan=ieee_value(0.d0,ieee_quiet_nan)
  do mode=1,18
    call fixture(sim,3,indices)
    call initialize_gk_neutral_ionization(sim,ion_group,electron_group,neutral_group,state)
    bad_rank=min(1,tasks-1) ! rank 0 has valid work when rank 1 rejects
    if(rank==bad_rank) then
      indices=[bad_rank+1]
      select case(mode)
      case(1); indices=[99]
      case(2); sim%groups(3)%particles(indices(1))%i_elm=0
      case(3); indices=[indices(1),indices(1)]
      case(4); sim%groups(1)%particles(3:5)%i_elm=-interior
      case(5); sim%groups(2)%particles(3:5)%i_elm=-interior
      case(6:12)
        select type(p=>sim%groups(3)%particles)
        type is(particle_kinetic_leapfrog)
          select case(mode)
          case(6); p(indices(1))%q=1
          case(7); p(indices(1))%weight=0.d0
          case(8); p(indices(1))%weight=-1.d0
          case(9); p(indices(1))%x(1)=nan
          case(10); p(indices(1))%v(1)=nan
          case(11); p(indices(1))%st(1)=1.1d0
          case(12)
            ! Source just inside the outer radial wall, positive radial GC shift.
            max_r=0.d0; boundary_element=0
            do e=1,elements%n_elements
              call interp_RZ(nodes,elements,e,0.5d0,0.5d0,r,z)
              if(r<=max_r) cycle
              max_r=r; boundary_element=e
            enddo
            p(indices(1))%i_elm=boundary_element
            p(indices(1))%st=[1.d0-1.d-6,0.5d0]
            call interp_RZ(nodes,elements,boundary_element,p(indices(1))%st(1),0.5d0,r,z)
            p(indices(1))%x(:2)=[r,z]
            p(indices(1))%v=[0.d0,1.d5,0.d0]
          end select
        end select
      case(13:14)
        select type(f=>sim%fields)
        type is(analytic_fields)
          f%vector=0.d0
          if(mode==14) f%vector=nan
        end select
      case(15); sim%groups(3)%particles(indices(1))%i_elm=elements%n_elements+1
      case(16); sim%groups(3)%particles(indices(1))%st(1)=nan
      case(17); sim%groups(3)%particles(indices(1))%weight=nan
      case(18); sim%groups(1)%particles(3)%i_life=huge(0)
      end select
    endif
    call initialize_particle_slot_registry(sim,[1,2,3],registry) ! finalize mutated test fixture
    call snapshot(sim,saved)
    call checked_ionize(sim,state,indices,reaction_time,ierr)
    write(label,'(A,I0)') 'global rollback case ',mode
    call assert_true(ierr/=0,trim(label)//' rejected on every rank')
    call unchanged(sim,saved,trim(label))
    call assert_equals(0,state%counts(2),'no committed conversions')
    call assert_true(all(state%totals==0.d0),'no committed conservation increments')
    call assert_equals(merge(1,0,rank==bad_rank),state%counts(3),'only originating rank counts failure')
  enddo
end subroutine

subroutine invalidate_group(sim,mode)
  type(particle_sim), intent(inout) :: sim
  integer, intent(in) :: mode
  integer :: g
  select case(mode)
  case(1:3)
    g=mode
    deallocate(sim%groups(g)%particles)
    if(g==3) then
      allocate(particle_gc_vpar::sim%groups(g)%particles(3))
    else
      allocate(particle_kinetic_leapfrog::sim%groups(g)%particles(5))
    endif
  case(4:5)
    select type(p=>sim%groups(mode-3)%particles)
    type is(particle_gc_vpar)
      p(1)%q=0
    end select
  case(6)
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      where(p%i_elm>0) p%q=1
    end select
  case(7); sim%groups(3)%mass=1.d0
  case(8); sim%groups(2)%mass=0.d0
  case(9); sim%groups(1)%mass=ieee_value(0.d0,ieee_quiet_nan)
  case(10); sim%groups(3)%coupling_scheme='ics'
  end select
end subroutine

subroutine invalid_groups()
  type(particle_sim) :: sim
  type(gk_neutral_ionization_state) :: state
  integer, allocatable :: indices(:)
  integer :: mode,ierr,eg,ng
  do mode=1,13
    call fixture(sim,3,indices)
    eg=electron_group; ng=neutral_group
    if(mode<=10) call invalidate_group(sim,mode)
    if(mode==11) eg=4
    if(mode==12) ng=0
    if(mode==13) ng=ion_group
    call initialize_gk_neutral_ionization(sim,ion_group,eg,ng,state,ierr)
    call assert_true(ierr/=0.and..not.state%initialized,'invalid initialization rejected globally')
  enddo
end subroutine

subroutine run_ionization_invalid(mode)
  character(len=*), intent(in) :: mode
  type(particle_sim) :: sim
  type(gk_neutral_ionization_state) :: state
  integer, allocatable :: indices(:)
  integer :: eg
  call setup_mesh()
  call fixture(sim,3,indices)
  eg=electron_group
  if(mode=='electron_type') call invalidate_group(sim,2)
  if(mode=='electron_absent') eg=4
  call initialize_gk_neutral_ionization(sim,ion_group,eg,neutral_group,state)
  if(rank==min(1,tasks-1)) sim%groups(2)%particles(3:5)%i_elm=-interior
  call initialize_particle_slot_registry(sim,[1,2,3],registry)
  call ionize_selected_gk_neutrals(sim,state,registry,indices,reaction_time)
end subroutine
subroutine run_ionization_event_time_tests(directory)
  character(len=*),intent(in)::directory
  event_test_data_directory=directory
  call setup_mesh()
  call run_test_case(curved_veto_transport,'outside GC veto followed by curved exit and reentry')
  call run_test_case(gc_boundary_selection_tests,'classified GC veto, mixed MPI events, errors and resumed wall transport')
  call run_test_case(waiting_time_statistics,'same-draw exponential waiting times and conditional CDF')
  call run_test_case(individual_event_times,'staged reaction states and individual physical birth/field times')
  call run_test_case(event_array_rollback,'per-event validation and MPI-global rollback')
  call run_test_case(ballistic_event_order,'wall starts, substeps, zero time and first wall precedence')
  call run_test_case(curved_event_order,'curved crossing with interior endpoint and existing reflection handoff')
  call run_test_case(duplicate_face_tests,'exact duplicate-node cubic faces and malformed connectivity rejection')
  call run_test_case(selection_transaction_tests,'selection RNG, source/pool/diagnostic rollback and wall-first handoff')
end subroutine

subroutine waiting_time_statistics()
  type(pcg32_rng) :: rng,replay
  integer, parameter :: samples=200000
  type(particle_kinetic_leapfrog) :: flight
  integer :: i,k,status,nselected,mismatches,cdf(3),position_cdf(3),bad
  real*8 :: nu,dt,u(1),v(1),tau,p,expected,se,thresholds(3),fraction,mean_tau
  logical :: selected
  call rng%initialize(1,20261103,tasks,100000000+rank)
  replay=rng
  do i=1,100
    call rng%next(u);call replay%next(v)
    call assert_true(all(u==v),'PCG sequence replay')
  enddo
  thresholds=[.25d0,.5d0,.75d0];dt=1.d0;nu=1.d0;p=ionization_probability(nu,dt)
  nselected=0;mismatches=0;cdf=0;position_cdf=0;mean_tau=0.d0;bad=0
  do i=1,samples
    call rng%next(u)
    call sample_ionization_time(nu,dt,u(1),tau,selected,status)
    if(status/=0)bad=bad+1
    if(selected.neqv.(u(1)<p))mismatches=mismatches+1
    if(.not.selected)cycle
    nselected=nselected+1;mean_tau=mean_tau+tau
    if(tau<0.d0.or.tau>=dt)bad=bad+1
    flight%x=[1.13d0,.1d0,.7d0];flight%v=[0.d0,.2d0,0.d0]
    call push_gk_neutral_ballistic(flight,tau)
    do k=1,3
      if(flight%x(2)<.1d0+.2d0*thresholds(k))position_cdf(k)=position_cdf(k)+1
      if(tau<thresholds(k))cdf(k)=cdf(k)+1
    enddo
  enddo
  call assert_true(all(position_cdf==cdf),'no-wall ballistic reaction-position CDF matches waiting-time CDF')
  call assert_equals(0,mismatches,'exact Stage1 Bernoulli decisions from same uniforms')
  call assert_equals(0,bad,'valid finite waiting times')
  se=sqrt(samples*p*(1.d0-p));call assert_true(abs(nselected-samples*p)<6.d0*se,'exponential survival count within 6 sigma')
  do k=1,3
    fraction=(1.d0-exp(-thresholds(k)))/p;se=sqrt(nselected*fraction*(1.d0-fraction))
    call assert_true(abs(cdf(k)-nselected*fraction)<6.d0*se,'conditional event-time CDF within 6 sigma')
  enddo
  expected=1.d0-1.d0/(exp(1.d0)-1.d0)
  call assert_true(abs(mean_tau/nselected-expected)<6.d0/sqrt(12.d0*nselected),'conditional mean waiting time')
  write(*,'(A,2I12,5ES24.15)')'WAIT_STATS rank/events/P/mean/expected/cdf25/cdf75=', &
    rank,nselected,p,mean_tau/nselected,expected,real(cdf(1),8)/nselected,real(cdf(3),8)/nselected
  do k=1,5
    select case(k)
    case(1);nu=0.d0
    case(2);nu=1.d-12
    case(3);nu=1.d0
    case(4);nu=1.d12
    case(5);nu=1.d300
    end select
    p=ionization_probability(nu,dt)
    do i=1,5
      select case(i)
      case(1);u=0.d0
      case(2);u=nearest(1.d0,-1.d0)
      case(3);u=tiny(1.d0)
      case(4);u=min(p,nearest(1.d0,-1.d0))
      case(5);u=max(0.d0,nearest(p,-1.d0))
      end select
      call sample_ionization_time(nu,dt,u(1),tau,selected,status)
      call assert_equals(0,status,'limiting waiting-time input')
      call assert_true(selected.eqv.(u(1)<p),'same decision at limiting uniforms')
      if(selected)call assert_true(ieee_is_finite(tau).and.tau>=0.d0.and.tau<dt,'limiting event stays in interval')
    enddo
  enddo
end subroutine

subroutine individual_event_times()
  type(particle_sim)::sim
  type(particle_snapshot)::saved
  type(gk_neutral_ionization_state)::state
  type(particle_kinetic_leapfrog),allocatable::sources(:)
  type(side_correction)::hit
  integer,allocatable::indices(:),slots(:)
  real*8,allocatable::times(:)
  real*8::tau,wall,B(3),E(3),psi,U
  integer::k,j,status
  call fixture(sim,3,indices)
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    if(owns(1))then
      p(1)%x=[1.13d0,-1.d0,.7d0];p(1)%v=[100.d0,10000.d0,50.d0]
      call locate_test_particle(p(1));p(1)%st(2)=-8.d-20
    endif
  end select
  call snapshot(sim,saved)
  call initialize_gk_neutral_ionization(sim,1,2,3,state)
  allocate(sources(size(indices)),times(size(indices)),slots(size(indices)))
  call peek_particle_slots(registry%pool(1),size(indices),slots)
  do k=1,size(indices)
    j=indices(k);tau=(j+1)*1.d-7;times(k)=.1d0+tau
    call stage_ionization_motion(sim,saved%neutrals(j),tau,5.d-8,sources(k),hit,wall)
    call assert_equals(INTERIOR_ENDPOINT,hit%status,'staged interior flight')
  enddo
  call unchanged(sim,saved,'staging does not touch particles')
  call ionize_selected_gk_neutrals(sim,state,registry,indices,.100001d0,status, &
    reaction_times=times,reaction_sources=sources)
  call assert_equals(0,status,'individual event collective commits')
  select type(p=>sim%groups(1)%particles)
  type is(particle_gc_vpar)
    do k=1,size(indices)
      j=indices(k)
      call assert_true(p(slots(k))%t_birth==real(times(k),4),'ion individual birth time')
      call sim%fields%calc_EBpsiU(times(k),p(slots(k))%i_elm,p(slots(k))%st,p(slots(k))%x(3),E,B,psi,U)
      call assert_equals(norm2(B),p(slots(k))%B_norm,1.d-12,'field evaluated at individual event time')
      call assert_equals(sources(k)%weight,p(slots(k))%weight,0.d0,'event weight preserved')
    enddo
  end select
  select type(p=>sim%groups(2)%particles)
  type is(particle_gc_vpar)
    do k=1,size(indices)
      call assert_true(p(slots(k))%t_birth==real(times(k),4),'electron individual birth time')
      call assert_true(all(p(slots(k))%x==sources(k)%x),'electron at physical reaction position')
      call assert_true(p(slots(k))%vpar==0.d0.and.p(slots(k))%mu==0.d0,'secondary remains zero-energy')
    enddo
  end select
  call validate_particle_slot_registry(sim,registry,[1,2,3],status)
  call assert_equals(0,status,'event slot ownership remains consistent')
  call assert_true(all(abs(state%totals(7:10))<1.d-12),'event charge and energy accounting')
end subroutine

subroutine event_array_rollback()
  type(particle_sim)::sim
  type(particle_snapshot)::saved
  type(gk_neutral_ionization_state)::state
  type(particle_slot_registry)::before
  type(particle_kinetic_leapfrog),allocatable::sources(:)
  integer,allocatable::indices(:)
  real*8,allocatable::times(:)
  integer::which,k,status,bad
  bad=min(1,tasks-1)
  do which=1,5
    call fixture(sim,3,indices);call snapshot(sim,saved);call snapshot_slot_registry(registry,before)
    call initialize_gk_neutral_ionization(sim,1,2,3,state)
    allocate(sources(size(indices)),times(size(indices)))
    do k=1,size(indices)
      sources(k)=saved%neutrals(indices(k));times(k)=.2d0
    enddo
    if(rank==bad)then
      select case(which)
      case(1);times(1)=ieee_value(1.d0,ieee_quiet_nan)
      case(2);sources(1)%weight=sources(1)%weight*2.d0
      case(3);sources(1)%v=sources(1)%v*2.d0
      case(4);sources(1)%i_elm=0
      case(5);times=[real*8::]
      end select
    endif
    call ionize_selected_gk_neutrals(sim,state,registry,indices,.25d0,status,reaction_times=times,reaction_sources=sources)
    call assert_true(status/=0,'invalid staged event rejected globally')
    call unchanged(sim,saved,'event-array global rollback')
    call assert_slot_registry_equal(registry,before,'event-array pool rollback')
    call assert_true(state%counts(2)==0.and.all(state%totals==0.d0),'no committed event diagnostics')
    deallocate(sources,times)
  enddo
end subroutine

subroutine locate_test_particle(p)
  type(particle_kinetic_leapfrog),intent(inout)::p
  integer::status
  call find_RZ_nearby(nodes,elements,0.d0,0.d0,0.d0,0.d0,0,p%x(1),p%x(2),p%st(1),p%st(2),p%i_elm,status,phi=p%x(3))
  call assert_true(status==0.and.p%i_elm>0,'independent test source location')
end subroutine

subroutine ballistic_event_order()
  type(particle_sim)::sim
  type(particle_kinetic_leapfrog)::original,source,subdivided,expected
  type(side_correction)::hit
  integer,allocatable::indices(:)
  integer::which
  real*8::wall,tau,cart0(3),cartv(3),cart1(3),phi
  call fixture(sim,0,indices)
  do which=1,4
    original%x=[.5d0,.1d0,7.3d0];original%v=[10000.d0,100.d0,30000.d0]
    if(which==2)original%x(1)=.5d0+1.d-12
    if(which==3)original%x=[1.d0,.1d0,7.3d0]
    call locate_test_particle(original);original%q=0;original%weight=1.d10
    tau=1.d-6;if(which==4)tau=0.d0
    call stage_ionization_motion(sim,original,tau,1.d-6,source,hit,wall)
    call assert_equals(INTERIOR_ENDPOINT,hit%status,'inward-wall/interior event occurs before wall')
    call stage_ionization_motion(sim,original,tau,2.5d-7,subdivided,hit,wall)
    call assert_equals(INTERIOR_ENDPOINT,hit%status,'several neutral substeps')
    call assert_true(norm2(source%x-subdivided%x)<2.d-14,'substeps reproduce reaction position')
    call assert_true(norm2(source%v-subdivided%v)<2.d-10,'substeps reproduce rotated velocity')
    phi=original%x(3);cart0=[original%x(1)*cos(phi),original%x(1)*sin(phi),original%x(2)]
    cartv=[original%v(1)*cos(phi)-original%v(3)*sin(phi),original%v(1)*sin(phi)+original%v(3)*cos(phi),original%v(2)]
    cart1=[source%x(1)*cos(source%x(3)),source%x(1)*sin(source%x(3)),source%x(2)]
    call assert_true(norm2(cart1-cart0-tau*cartv)<2.d-14,'physical Cartesian reaction position')
  enddo
  original%x=[1.13d0,.99d0,0.d0];original%v=[0.d0,10000.d0,0.d0];call locate_test_particle(original)
  call stage_ionization_motion(sim,original,1.d-5,3.d-6,source,hit,wall)
  call assert_equals(WALL_CROSSING,hit%status,'wall precedes sampled event')
  call assert_equals(1.d-6,wall,1.d-14,'first wall contact time')
  original%x=[1.13d0,1.d0,0.d0];call locate_test_particle(original)
  call stage_ionization_motion(sim,original,0.d0,1.d-6,source,hit,wall)
  call assert_equals(WALL_CROSSING,hit%status,'outward wall and zero-time tie belongs to wall')
end subroutine

subroutine curved_event_order()
  type(particle_sim)::sim
  type(gk_neutral_state)::passive
  type(gk_neutral_reflection_state)::reflection
  type(particle_kinetic_leapfrog)::original,source,endpoint
  type(side_correction)::hit
  integer,allocatable::indices(:),forced(:)
  real*8,allocatable::saved_x(:,:,:,:)
  real*8::z,wall,dt,alpha,energy
  integer::j,k,status
  allocate(saved_x(nodes%n_nodes,n_coord_tor,n_degrees,2));alpha=.2d0
  do j=1,nodes%n_nodes
    saved_x(j,:,:,:)=nodes%node(j)%x
    ! R=Rbase+alpha*Z^2, with exact chain-rule nodal derivatives.
    z=nodes%node(j)%x(1,1,2)
    nodes%node(j)%x(1,1,1)=nodes%node(j)%x(1,1,1)+alpha*z*z
    do k=2,3
      nodes%node(j)%x(1,k,1)=nodes%node(j)%x(1,k,1)+2.d0*alpha*z*nodes%node(j)%x(1,k,2)
    enddo
    nodes%node(j)%x(1,4,1)=nodes%node(j)%x(1,4,1)+2.d0*alpha* &
      (nodes%node(j)%x(1,2,2)*nodes%node(j)%x(1,3,2)+z*nodes%node(j)%x(1,4,2))
  enddo
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
  call fixture(sim,1,indices)
  original%x=[2.52d0,-.4d0,0.d0];original%v=[0.d0,8000.d0,0.d0]
  original%q=0;original%weight=1.d10;original%i_life=4;original%t_birth=-1.e0
  call locate_test_particle(original)
  dt=1.d-4;endpoint=original;call push_gk_neutral_ballistic(endpoint,dt)
  call locate_test_particle(endpoint)
  call assert_true(endpoint%i_elm>0,'curved-exit test endpoint is inside')
  call stage_ionization_motion(sim,original,dt,dt,source,hit,wall)
  call assert_equals(WALL_CROSSING,hit%status,'first curved wall crossing despite inside endpoint')
  call assert_equals((.4d0-sqrt(.1d0))/8000.d0,wall,2.d-13,'analytic curved-wall first contact')
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    if(owns(1))p(1)=original
  end select
  forced=indices
  call initialize_gk_neutrals(sim,3,passive)
  use_manual_random_seed=.true.;manual_seed=20261103
  call initialize_gk_neutral_reflection(sim,3,.true.,reflection)
  call advance_gk_neutrals(sim,passive,dt,1,verify_indices=forced)
  if(owns(1))call assert_equals(1,passive%n_newly_lost_indices,'wall-first candidate reaches existing loss buffer')
  call reflect_newly_lost_gk_neutrals(sim,passive,reflection,dt,dt,status)
  call assert_equals(0,status,'existing curved-wall reflection commits')
  if(owns(1))then
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      call assert_true(p(1)%i_elm>0,'reflected neutral stays active')
      call assert_equals(original%weight,p(1)%weight,0.d0,'wall-first weight unchanged')
      call assert_true(abs(sum(p(1)%v**2)/sum(original%v**2)-1.d0)<1.d-13,'wall-first reflection energy')
      call assert_true(norm2(p(1)%x(:2)-hit%wall)<2.d-11,'unused time follows existing reflection convention')
    end select
  endif
  do j=1,nodes%n_nodes
    nodes%node(j)%x=saved_x(j,:,:,:)
  enddo
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
end subroutine

subroutine duplicate_face_tests()
  type(particle_sim)::sim
  type(particle_kinetic_leapfrog)::source,endpoint
  type(side_correction)::hit
  integer,allocatable::indices(:)
  integer::e,next,entry,old_node,duplicate,which
  real*8::right_r,right_z,tau
  call fixture(sim,0,indices)
  e=interior;next=elements%element(e)%neighbours(2)
  entry=findloc(elements%element(next)%neighbours,e,dim=1)
  call assert_equals(4,entry,'square mesh reciprocal face')
  old_node=elements%element(next)%vertex(1);duplicate=nodes%n_nodes+1
  call assert_true(duplicate<=size(nodes%node),'test duplicate node capacity')
  nodes%n_nodes=duplicate
  source%i_elm=e;source%st=[.8d0,.4d0]
  call interp_RZ(nodes,elements,e,source%st(1),source%st(2),source%x(1),source%x(2))
  source%x(3)=.3d0;source%v=[10000.d0,0.d0,0.d0]
  call interp_RZ(nodes,elements,e,1.d0,.4d0,right_r,right_z)
  tau=1.5d0*(right_r-source%x(1))/source%v(1)
  elements%element(next)%vertex(1)=duplicate
  do which=1,4
    nodes%node(duplicate)=nodes%node(old_node)
    if(which==2)nodes%node(duplicate)%x(1,1,1)=nodes%node(duplicate)%x(1,1,1)+1.d-8
    if(which==3)nodes%node(duplicate)%x(1,3,1)=nodes%node(duplicate)%x(1,3,1)+1.d-8
    if(which==4)elements%element(next)%neighbours(entry)=0
    endpoint=source;call push_gk_neutral_ballistic(endpoint,tau)
    call verify_neutral_trajectory(nodes,elements,e,source%st,source%x,source%v,tau,hit,endpoint%x,endpoint%v)
    if(which==1)then
      call assert_equals(INTERIOR_ENDPOINT,hit%status,'identical cubic duplicate-node face traverses')
      call assert_equals(next,hit%element,'duplicate-node crossing reaches actual neighbour')
    else
      call assert_equals(UNRESOLVED,hit%status,'displaced/tangent-mismatched/nonreciprocal faces reject')
    endif
  enddo
  elements%element(next)%neighbours(entry)=e
  elements%element(next)%vertex(1)=old_node;nodes%n_nodes=duplicate-1
end subroutine

subroutine selection_transaction_tests()
  type(particle_sim)::sim
  type(particle_snapshot)::saved
  type(particle_slot_registry)::pool_before
  type(gk_adas_ionization_state)::state
  type(gk_neutral_ionization_state)::conversion_before
  type(type_node_list)::profiles
  type(pcg32_rng)::expected_rng
  type(gk_neutral_state)::passive
  type(gk_neutral_reflection_state)::reflection
  integer,allocatable::indices(:),walls(:)
  integer::which,k,j,status,bad,nactive
  real*8::expected(1),actual(1),pressure,dt
  bad=min(1,tasks-1)
  do which=1,4
    call fixture(sim,merge(1,3,which==3),indices)
    if(which==4)then
      do j=1,3
        sim%groups(3)%particles(j)%i_elm=0
      enddo
      indices=[integer::]
    endif
    select type(f=>sim%fields)
    type is(analytic_fields)
      f%vector=[0.d0,0.d0,1.d0]
    end select
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      if(which==2.and.rank==bad)then
        j=indices(1);p(j)%x=[2.5d0-1.d-6,.1d0,0.d0];p(j)%v=[0.d0,1.d5,0.d0]
        call locate_test_particle(p(j))
      endif
      if(which==3.and.owns(1))then
        p(1)%x=[1.13d0,1.d0-1.d-7,0.d0];p(1)%v=[0.d0,10000.d0,0.d0]
        call locate_test_particle(p(1))
      endif
    end select
    call initialize_particle_slot_registry(sim,[1,2,3],registry)
    call snapshot(sim,saved);call snapshot_slot_registry(registry,pool_before)
    profiles=nodes;profiles%n_values=7
    pressure=1.d19*EL_CHG*100.d0
    if(which==1)pressure=1.d19*EL_CHG ! cold guard: draw but never select
    do j=1,profiles%n_nodes
      call init_node(profiles%node(j),7)
      profiles%node(j)%deltas=0.d0
      profiles%node(j)%values=0.d0
      profiles%node(j)%values(1,1,2)=1.d19
      profiles%node(j)%values(1,1,4)=pressure
    enddo
    use_manual_random_seed=.true.;manual_seed=20261103
    call initialize_gk_adas_ionization(sim,state,.true.,1,2,3,'12_h',trim(event_test_data_directory),.false.)
    call cache_ionization_profiles(sim,state,profiles,0.d0,0.d0)
    expected_rng=state%rng;conversion_before=state%conversion
    nactive=count(sim%groups(3)%particles%i_elm>0)
    do j=1,nactive
      call expected_rng%next(expected)
    enddo
    ! A controlled existing RNG override exercises a rare short positive flight.
    use_fixed_rng_value=which==2.or.which==3
    fixed_rng_value=1.d-8;if(which==3)fixed_rng_value=.5d0
    dt=1.d-5
    call advance_gk_adas_ionization(sim,state,registry,1,0.d0,dt,dt,walls,substep_dt=dt/4.d0,ierr=status)
    use_fixed_rng_value=.false.
    if(which==2)then
      call assert_true(status/=0,'inside source/outside GC rejected on all ranks')
      call unchanged(sim,saved,'selection-layer MPI rollback')
      call assert_slot_registry_equal(registry,pool_before,'selection-layer pools rolled back')
      call assert_true(all(state%conversion%counts==conversion_before%counts),'selection-layer diagnostics rolled back')
      call assert_true(all(state%conversion%totals==conversion_before%totals),'selection-layer physical totals rolled back')
      call assert_equals(0,size(walls),'failed selection emits no wall handoff')
      call assert_equals(0.d0,state%elapsed,0.d0,'failed step has no committed timing')
      call expected_rng%initialize(1,20261103,tasks,100000000+rank)
    else
      call assert_equals(0,status,'zero-event or wall-first step commits')
      call unchanged(sim,saved,'no ionizing event leaves all source/product states unchanged')
      call assert_slot_registry_equal(registry,pool_before,'no event leaves pools unchanged')
      if(which==3)then
        call assert_equals(size(indices),size(walls),'all sampled events are wall-first')
        call initialize_gk_neutrals(sim,3,passive)
        call initialize_gk_neutral_reflection(sim,3,.true.,reflection)
        call advance_gk_neutrals(sim,passive,dt/4.d0,4,verify_indices=walls)
        call reflect_newly_lost_gk_neutrals(sim,passive,reflection,dt,dt/4.d0,status)
        call assert_equals(0,status,'wall-first selection uses normal reflection transaction')
        if(owns(1))call assert_true(sim%groups(3)%particles(1)%i_elm>0,'wall-first source reflected')
      else
        call assert_equals(0,size(walls),'no wall candidates without Poisson events')
      endif
    endif
    call expected_rng%next(expected);call state%rng%next(actual)
    call assert_true(all(expected==actual),'exact one-draw consumption or global RNG rollback')
    close(state%event_unit)
    if(rank==0)close(state%step_unit)
  enddo
end subroutine

subroutine copy_test_nodes(source,destination)
  type(type_node_list),intent(in)::source
  type(type_node_list),intent(out)::destination
  destination=source
end subroutine

subroutine gc_boundary_selection_tests()
  type(particle_sim)::sim
  type(particle_snapshot)::saved
  type(particle_slot_registry)::pool_before
  type(gk_adas_ionization_state)::state
  type(type_node_list)::profiles
  type(pcg32_rng)::expected_rng
  type(gk_neutral_state)::passive
  type(gk_neutral_reflection_state)::reflection
  type(particle_kinetic_leapfrog)::endpoint
  integer,allocatable::indices(:),walls(:)
  integer::which,j,k,status,nsource,nvalid,noutside,broken,old_neighbor,e,mode,steps
  real*8::expected(1),actual(1),dt,r,z,wvalid,woutside
  logical::fatal,keep
  character(len=16)::action
  do which=1,8
    ! 1 strict outside; 2 keep outside; 3 mixed; 4 mixed numeric error;
    ! 5/6 numeric-only strict/keep; 7 bad capacity with outside; 8 wall after veto.
    nsource=1;if(which==3.or.which==4.or.which==7)nsource=3
    call fixture(sim,nsource,indices)
    keep=which/=1.and.which/=5;action='abort';if(keep)action='keep_neutral'
    fatal=which==1.or.which==4.or.which==5.or.which==6.or.which==7
    select type(f=>sim%fields)
    type is(analytic_fields)
      f%vector=[0.d0,0.d0,1.d0]
    end select
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      j=1;if(nsource==3)j=2
      if(owns(j))then
        p(j)%x=[2.5d0-1.d-6,.1d0,0.d0];p(j)%v=[0.d0,1.d5,0.d0]
        if(which==8)p(j)%x(:2)=[2.499d0,.99d0]
        call locate_test_particle(p(j))
      endif
      if(which==4.or.which==5.or.which==6)then
        j=1;if(which==4)j=3
        if(owns(j))then
          ! GC crosses an INTERNAL face. Damaged neighbour cannot authorize a veto.
          call interp_RZ(nodes,elements,interior,1.d0,.5d0,r,z)
          p(j)%x=[r-1.d-6,z,0.d0];p(j)%v=[0.d0,1.d5,0.d0]
          call locate_test_particle(p(j))
        endif
      endif
    end select
    call initialize_particle_slot_registry(sim,[1,2,3],registry)
    if(which==7.and.owns(1))then
      ! Exhaust local product capacity with legitimate occupied slots.
      sim%groups(1)%particles(3:5)%i_elm=interior
      call initialize_particle_slot_registry(sim,[1,2,3],registry)
    endif
    call snapshot(sim,saved);call snapshot_slot_registry(registry,pool_before)
    call copy_test_nodes(nodes,profiles);profiles%n_values=7
    do j=1,profiles%n_nodes
      call init_node(profiles%node(j),7);profiles%node(j)%deltas=0.d0;profiles%node(j)%values=0.d0
      profiles%node(j)%values(1,1,2)=1.d19;profiles%node(j)%values(1,1,4)=1.d19*EL_CHG*100.d0
    enddo
    use_manual_random_seed=.true.;manual_seed=20261103
    call initialize_gk_adas_ionization(sim,state,.true.,1,2,3,'12_h',trim(event_test_data_directory), &
      .false.,action,.true.)
    call assert_true(state%boundary%valid,'closed square wall certified')
    call cache_ionization_profiles(sim,state,profiles,0.d0,0.d0)
    broken=0;old_neighbor=0
    if(which==4.or.which==5.or.which==6)then
      j=1;if(which==4)j=3
      if(owns(j))then
        broken=interior;old_neighbor=elements%element(broken)%neighbours(2)
        elements%element(broken)%neighbours(2)=0
      endif
    endif
    expected_rng=state%rng
    if(.not.fatal)then
      do j=1,size(indices)
        call expected_rng%next(expected)
      enddo
    endif
    use_fixed_rng_value=.true.;fixed_rng_value=1.d-8;dt=1.d-8
    if(which==8)dt=1.d-6
    call advance_gk_adas_ionization(sim,state,registry,1,0.d0,dt,dt,walls,substep_dt=dt,ierr=status)
    use_fixed_rng_value=.false.
    if(broken>0)elements%element(broken)%neighbours(2)=old_neighbor
    if(fatal)then
      call assert_true(status/=0,'strict physical veto, numeric defect or capacity failure stays global error')
      call unchanged(sim,saved,'fatal mixed step rolls back all particles')
      call assert_slot_registry_equal(registry,pool_before,'fatal mixed step rolls back pools')
      call assert_true(all(state%conversion%totals==0.d0),'fatal step has no binding/energy source')
      call assert_true(all(state%boundary_counts==0),'fatal step has no committed veto counts')
      call assert_true(all(state%boundary_weights==0.d0),'fatal step has no committed veto weight')
      call assert_equals(0,size(walls),'fatal step emits no passive handoff')
    else
      call assert_equals(0,status,'physical veto is nonfatal')
      nvalid=0;wvalid=0.d0;noutside=1;woutside=1.d10
      if(which==3)then
        nvalid=2;wvalid=4.d10;woutside=2.d10
      endif
      call assert_true(all(state%boundary_counts==[nsource,0,nvalid,noutside,0]),'global filtered counts with empty ranks')
      call assert_true(all(state%boundary_weights==[wvalid+woutside,0.d0,wvalid,woutside,0.d0]),'global physical weights')
      do j=1,nsource
        if(.not.owns(j))cycle
        if(which==3.and.j/=2)then
          call assert_equals(0,sim%groups(3)%particles(j)%i_elm,'valid mixed source converted')
        else
          select type(p=>sim%groups(3)%particles)
          type is(particle_kinetic_leapfrog)
            call assert_true(same_neutral(p(j),saved%neutrals(j)),'veto leaves exact original neutral before passive advance')
          end select
          call assert_true(any(walls==j),'veto receives forced curved-wall verification')
        endif
      enddo
      if(which/=3)then
        call unchanged(sim,saved,'veto creates no products')
        call assert_slot_registry_equal(registry,pool_before,'veto consumes/releases no slots')
        call assert_true(all(state%conversion%totals==0.d0),'veto adds no energy or binding obligation')
        call initialize_gk_neutrals(sim,3,passive)
        call initialize_gk_neutral_reflection(sim,3,.true.,reflection)
        call advance_gk_neutrals(sim,passive,dt,1,verify_indices=walls)
        if(owns(1))then
          if(which==8)then
            call assert_equals(1,passive%n_newly_lost_indices,'wall AFTER veto still detected')
          else
            endpoint=saved%neutrals(1);call push_gk_neutral_ballistic(endpoint,dt)
            call assert_true(all(sim%groups(3)%particles(1)%x==endpoint%x),'veto takes full outer flight from original')
          endif
        endif
        call reflect_newly_lost_gk_neutrals(sim,passive,reflection,dt,dt,status)
        call assert_equals(0,status,'normal post-veto reflection succeeds')
        if(owns(1))then
          select type(p=>sim%groups(3)%particles)
          type is(particle_kinetic_leapfrog)
            call assert_true(p(1)%i_elm>0.and.p(1)%q==0,'vetoed neutral stays active after transport')
            call assert_equals(saved%neutrals(1)%weight,p(1)%weight,0.d0,'post-veto transport preserves weight')
            call assert_true(abs(sum(p(1)%v**2)/sum(saved%neutrals(1)%v**2)-1.d0)<1.d-13,'post-veto energy')
          end select
        endif
        if(which==2)then
          ! Same retained neutral can be attempted again; no residual-time resampling.
          do j=1,size(indices)
            call expected_rng%next(expected)
          enddo
          use_fixed_rng_value=.true.;fixed_rng_value=1.d-8
          call advance_gk_adas_ionization(sim,state,registry,2,dt,dt,2.d0*dt,walls,ierr=status)
          use_fixed_rng_value=.false.
          call assert_equals(0,status,'later outer timestep can retry retained neutral')
          call assert_equals(2,int(state%boundary_counts(4)),'two rejected attempts, one physical neutral')
        endif
      endif
    endif
    call expected_rng%next(expected);call state%rng%next(actual)
    call assert_true(all(expected==actual),'exact one draw per active marker or global fatal rollback')
    close(state%event_unit);close(state%rejection_unit)
    if(rank==0)then
      close(state%step_unit);close(state%boundary_unit)
    endif
  enddo
end subroutine

subroutine curved_veto_transport()
  type(particle_sim)::sim
  type(gk_adas_ionization_state)::state
  type(type_node_list)::profiles
  type(pcg32_rng)::expected_rng
  real*8::expected(1),actual(1)
  type(gk_neutral_state)::passive
  type(gk_neutral_reflection_state)::reflection
  type(particle_kinetic_leapfrog)::original,source,endpoint
  type(side_correction)::hit
  integer,allocatable::indices(:),forced(:)
  real*8,allocatable::saved_x(:,:,:,:)
  real*8::z,wall,dt,alpha,energy
  integer::j,k,status
  allocate(saved_x(nodes%n_nodes,n_coord_tor,n_degrees,2));alpha=.2d0
  do j=1,nodes%n_nodes
    saved_x(j,:,:,:)=nodes%node(j)%x
    ! R=Rbase+alpha*Z^2, with exact chain-rule nodal derivatives.
    z=nodes%node(j)%x(1,1,2)
    nodes%node(j)%x(1,1,1)=nodes%node(j)%x(1,1,1)+alpha*z*z
    do k=2,3
      nodes%node(j)%x(1,k,1)=nodes%node(j)%x(1,k,1)+2.d0*alpha*z*nodes%node(j)%x(1,k,2)
    enddo
    nodes%node(j)%x(1,4,1)=nodes%node(j)%x(1,4,1)+2.d0*alpha* &
      (nodes%node(j)%x(1,2,2)*nodes%node(j)%x(1,3,2)+z*nodes%node(j)%x(1,4,2))
  enddo
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
  call fixture(sim,1,indices)
  original%x=[2.52d0,-.4d0,0.d0];original%v=[0.d0,8000.d0,0.d0]
  original%q=0;original%weight=1.d10;original%i_life=4;original%t_birth=-1.e0
  call locate_test_particle(original)
  dt=1.d-4;endpoint=original;call push_gk_neutral_ballistic(endpoint,dt)
  call locate_test_particle(endpoint)
  call assert_true(endpoint%i_elm>0,'curved-exit test endpoint is inside')
  call stage_ionization_motion(sim,original,dt,dt,source,hit,wall)
  call assert_equals(WALL_CROSSING,hit%status,'first curved wall crossing despite inside endpoint')
  call assert_equals((.4d0-sqrt(.1d0))/8000.d0,wall,2.d-13,'analytic curved-wall first contact')
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    if(owns(1))p(1)=original
  end select
  select type(f=>sim%fields)
  type is(analytic_fields)
    f%vector=[0.d0,0.d0,.005d0]
  end select
  call copy_test_nodes(nodes,profiles);profiles%n_values=7
  do j=1,profiles%n_nodes
    call init_node(profiles%node(j),7);profiles%node(j)%deltas=0.d0;profiles%node(j)%values=0.d0
    profiles%node(j)%values(1,1,2)=1.d19;profiles%node(j)%values(1,1,4)=1.d19*EL_CHG*100.d0
  enddo
  use_manual_random_seed=.true.;manual_seed=20261103
  call initialize_gk_adas_ionization(sim,state,.true.,1,2,3,'12_h',trim(event_test_data_directory), &
    .false.,'keep_neutral',.true.)
  call assert_true(state%boundary%valid,'curved closed boundary certificate')
  call cache_ionization_profiles(sim,state,profiles,0.d0,0.d0)
  expected_rng=state%rng
  if(owns(1))call expected_rng%next(expected)
  use_fixed_rng_value=.true.;fixed_rng_value=1.d-8
  call advance_gk_adas_ionization(sim,state,registry,1,0.d0,dt,dt,forced,ierr=status)
  use_fixed_rng_value=.false.
  call assert_equals(0,status,'curved physical outside GC is nonfatal')
  call assert_true(all(state%boundary_counts==[1,0,0,1,0]),'curved candidate vetoed BEFORE first wall')
  call assert_equals(size(indices),size(forced),'veto is handed to full-interval verifier')
  if(owns(1))then
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      call assert_true(same_neutral(p(1),original),'curved veto retains original source')
    end select
  endif
  call initialize_gk_neutrals(sim,3,passive)
  use_manual_random_seed=.true.;manual_seed=20261103
  call initialize_gk_neutral_reflection(sim,3,.true.,reflection)
  call advance_gk_neutrals(sim,passive,dt,1,verify_indices=forced)
  if(owns(1))call assert_equals(1,passive%n_newly_lost_indices,'wall-first candidate reaches existing loss buffer')
  call reflect_newly_lost_gk_neutrals(sim,passive,reflection,dt,dt,status)
  call assert_equals(0,status,'existing curved-wall reflection commits')
  if(owns(1))then
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      call assert_true(p(1)%i_elm>0,'reflected neutral stays active')
      call assert_equals(original%weight,p(1)%weight,0.d0,'wall-first weight unchanged')
      call assert_true(abs(sum(p(1)%v**2)/sum(original%v**2)-1.d0)<1.d-13,'wall-first reflection energy')
      call assert_true(norm2(p(1)%x(:2)-hit%wall)<2.d-11,'unused time follows existing reflection convention')
    end select
  endif
  call expected_rng%next(expected);call state%rng%next(actual)
  call assert_true(all(expected==actual),'post-veto reflection does not draw ionization stream')
  close(state%event_unit);close(state%rejection_unit)
  if(rank==0)then
    close(state%step_unit);close(state%boundary_unit)
  endif
  do j=1,nodes%n_nodes
    nodes%node(j)%x=saved_x(j,:,:,:)
  enddo
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
end subroutine

end module
