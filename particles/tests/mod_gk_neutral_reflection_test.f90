module mod_gk_neutral_reflection_test
use mpi
#ifdef GK_NEUTRAL_REFINEMENT_TEST
! The isolated test build exposes these private names without changing production visibility.
use mod_gk_neutral_trajectory, only: trajectory_origin,refine_side_intersection,side_geometry
#endif
use fruit
use data_structure
use, intrinsic :: ieee_arithmetic
use mod_particle_sim, only: particle_sim,configure_particle_groups
use mod_particle_types
use mod_gk_neutrals
use mod_gk_neutral_trajectory, only: verify_neutral_trajectory,side_correction, &
  WALL_CROSSING,INTERIOR_ENDPOINT,UNRESOLVED,max_trace_elements
use mod_gk_neutral_reflection
use mod_particle_slot_pool
use mod_particle_slot_pool_test, only: snapshot_slot_registry,assert_slot_registry_equal
use mod_fields_linear, only: jorek_fields_interp_linear
use mod_projection_helpers_test_tools, only: default_square_grid
use mod_neighbours, only: update_neighbours
use mod_boundary, only: wall_normal_vector
use mod_interp, only: interp_RZ
use mod_pcg32_rng, only: pcg32_rng
use mod_sampling, only: sample_cosine
use constants, only: ATOMIC_MASS_UNIT
implicit none
private
public :: run_reflection_tests,run_reflection_invalid
integer :: rank,tasks,wall_element,interior_element,corner_element
real*8 :: wall_st(2),corner_st(2),nan_value
real*8, parameter :: dt=1.d-4
integer, parameter :: nslots=6
integer, parameter :: n_invalid=22
character(len=20), parameter :: invalid_names(n_invalid)=[character(len=20):: &
 'index','duplicate','active','free','charge','weight_zero','weight_negative','weight_nan', &
 'position_nan','velocity_nan','st_nan','st_range','element','interior','opposite', &
 'corner','energy','inward','mass','mass_nan','event_count','staging_overflow']
type(type_node_list), target :: nodes
type(type_element_list), target :: elements
contains
subroutine mesh(nradial)
  integer, optional, intent(in) :: nradial
  integer :: error,e,s,n1,n2,k,nr
  logical :: sides(4)
  real*8 :: r,z
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,tasks,error)
  nan_value=ieee_value(0.d0,ieee_quiet_nan)
  call init_node_list(nodes,n_nodes_max,nodes%n_dof,n_var)
  nr=5;if(present(nradial)) nr=nradial
  call default_square_grid(rank,tasks,nr,5,nodes,elements,error,0.5d0,2.5d0,-1.d0,1.d0)
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
  wall_element=0
  do e=1,elements%n_elements
    do s=1,4
      n1=elements%element(e)%vertex(s); n2=elements%element(e)%vertex(mod(s,4)+1)
      sides(s)=nodes%node(n1)%boundary>0.and.nodes%node(n2)%boundary>0
    enddo
    if(count(sides)==0) interior_element=e
    if(count(sides)==1.and.sides(2)) then
      call interp_RZ(nodes,elements,e,1.d0,0.4d0,r,z)
      if(abs(r-2.5d0)<1.d-12) then
        wall_element=e; wall_st=[1.d0,0.4d0]
      endif
    endif
    if(count(sides)==2) then
      corner_element=e
      do k=1,4
        if(.not.(sides(k).and.sides(mod(k+2,4)+1))) cycle
        select case(k)
        case(1); corner_st=[0.d0,0.d0]
        case(2); corner_st=[1.d0,0.d0]
        case(3); corner_st=[1.d0,1.d0]
        case(4); corner_st=[0.d0,1.d0]
        end select
      enddo
    endif
  enddo
  if(wall_element==0) call MPI_ABORT(MPI_COMM_WORLD,2,error)
end subroutine
logical function owns(j)
  integer, intent(in) :: j
  owns=rank==mod(j-1,tasks)
end function
subroutine fixture(sim,passive,state,number,enabled)
  use phys_module, only: use_manual_random_seed,manual_seed
  type(particle_sim), intent(out) :: sim
  type(gk_neutral_state), intent(out) :: passive
  type(gk_neutral_reflection_state), intent(out) :: state
  integer, intent(in) :: number
  logical, optional, intent(in) :: enabled
  integer :: j
  logical :: on
  real*8 :: r,z
  on=.true.; if(present(enabled)) on=enabled
  use_manual_random_seed=.true.; manual_seed=12345
  sim%my_id=rank; sim%n_mpi=tasks; sim%time=0.125d0
  allocate(jorek_fields_interp_linear::sim%fields)
  sim%fields%node_list=>nodes; sim%fields%element_list=>elements
  allocate(sim%groups(3))
  sim%groups%mass=2.01410178d0; sim%groups%coupling_scheme='non'; sim%groups%ad%suffix=''
  sim%groups(1)%id='ion'; sim%groups(2)%id='ele'; sim%groups(3)%id='neu'
  allocate(particle_gc_vpar::sim%groups(1)%particles(2),sim%groups(2)%particles(2))
  do j=1,2
    sim%groups(j)%particles%i_elm=interior_element
    sim%groups(j)%particles%weight=123.d0+j
    select type(p=>sim%groups(j)%particles)
    type is(particle_gc_vpar)
      p%q=3-2*j
    end select
  enddo
  allocate(particle_kinetic_leapfrog::sim%groups(3)%particles(nslots))
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,number
      if(.not.owns(j)) cycle
      p(j)%i_elm=wall_element; p(j)%st=[0.999d0,0.2d0+0.15d0*j]
      call interp_RZ(nodes,elements,wall_element,p(j)%st(1),p(j)%st(2),r,z)
      p(j)%x=[r,z,7.d0+0.2d0*j]; p(j)%v=[100.d0*j,5.d0*j,7.d0*j]
      p(j)%q=0; p(j)%weight=1.d10*j; p(j)%i_life=20+j; p(j)%t_birth=0.0625_4
    enddo
    p(4)%i_elm=-wall_element; p(4)%x=[2.51d0,0.1d0,0.7d0]; p(4)%v=[2.d0,0.d0,0.d0]
    p(4)%weight=3.d13; p(4)%st=wall_st
    p(6)%i_elm=interior_element; p(6)%st=0.5d0
    call interp_RZ(nodes,elements,interior_element,0.5d0,0.5d0,r,z)
    p(6)%x=[r,z,0.d0]; p(6)%v=0.d0; p(6)%weight=7.d0
  end select
  call initialize_gk_neutrals(sim,3,passive)
  call initialize_gk_neutral_reflection(sim,3,on,state)
end subroutine
logical function bits(a,b)
  real*8, intent(in) :: a(:),b(:)
  bits=all(transfer(a,[0_8],size(a))==transfer(b,[0_8],size(b)))
end function
subroutine same_particles(a,b,label)
  class(particle_base), intent(in) :: a(:),b(:)
  character(len=*), intent(in) :: label
  integer :: j
  call assert_equals(size(a),size(b),label//' size')
  if(size(a)>0) call assert_equals(codify_particle_type(a(1)),codify_particle_type(b(1)),label//' type')
  do j=1,size(a)
    call assert_true(bits(a(j)%x,b(j)%x).and.bits(a(j)%st,b(j)%st),label//' position bits')
    call assert_true(bits([a(j)%weight],[b(j)%weight]),label//' weight bits')
    call assert_true(a(j)%i_elm==b(j)%i_elm.and.particle_get_q(a(j))==particle_get_q(b(j)).and.a(j)%i_life==b(j)%i_life.and. &
                     transfer(a(j)%t_birth,0)==transfer(b(j)%t_birth,0),label//' integer/birth bits')
  enddo
  select type(a)
  type is(particle_kinetic_leapfrog)
    select type(b)
    type is(particle_kinetic_leapfrog)
      do j=1,size(a)
        call assert_true(bits(a(j)%v,b(j)%v),label//' velocity bits')
      enddo
    end select
  type is(particle_gc_vpar)
    select type(b)
    type is(particle_gc_vpar)
      do j=1,size(a)
        call assert_true(bits([a(j)%vpar,a(j)%mu,a(j)%B_norm],[b(j)%vpar,b(j)%mu,b(j)%B_norm]),label//' GC bits')
      enddo
    end select
  end select
end subroutine
subroutine same_state(a,b,label)
  type(gk_neutral_reflection_state), intent(in) :: a,b
  character(len=*), intent(in) :: label
  call assert_true(all(a%counts==b%counts).and.all(a%global_counts==b%global_counts),label//' counts')
  call assert_true(a%trajectory_count==b%trajectory_count.and. &
                   a%global_trajectory_count==b%global_trajectory_count,label//' trajectory counts')
  call assert_true(bits(a%totals,b%totals).and.bits(a%global_totals,b%global_totals),label//' totals')
  call assert_true(all(transfer(a%rng,[0_1],storage_size(a%rng)/8)== &
                       transfer(b%rng,[0_1],storage_size(b%rng)/8)),label//' RNG bits')
  call assert_true((a%initialized.eqv.b%initialized).and.(a%enabled.eqv.b%enabled).and. &
                   a%neutral_group==b%neutral_group,label//' flags')
end subroutine
! Snapshot every particle group, pool array, RNG, local/global diagnostics and events.
subroutine checked_reflect(sim,passive,state,error,expect_failure)
  type(particle_sim), intent(inout) :: sim
  type(gk_neutral_state), intent(inout) :: passive
  type(gk_neutral_reflection_state), intent(inout) :: state
  integer, intent(out) :: error
  logical, intent(in) :: expect_failure
  class(particle_base), allocatable :: ions(:),electrons(:),neutrals(:)
  type(particle_slot_registry) :: pool,before
  type(gk_neutral_reflection_state) :: old
  integer, allocatable :: events(:)
  type(gk_neutral_event), allocatable :: records(:)
  integer :: count_before,valid,j
  allocate(ions,source=sim%groups(1)%particles); allocate(electrons,source=sim%groups(2)%particles)
  allocate(neutrals,source=sim%groups(3)%particles)
  call initialize_particle_slot_registry(sim,[1,2,3],pool)
  call snapshot_slot_registry(pool,before)
  old=state; events=passive%newly_lost_indices; count_before=passive%n_newly_lost_indices
  records=passive%newly_lost_events
  call reflect_newly_lost_gk_neutrals(sim,passive,state,sim%time,dt,error)
  call assert_true((error/=0).eqv.expect_failure,'expected reflection status')
  call assert_slot_registry_equal(pool,before,'reflection never touches pool')
  call validate_particle_slot_registry(sim,pool,[1,2,3],valid)
  call assert_equals(0,valid,'pool remains physically valid')
  call same_particles(sim%groups(1)%particles,ions,'ions protected')
  call same_particles(sim%groups(2)%particles,electrons,'electrons protected')
  if(error/=0) then
    call same_particles(sim%groups(3)%particles,neutrals,'global rollback')
    call same_state(state,old,'global rollback')
    call assert_equals(count_before,passive%n_newly_lost_indices,'rollback event count')
    call assert_true(all(events==passive%newly_lost_indices),'rollback entire event buffer')
    call assert_equals(size(records),size(passive%newly_lost_events),'rollback record capacity')
    do j=1,size(records)
      call assert_true(records(j)%marker==passive%newly_lost_events(j)%marker.and. &
        records(j)%element==passive%newly_lost_events(j)%element.and. &
        records(j)%substep==passive%newly_lost_events(j)%substep,'rollback event integers')
      call assert_true(bits(records(j)%x_old,passive%newly_lost_events(j)%x_old).and. &
        bits(records(j)%v_old,passive%newly_lost_events(j)%v_old).and. &
        bits(records(j)%st_old,passive%newly_lost_events(j)%st_old).and. &
        bits([records(j)%dt],[passive%newly_lost_events(j)%dt]),'rollback event provenance bits')
    enddo
  endif
end subroutine
subroutine noops_and_absorption()
  type(particle_sim) :: sim,empty
  type(gk_neutral_state) :: passive,absent
  type(gk_neutral_reflection_state) :: state,old
  class(particle_base), allocatable :: saved(:)
  integer :: error
  call fixture(sim,passive,state,1,.false.)
  call advance_gk_neutrals(sim,passive,dt,3)
  allocate(saved,source=sim%groups(3)%particles); old=state
  call checked_reflect(sim,passive,state,error,.false.)
  call same_particles(saved,sim%groups(3)%particles,'disabled reflection')
  call same_state(state,old,'disabled state')
  call advance_gk_neutrals(sim,passive,dt,3)
  call same_particles(saved,sim%groups(3)%particles,'absorbed marker frozen')
  call assert_equals(0,passive%n_newly_lost_indices,'old negative never re-enqueued')
  if(owns(1)) then
    call assert_true(saved(1)%i_elm<0.and.saved(1)%x(1)>2.5d0,'disabled overshot position remains')
    select type(p=>saved)
    type is(particle_kinetic_leapfrog)
      call assert_true(p(1)%v(1)>0.d0,'disabled outward velocity remains')
    end select
  endif
  empty%my_id=rank; empty%n_mpi=tasks; allocate(empty%groups(2))
  call initialize_gk_neutrals(empty,3,absent)
  call initialize_gk_neutral_reflection(empty,3,.true.,state,error); old=state
  call reflect_newly_lost_gk_neutrals(empty,absent,state,0.d0,dt,error)
  call report_gk_neutral_reflection(empty,state)
  call same_state(state,old,'absent group no-op')
  call fixture(sim,passive,state,0)
  deallocate(sim%groups(3)%particles)
  allocate(particle_kinetic_leapfrog::sim%groups(3)%particles(0))
  call advance_gk_neutrals(sim,passive,dt,1)
  call assert_equals(0,size(passive%newly_lost_indices),'event buffer resizes to empty array')
  call assert_equals(0,size(passive%newly_lost_events),'metadata buffer resizes to empty array')
  old=state
  call checked_reflect(sim,passive,state,error,.false.)
  call same_state(state,old,'empty events no-op')
end subroutine
subroutine reflections()
  type(particle_sim) :: sim,replay
  type(gk_neutral_state) :: passive,again
  type(gk_neutral_reflection_state) :: state,other,before
  type(particle_kinetic_leapfrog) :: pre
  type(pcg32_rng) :: expected_rng
  class(particle_base), allocatable :: hit(:),saved(:)
  real*8 :: r,z,normal(3),u(2),expected(3),expected_energy
  integer :: j,error,n,k
  call fixture(sim,passive,state,3)
  call fixture(replay,again,other,3)
  call advance_gk_neutrals(sim,passive,dt,3)
  call advance_gk_neutrals(replay,again,dt,3)
  allocate(hit,source=sim%groups(3)%particles)
  n=count([(owns(j),j=1,3)])
  call assert_equals(n,passive%n_newly_lost_indices,'one event per physical current loss')
  call assert_true(all(passive%newly_lost_indices(:n)<=3),'free/old-negative/interior protected')
  do k=1,n
    j=passive%newly_lost_indices(k)
    pre%x=passive%newly_lost_events(k)%x_old;pre%v=passive%newly_lost_events(k)%v_old
    call push_gk_neutral_ballistic(pre,passive%newly_lost_events(k)%dt)
    select type(h=>hit)
    type is(particle_kinetic_leapfrog)
      call assert_true(bits(pre%x,h(j)%x).and.bits(pre%v,h(j)%v),'producer provenance reproduces full endpoint')
    end select
    call assert_equals(j,passive%newly_lost_events(k)%marker,'producer provenance marker')
    call assert_true(passive%newly_lost_events(k)%element>0.and. &
                     passive%newly_lost_events(k)%substep>=1,'producer positive pre-element and substep')
  enddo
  expected_rng=state%rng
  call checked_reflect(sim,passive,state,error,.false.)
  call checked_reflect(replay,again,other,error,.false.)
  call same_particles(sim%groups(3)%particles,replay%groups(3)%particles,'fixed-seed replay')
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
  select type(h=>hit)
  type is(particle_kinetic_leapfrog)
    do j=1,3
      if(.not.owns(j)) cycle
      call interp_RZ(nodes,elements,-h(j)%i_elm,h(j)%st(1),h(j)%st(2),r,z)
      normal=wall_normal_vector(nodes,elements,-h(j)%i_elm,h(j)%st(1),h(j)%st(2))
      call expected_rng%next(u); expected=sqrt(sum(h(j)%v**2))*sample_cosine(u,normal)
      call assert_true(bits(p(j)%v,expected),'exact cosine sample; two uniforms')
      call assert_true(bits(p(j)%x,[r,z,h(j)%x(3)]),'exact wall RZ; retained continuous end phi')
      call assert_true(bits(p(j)%st,h(j)%st),'retained st unchanged')
      call assert_equals(-h(j)%i_elm,p(j)%i_elm,'positive last element')
      call assert_true(p(j)%weight==h(j)%weight.and.p(j)%q==h(j)%q.and.p(j)%i_life==h(j)%i_life.and. &
                       p(j)%t_birth==h(j)%t_birth,'marker identity preserved')
      call assert_true(dot_product(p(j)%v,normal)>0.d0,'outgoing inward')
      call assert_true(abs(sum(p(j)%v**2)/sum(h(j)%v**2)-1.d0)<128*epsilon(1.d0),'speed squared conserved')
    enddo
  end select
  end select
  call same_particles(sim%groups(3)%particles(4:6),hit(4:6),'unrelated slots protected')
  call assert_equals(0,passive%n_newly_lost_indices,'success consumes current event list')
  allocate(saved,source=sim%groups(3)%particles); before=state
  call checked_reflect(sim,passive,state,error,.false.)
  call same_particles(sim%groups(3)%particles,saved,'second reflection no-op')
  call same_state(state,before,'second reflection no-op')
  call report_gk_neutral_reflection(sim,state); before=state
  call report_gk_neutral_reflection(sim,state)
  call same_state(state,before,'repeated report no double reduction')
  call same_particles(sim%groups(3)%particles,saved,'report immutable particles')
  call assert_equals(3,state%global_counts(1),'MPI-independent physical count')
  call assert_equals(3,state%global_counts(2),'all incident markers reflected')
  call assert_true(all(state%global_totals(1:2)==6.d10),'MPI-independent weight')
  expected_energy=0.5d0*2.01410178d0*ATOMIC_MASS_UNIT*1.d10*10074.d0*36.d0
  call assert_true(maxval(abs(state%global_totals(3:4)/expected_energy-1.d0))<128*epsilon(1.d0),'global energy')
  call restart_roundtrip(sim)
  call advance_gk_neutrals(sim,passive,dt,1)
  do j=1,3
    if(.not.owns(j)) cycle
    call assert_true(sim%groups(3)%particles(j)%i_elm>0,'exact-wall start moves inward next step')
    call assert_true(sim%groups(3)%particles(j)%x(1)<saved(j)%x(1),'moves inward without displacement')
  enddo
  call assert_equals(0,passive%n_newly_lost_indices,'no immediate new wall hit')
  ! Arrange a second exit with the same speed and identity; no new marker/slot.
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,3
      if(.not.owns(j)) cycle
      p(j)%i_elm=wall_element; p(j)%st=[0.999d0,0.4d0]
      call interp_RZ(nodes,elements,wall_element,p(j)%st(1),p(j)%st(2),r,z)
      p(j)%x(1:2)=[r,z]; p(j)%v=[sqrt(sum(p(j)%v**2)),0.d0,0.d0]
    enddo
  end select
  call advance_gk_neutrals(sim,passive,dt,1)
  call checked_reflect(sim,passive,state,error,.false.)
  call report_gk_neutral_reflection(sim,state)
  call assert_equals(6,state%global_counts(2),'second physical impact counted once')
  call assert_true(all(state%global_totals(1:2)==12.d10),'repeated weight conserved')
end subroutine
subroutine grazing_and_zero()
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state,before
  integer :: error,j
  call fixture(sim,passive,state,1)
  call advance_gk_neutrals(sim,passive,dt,1); before=state
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    if(owns(1)) p(1)%v=0.d0
  end select
  call checked_reflect(sim,passive,state,error,.false.)
  call assert_true(all(transfer(before%rng,[0_1],storage_size(before%rng)/8)== &
                       transfer(state%rng,[0_1],storage_size(state%rng)/8)),'zero speed consumes no RNG')
  call assert_true(all(state%totals(3:4)==0.d0),'zero speed creates no energy')
  call fixture(sim,passive,state,1)
  call advance_gk_neutrals(sim,passive,dt,1)
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    if(owns(1)) p(1)%v=[-1.d-14,1.d0,0.d0]
  end select
  call checked_reflect(sim,passive,state,error,.false.)
  call report_gk_neutral_reflection(sim,state)
  call assert_equals(1,state%global_counts(3),'roundoff inward grazing accepted and counted')
end subroutine
subroutine corrupt(sim,passive,state,name)
  type(particle_sim), intent(inout) :: sim
  type(gk_neutral_state), intent(inout) :: passive
  type(gk_neutral_reflection_state), intent(inout) :: state
  character(len=*), intent(in) :: name
  integer :: j,k,e
  if(rank/=min(1,tasks-1)) return
  j=rank+1
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    select case(name)
    case('index'); passive%newly_lost_indices(1)=nslots+1
    case('duplicate'); passive%n_newly_lost_indices=2; passive%newly_lost_indices(1:2)=j
    case('active'); p(j)%i_elm=wall_element
    case('free'); p(j)%i_elm=0
    case('charge'); p(j)%q=1
    case('weight_zero'); p(j)%weight=0.d0
    case('weight_negative'); p(j)%weight=-1.d0
    case('weight_nan'); p(j)%weight=nan_value
    case('position_nan'); p(j)%x(1)=nan_value
    case('velocity_nan'); p(j)%v(1)=nan_value
    case('st_nan'); p(j)%st(1)=nan_value
    case('st_range'); p(j)%st(1)=1.1d0
    case('element'); p(j)%i_elm=-huge(0)-1
    case('interior'); p(j)%i_elm=-interior_element
    case('opposite')
      e=interior_element; p(j)%i_elm=-e
      do k=1,4
        nodes%node(elements%element(e)%vertex(k))%boundary=mod(k,2)
      enddo
    case('corner'); p(j)%i_elm=-corner_element; p(j)%st=corner_st
    case('energy'); p(j)%v=huge(1.d0)/2.d0
    case('inward'); p(j)%v=-p(j)%v
    case('mass'); sim%groups(3)%mass=0.d0
    case('mass_nan'); sim%groups(3)%mass=nan_value
    case('event_count'); passive%n_newly_lost_indices=nslots+1
    case('staging_overflow'); state%totals(1)=huge(1.d0); p(j)%weight=huge(1.d0)/2.d0
    end select
  end select
end subroutine
subroutine rejections()
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state
  integer :: k,error
  integer, allocatable :: boundaries(:)
  boundaries=nodes%node%boundary
  do k=1,n_invalid
    call fixture(sim,passive,state,3)
    call advance_gk_neutrals(sim,passive,dt,1)
    call corrupt(sim,passive,state,trim(invalid_names(k)))
    call checked_reflect(sim,passive,state,error,.true.)
    nodes%node%boundary=boundaries
  enddo
  call fixture(sim,passive,state,0)
  deallocate(sim%groups(3)%particles)
  allocate(particle_gc_vpar::sim%groups(3)%particles(nslots))
  call initialize_gk_neutral_reflection(sim,3,.true.,state,error)
  call assert_true(error/=0,'wrong dynamic neutral type rejected collectively')
end subroutine
subroutine restart_roundtrip(sim)
  use hdf5, only: H5open_f
  use mod_particle_io, only: write_simulation_hdf5,read_simulation_hdf5
  use phys_module, only: n_part_groups,part_groups_in_use,part_group_configs
  use mod_particle_group_id, only: matching_part_config_indices
  type(particle_sim), intent(in) :: sim
  type(particle_sim) :: restored
  type(particle_slot_registry) :: pool
  integer :: error,g,unit
  integer*8 :: tick
  character(len=160) :: filename
  if(rank==0) then
    call system_clock(count=tick)
    write(filename,'(A,I0,A)') 'reflection_',tick,'.h5'
  endif
  call MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,error)
  n_part_groups=3; part_groups_in_use(1:3)=sim%groups%id
  call H5open_f(error)
  call write_simulation_hdf5(sim,trim(filename))
  call MPI_BARRIER(MPI_COMM_WORLD,error)
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
  call initialize_particle_slot_registry(restored,[1,2,3],pool)
  call validate_particle_slot_registry(restored,pool,[1,2,3],error)
  call assert_equals(0,error,'restart reconstructed pools valid')
  do g=1,3
    call assert_equals(sim%groups(g)%id,restored%groups(g)%id,'restart group ID')
    call same_particles(sim%groups(g)%particles,restored%groups(g)%particles,'exact common restart')
    call assert_equals(count(restored%groups(g)%particles%i_elm==0),pool%pool(g)%n_free,'restart pool membership')
  enddo
  call MPI_BARRIER(MPI_COMM_WORLD,error)
  if(rank==0) then
    open(newunit=unit,file=trim(filename),status='old'); close(unit,status='delete')
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,error)
end subroutine
subroutine multiple_side_and_initialization()
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state
  integer :: k,error,side,n1,n2
  real*8 :: st(2),normal(3)
  call fixture(sim,passive,state,1)
  call advance_gk_neutrals(sim,passive,dt,1)
  ! A corner ELEMENT is legal when retained st identifies one wall side.
  side=0
  do k=1,4
    n1=elements%element(corner_element)%vertex(k)
    n2=elements%element(corner_element)%vertex(mod(k,4)+1)
    if(nodes%node(n1)%boundary>0.and.nodes%node(n2)%boundary>0) side=k
  enddo
  select case(side)
  case(1); st=[0.4d0,0.d0]
  case(2); st=[1.d0,0.4d0]
  case(3); st=[0.4d0,1.d0]
  case(4); st=[0.d0,0.4d0]
  end select
  normal=wall_normal_vector(nodes,elements,corner_element,st(1),st(2))
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    if(owns(1)) then
      p(1)%i_elm=-corner_element; p(1)%st=st; p(1)%v=-100.d0*normal
    endif
  end select
  call checked_reflect(sim,passive,state,error,.false.)
  call report_gk_neutral_reflection(sim,state)
  call assert_equals(1,state%global_counts(5),'unambiguous multiple-side element accepted/counted')
  do k=1,5
    call fixture(sim,passive,state,1)
    if(rank==0) then
      select case(k)
      case(1); sim%groups(3)%coupling_scheme='pic'
      case(2); sim%groups(3)%mass=-1.d0
      case(3); sim%groups(3)%mass=nan_value
      case(4)
        select type(p=>sim%groups(3)%particles)
        type is(particle_kinetic_leapfrog)
          p(1)%q=1
        end select
      case(5)
        deallocate(sim%groups(3)%particles)
      end select
    endif
    call initialize_gk_neutral_reflection(sim,3,.true.,state,error)
    call assert_true(error/=0,'invalid initialization rejected collectively')
  enddo
  call fixture(sim,passive,state,0)
  call initialize_gk_neutral_reflection(sim,0,.true.,state,error)
  call assert_true(error/=0,'nonpositive explicit group rejected')
end subroutine
! Analytic cubic side 2 used to test the production interp_RZ intersection.
! Other faces stay internal; the whole map keeps positive Jacobian.
subroutine wall_case(which,amplitude,roots,u0,du,target,vphi,reported_u)
  integer, intent(in) :: which
  real*8, intent(out) :: amplitude,roots(3),u0,du,target,vphi,reported_u
  amplitude=0.2d0;roots=[0.2d0,0.5d0,0.8d0]
  u0=0.3d0;du=0.6d0;target=2.5d0;vphi=300.d0;reported_u=0.003d0
  select case(which)
  case(2);vphi=0.d0
  case(3);amplitude=1.d0;vphi=10000.d0
  case(4)
    roots(3)=0.999d0;u0=0.9995d0;du=-0.1d0;vphi=0.d0;reported_u=0.5d0
  case(5,18)
    amplitude=-0.2d0;u0=0.1d0;du=0.8d0;vphi=0.d0;reported_u=0.5d0
  case(10)
    roots(3)=0.999d0;u0=0.8d0;du=0.4d0;target=2.48d0;vphi=0.d0
  case(11)
    u0=0.6d0;du=0.15d0;vphi=0.d0
    target=wall_radius((3.d0+sqrt(1.08d0))/6.d0,amplitude,roots)
  end select
end subroutine
pure function wall_radius(u,amplitude,roots) result(r)
  real*8, intent(in) :: u,amplitude,roots(3)
  real*8 :: r
  r=2.5d0+amplitude*product(u-roots)
end function
pure function wall_derivative(u,amplitude,roots) result(dr)
  real*8, intent(in) :: u,amplitude,roots(3)
  real*8 :: dr
  dr=amplitude*((u-roots(1))*(u-roots(2))+(u-roots(1))*(u-roots(3))+(u-roots(2))*(u-roots(3)))
end function
subroutine selected_fixture(sim,passive,state,which,number)
  type(particle_sim), intent(out) :: sim
  type(gk_neutral_state), intent(out) :: passive
  type(gk_neutral_reflection_state), intent(out) :: state
  integer, intent(in) :: which,number
  type(particle_kinetic_leapfrog) :: old
  integer :: e,n1,n2,j,k,slot
  real*8 :: amplitude,roots(3),u0,du,target,vphi,reported_u,a,b,m,r,z,z0,z1
  call fixture(sim,passive,state,0)
  call wall_case(which,amplitude,roots,u0,du,target,vphi,reported_u)
  e=wall_element;n1=elements%element(e)%vertex(2);n2=elements%element(e)%vertex(3)
  nodes%node(n1)%x(1,1,1)=wall_radius(0.d0,amplitude,roots)
  nodes%node(n2)%x(1,1,1)=wall_radius(1.d0,amplitude,roots)
  nodes%node(n1)%x(1,3,1)=wall_derivative(0.d0,amplitude,roots)/(3.d0*elements%element(e)%size(2,3))
  nodes%node(n2)%x(1,3,1)=-wall_derivative(1.d0,amplitude,roots)/(3.d0*elements%element(e)%size(3,3))
  call interp_RZ(nodes,elements,e,1.d0,0.d0,r,z0)
  call interp_RZ(nodes,elements,e,1.d0,1.d0,r,z1)
  do j=1,number
    if(.not.owns(j)) cycle
    old%i_elm=e;old%st(2)=u0;a=0.d0;b=1.d0
    do k=1,60
      m=0.5d0*(a+b);call interp_RZ(nodes,elements,e,m,u0,r,z)
      if(r<target) then
        a=m
      else
        b=m
      endif
    enddo
    old%st(1)=0.5d0*(a+b)
    call interp_RZ(nodes,elements,e,old%st(1),old%st(2),r,z)
    old%x=[r,z,7.2d0+0.1d0*j];old%v=[0.d0,du*(z1-z0)/dt,vphi*j]
    if(which==10) old%v(1)=(wall_radius(1.d0,amplitude,roots)-r)/(0.5d0*dt)
    old%q=0;old%weight=1.d10*j;old%i_life=70+j;old%t_birth=0.03125_4
    ! Simulate locator-returned loss metadata; trajectory provenance itself is exact.
    slot=passive%n_newly_lost_indices+1;passive%n_newly_lost_indices=slot
    passive%newly_lost_indices(slot)=j
    passive%newly_lost_events(slot)=gk_neutral_event(j,e,1,old%x,old%v,old%st,dt)
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      p(j)=old;call push_gk_neutral_ballistic(p(j),dt)
      p(j)%i_elm=-e;p(j)%st=[1.d0,reported_u]
    end select
  enddo
end subroutine
subroutine corrupt_selected(sim,passive,state,which)
  type(particle_sim), intent(inout) :: sim
  type(gk_neutral_state), intent(inout) :: passive
  type(gk_neutral_reflection_state), intent(inout) :: state
  integer, intent(in) :: which
  integer :: k,j,e
  if(rank/=min(1,tasks-1)) return
  k=1;j=passive%newly_lost_indices(k);e=-sim%groups(3)%particles(j)%i_elm
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    select case(which)
    case(8)
      passive%newly_lost_events(k)%v_old=[-100.d0,1.d0,0.d0]
      p(j)%x=passive%newly_lost_events(k)%x_old;p(j)%v=passive%newly_lost_events(k)%v_old
      call push_gk_neutral_ballistic(p(j),dt)
    case(9);elements%element(e)%neighbours(2)=interior_element
    case(10,11) ! Corner and tangency are supplied by the geometric fixture.
    case(12);passive%newly_lost_events(k)%marker=0
    case(13);passive%newly_lost_events(k)%x_old(1)=nan_value
    case(14);passive%newly_lost_events(k)%dt=2.d0*dt
    case(15);p(j)%x(1)=p(j)%x(1)+0.1d0
    case(16);state%trajectory_count=huge(0)
    end select
  end select
end subroutine
subroutine selected_wall_cases()
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state,before
  type(gk_neutral_event) :: event
  type(pcg32_rng) :: expected_rng
  class(particle_base), allocatable :: saved(:)
  integer, parameter :: cases(16)=[1,2,3,4,5,8,9,10,11,12,13,14,15,16,17,18]
  integer :: icase,which,j,k,e,error,number,it
  real*8, allocatable :: saved_x(:,:,:,:)
  integer :: saved_neighbours(4,elements%n_elements)
  real*8 :: amplitude,roots(3),u0,du,target,vphi,reported_u,a,b,lambda,u,r,z,chord
  real*8 :: normal(3),draw(2),expected_velocity(3),speed2
  allocate(saved_x(nodes%n_nodes,n_coord_tor,n_degrees,2))
  do j=1,nodes%n_nodes
    saved_x(j,:,:,:)=nodes%node(j)%x
  enddo
  do j=1,elements%n_elements
    saved_neighbours(:,j)=elements%element(j)%neighbours
  enddo
  do icase=1,size(cases)
    which=cases(icase)
    do j=1,nodes%n_nodes
      nodes%node(j)%x=saved_x(j,:,:,:)
    enddo
    do j=1,elements%n_elements
      elements%element(j)%neighbours=saved_neighbours(:,j)
    enddo
    number=1;if(which==1.or.which>=8) number=3
    call selected_fixture(sim,passive,state,which,number)
    if(which==17.or.which==18) then
      do j=1,number
        if(.not.owns(j)) cycle
        sim%groups(3)%particles(j)%i_elm=-(wall_element-3)
        sim%groups(3)%particles(j)%st=[0.d0,0.5d0]
        if(which==18) then
          sim%groups(3)%particles(j)%i_elm=-1
          sim%groups(3)%particles(j)%st=[0.5d0,0.d0]
        endif
      enddo
    endif
    if(which>=8.and.which<=16) then
      call corrupt_selected(sim,passive,state,which)
      call checked_reflect(sim,passive,state,error,.true.)
      cycle
    endif
    expected_rng=state%rng
    allocate(saved,source=sim%groups(3)%particles)
    call checked_reflect(sim,passive,state,error,.false.)
    if(error/=0) then
      deallocate(saved);cycle
    endif
    call report_gk_neutral_reflection(sim,state)
    call assert_equals(number,state%global_trajectory_count,'selected-wall correction count reduced once')
    call assert_equals(number,state%global_counts(2),'corrected incident markers reflected once')
    call assert_equals(0.5d0*number*(number+1)*1.d10,state%global_totals(1),0.d0,'global incident weight')
    call assert_equals(state%global_totals(1),state%global_totals(2),0.d0,'global corrected weight conserved')
    call assert_true(abs(state%global_totals(3)-state%global_totals(4))<= &
                     128*epsilon(1.d0)*state%global_totals(3),'global corrected energy conserved')
    call assert_equals(0,passive%n_newly_lost_indices,'corrected events consumed')
    call wall_case(which,amplitude,roots,u0,du,target,vphi,reported_u)
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      do k=1,count([(owns(j),j=1,number)])
        event=passive%newly_lost_events(k);j=event%marker;e=p(j)%i_elm
        a=0.d0;b=0.5d0
        if(which==4) b=0.05d0
        if(which==5.or.which==18) b=0.25d0
        ! Independent bisection of analytic wall minus curved Cartesian trajectory.
        do it=1,64
          lambda=0.5d0*(a+b);u=u0+du*lambda
          r=sqrt((event%x_old(1)+lambda*dt*event%v_old(1))**2+(lambda*dt*event%v_old(3))**2)
          if(wall_radius(u,amplitude,roots)>r) then
            a=lambda
          else
            b=lambda
          endif
        enddo
        lambda=0.5d0*(a+b);u=u0+du*lambda
        r=sqrt((event%x_old(1)+lambda*dt*event%v_old(1))**2+(lambda*dt*event%v_old(3))**2)
        z=event%x_old(2)+lambda*dt*event%v_old(2)
        call assert_true(norm2(p(j)%x(:2)-[r,z])<2.d-11,'first true ballistic/curved-wall intersection')
        call assert_equals(wall_element,e,'verified starting element becomes positive destination')
        call assert_equals(1.d0,p(j)%st(1),0.d0,'selected physical side 2 retained')
        call assert_equals(u,p(j)%st(2),2.d-10,'corrected local wall coordinate')
        if(which==3) then
          chord=event%x_old(1)+lambda*(saved(j)%x(1)-event%x_old(1))
          call assert_true(abs(chord-r)>1.d-3,'large toroidal velocity differs from straight RZ chord')
        endif
        if(which==4) call assert_true(p(j)%st(2)>0.998d0,'near-vertex crossing is unambiguous')
        if(which==5.or.which==18) &
          call assert_equals(0.2d0,p(j)%st(2),2.d-10,'first of two outward crossings chosen')
        normal=wall_normal_vector(nodes,elements,e,p(j)%st(1),p(j)%st(2))
        call expected_rng%next(draw)
        select type(h=>saved)
        type is(particle_kinetic_leapfrog)
          speed2=sum(h(j)%v**2);expected_velocity=sqrt(speed2)*sample_cosine(draw,normal)
          call assert_true(norm2(p(j)%v-expected_velocity)<1.d-10,'unchanged cosine law and two uniforms')
          call assert_true(abs(sum(p(j)%v**2)/speed2-1.d0)<128*epsilon(1.d0),'corrected energy conserved')
          call assert_true(p(j)%weight==h(j)%weight.and.p(j)%q==h(j)%q.and.p(j)%i_life==h(j)%i_life.and. &
                          p(j)%t_birth==h(j)%t_birth,'same marker weight/charge/life/birth')
          call assert_equals(h(j)%x(3),p(j)%x(3),0.d0,'continuous end phi convention unchanged')
        end select
        call assert_true(dot_product(p(j)%v,normal)>0.d0,'outgoing velocity inward at corrected wall')
      enddo
    end select
    call assert_true(all(transfer(expected_rng,[0_1],storage_size(expected_rng)/8)== &
                        transfer(state%rng,[0_1],storage_size(state%rng)/8)),'no extra intersection RNG draws')
    deallocate(saved);allocate(saved,source=sim%groups(3)%particles);before=state
    call checked_reflect(sim,passive,state,error,.false.)
    call same_particles(sim%groups(3)%particles,saved,'consumed intersection cannot reflect twice')
    call same_state(state,before,'consumed correction diagnostics and RNG unchanged')
    call report_gk_neutral_reflection(sim,state)
    call same_state(state,before,'corrected report does not double-reduce')
    if(which==1.or.which==5.or.which>=17) call restart_roundtrip(sim)
    if(which==2) then
      call advance_gk_neutrals(sim,passive,dt,1)
      if(owns(1)) call assert_true(sim%groups(3)%particles(1)%i_elm>0,'corrected neutral next inward outer step')
    endif
    deallocate(saved)
  enddo
  do j=1,nodes%n_nodes
    nodes%node(j)%x=saved_x(j,:,:,:)
  enddo
  do j=1,elements%n_elements
    elements%element(j)%neighbours=saved_neighbours(:,j)
  enddo
end subroutine

! Rejected starting-cell paths must cancel valid corrections on other ranks.
subroutine starting_element_rejections()
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state
  real*8, allocatable :: saved_x(:,:,:,:)
  integer :: which,fixture_kind,j,k,error,bad,returned
  integer :: saved_neighbours(4,elements%n_elements)
  bad=min(1,tasks-1)
  allocate(saved_x(nodes%n_nodes,n_coord_tor,n_degrees,2))
  do j=1,nodes%n_nodes
    saved_x(j,:,:,:)=nodes%node(j)%x
  enddo
  do j=1,elements%n_elements
    saved_neighbours(:,j)=elements%element(j)%neighbours
  enddo
  do which=1,7
    do j=1,nodes%n_nodes
      nodes%node(j)%x=saved_x(j,:,:,:)
    enddo
    do j=1,elements%n_elements
      elements%element(j)%neighbours=saved_neighbours(:,j)
    enddo
    fixture_kind=1
    if(rank==bad.and.which==4) fixture_kind=10
    if(rank==bad.and.which==5) fixture_kind=11
    call selected_fixture(sim,passive,state,fixture_kind,3)
    if(which==1) call corrupt_selected(sim,passive,state,12)
    if(which==2) call corrupt_selected(sim,passive,state,16)
    if(which==3) call corrupt_selected(sim,passive,state,8)
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      do k=1,passive%n_newly_lost_indices
        j=passive%newly_lost_indices(k);returned=wall_element-3
        p(j)%st=[0.d0,0.5d0]
        if(rank==bad.and.(which==3.or.which==5)) then
          returned=elements%element(wall_element)%neighbours(1);p(j)%st=[1.d0,0.997d0]
        endif
        p(j)%i_elm=-returned
        if(rank==bad.and.which==6) then
          ! A cyclic topology must reject even though the orbit has later wall roots.
          elements%element(wall_element)%neighbours(4)=wall_element
          passive%newly_lost_events(k)%st_old=[0.5d0,0.3d0]
          call interp_RZ(nodes,elements,wall_element,0.5d0,0.3d0, &
                         passive%newly_lost_events(k)%x_old(1),passive%newly_lost_events(k)%x_old(2))
          passive%newly_lost_events(k)%v_old=[-46000.d0,500.d0,10000.d0]
          p(j)%x=passive%newly_lost_events(k)%x_old;p(j)%v=passive%newly_lost_events(k)%v_old
          call push_gk_neutral_ballistic(p(j),dt)
        endif
        if(rank==bad.and.which==7) passive%newly_lost_events(k)%element=elements%n_elements+1
      enddo
    end select
    call checked_reflect(sim,passive,state,error,.true.)
  enddo
  do j=1,nodes%n_nodes
    nodes%node(j)%x=saved_x(j,:,:,:)
  enddo
  do j=1,elements%n_elements
    elements%element(j)%neighbours=saved_neighbours(:,j)
  enddo
end subroutine

! Analytic maps R=Rbase+alpha*Z**2 give curved, conforming interior faces.
subroutine traversal_cases()
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state,before
  type(particle_kinetic_leapfrog) :: source,incident
  type(pcg32_rng) :: expected_rng
  class(particle_base), allocatable :: saved(:)
  integer :: which,nr,j,k,e,left,next,error,expected_element,bad
  real*8 :: alpha,z,r,a,b,lambda,normal(3),draw(2),vexpected(3),speed2,limit
  logical :: reject_case
  bad=min(1,tasks-1)
  ! One/several/curved crossings, radial turn, malformed topology, bounds and cycles.
  do which=1,14
    nr=5;if(which==12) nr=max_trace_elements+6
    call mesh(nr)
    alpha=0.d0;if(which==3.or.which==9.or.which==14) alpha=0.2d0
    do j=1,nodes%n_nodes
      z=nodes%node(j)%x(1,1,2)
      nodes%node(j)%x(1,1,1)=nodes%node(j)%x(1,1,1)+alpha*z*z
      nodes%node(j)%x(1,2:3,1)=nodes%node(j)%x(1,2:3,1)+2*alpha*z*nodes%node(j)%x(1,2:3,2)
    enddo
    call fixture(sim,passive,state,0)
    left=wall_element-(nr-2)
    do j=1,3
      if(.not.owns(j)) cycle
      e=left;if(which==1.or.(which==12.and.rank/=bad)) e=wall_element-1
      if(which==4) e=wall_element
      if((which==13.and.rank==bad).or.which==14) e=left+1
      source%i_elm=e;source%st=[0.5d0,0.4d0]
      call interp_RZ(nodes,elements,e,source%st(1),source%st(2),r,z)
      source%x=[r,z,7.3d0];source%v=[22000.d0,1000.d0,4000.d0]
      if(which==4) source%v=[-46000.d0,500.d0,10000.d0]
      ! A radial turn returns through the shared face into an already visited cell.
      if(which==13.and.rank==bad) source%v=[-10000.d0,500.d0,10000.d0]
      if(which==14) source%v=[-30000.d0,1500.d0,30000.d0]
      if(which==10.and.rank==bad) then
        source%st=[0.5d0,0.5d0]
        call interp_RZ(nodes,elements,e,0.5d0,0.5d0,source%x(1),source%x(2))
        source%v=[10000.d0,10000.d0,0.d0]
      endif
      source%q=0;source%weight=j*1.d10;source%i_life=83+j;source%t_birth=0.0625_4
      k=passive%n_newly_lost_indices+1;passive%n_newly_lost_indices=k
      passive%newly_lost_indices(k)=j
      passive%newly_lost_events(k)=gk_neutral_event(j,e,1,source%x,source%v,source%st,dt)
      select type(p=>sim%groups(3)%particles)
      type is(particle_kinetic_leapfrog)
        p(j)=source;call push_gk_neutral_ballistic(p(j),dt)
        p(j)%i_elm=-left;p(j)%st=[0.d0,0.4d0]
      end select
    enddo
    if(rank==bad) then
      next=elements%element(left)%neighbours(2)
      select case(which)
      case(5);elements%element(left)%neighbours(2)=left
      case(6);elements%element(next)%neighbours(2)=left
      case(7);elements%element(left)%neighbours(2)=elements%n_elements+1
      case(8);elements%element(next)%neighbours(4)=0
      case(9);elements%element(next)%size(1,1)=elements%element(next)%size(1,1)+1.d-3
      case(11);state%trajectory_count=huge(0)
      end select
    endif
    reject_case=which>=5.and.which<=13
    expected_rng=state%rng
    allocate(saved,source=sim%groups(3)%particles)
    call checked_reflect(sim,passive,state,error,reject_case)
    if(reject_case.or.error/=0) then
      deallocate(saved);cycle
    endif
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      call assert_equals(0,passive%n_newly_lost_indices,'successful traversal consumes all events')
      do k=1,count([(owns(j),j=1,3)])
        j=passive%newly_lost_indices(k)
        source%x=passive%newly_lost_events(k)%x_old;source%v=passive%newly_lost_events(k)%v_old
        a=0.d0;b=1.d0
        if(which==4) b=-source%x(1)*source%v(1)/(dt*sum(source%v([1,3])**2))
        do e=1,64
          lambda=0.5d0*(a+b);z=source%x(2)+lambda*dt*source%v(2)
          r=sqrt((source%x(1)+lambda*dt*source%v(1))**2+(lambda*dt*source%v(3))**2)
          limit=2.5d0+alpha*z*z;if(which==4) limit=0.5d0
          if((which/=4.and.r<limit).or.(which==4.and.r>limit)) then
            a=lambda
          else
            b=lambda
          endif
        enddo
        expected_element=wall_element;normal=[-1.d0,2*alpha*z,0.d0]
        if(which==4) then
          expected_element=left;normal=[1.d0,0.d0,0.d0]
        endif
        normal=normal/norm2(normal)
        call assert_equals(expected_element,p(j)%i_elm,'earliest physical wall after interior crossings')
        call assert_true(norm2(p(j)%x(:2)-[r,z])<2.d-11,'independent curved traversal intersection')
        call assert_true(dot_product(p(j)%v,normal)>0.d0,'traced reflection points inward')
        select type(q=>saved)
        type is(particle_kinetic_leapfrog)
          speed2=sum(q(j)%v**2)
          call expected_rng%next(draw);vexpected=sqrt(speed2)*sample_cosine(draw,normal)
          call assert_true(norm2(vexpected-p(j)%v)<2.d-10,'traversal draws no extra random numbers')
          call assert_true(abs(sum(p(j)%v**2)/speed2-1.d0)<128*epsilon(1.d0),'traversal energy conserved')
          call assert_true(p(j)%weight==q(j)%weight.and.p(j)%q==q(j)%q.and. &
            p(j)%i_life==q(j)%i_life.and.p(j)%t_birth==q(j)%t_birth,'traversal marker identity preserved')
          call assert_equals(q(j)%x(3),p(j)%x(3),0.d0,'traversal keeps existing continuous phi convention')
        end select
      enddo
    end select
    call assert_true(all(transfer(expected_rng,[0_1])==transfer(state%rng,[0_1])),'traversal RNG stream')
    call report_gk_neutral_reflection(sim,state)
    call assert_equals(3,state%global_counts(2),'three reflected physical markers across MPI')
    call assert_equals(6.d10,state%global_totals(1),0.d0,'traversal incident weight')
    call assert_equals(state%global_totals(1),state%global_totals(2),0.d0,'traversal reflected weight')
    call assert_true(abs(state%global_totals(3)-state%global_totals(4))<= &
      128*epsilon(1.d0)*state%global_totals(3),'traversal global energy conserved')
    before=state;deallocate(saved);allocate(saved,source=sim%groups(3)%particles)
    call checked_reflect(sim,passive,state,error,.false.)
    call same_particles(saved,sim%groups(3)%particles,'no repeated traversal/reflection')
    call same_state(state,before,'no repeat traversal RNG or diagnostics')
    if(which==3.or.which==4.or.which==14) call restart_roundtrip(sim)
    deallocate(saved)
  enddo
  call mesh()
end subroutine

! Analytic warped square: R=R_linear+alpha*Z**2. Locate sources independently
! of find_RZ_nearby; the same mapping supplies reference endpoint coordinates.
subroutine recovery_mesh(alpha)
  real*8, intent(in) :: alpha
  integer :: j
  real*8 :: z
  call mesh()
  do j=1,nodes%n_nodes
    z=nodes%node(j)%x(1,1,2)
    nodes%node(j)%x(1,1,1)=nodes%node(j)%x(1,1,1)+alpha*z*z
    nodes%node(j)%x(1,2:3,1)=nodes%node(j)%x(1,2:3,1)+2*alpha*z*nodes%node(j)%x(1,2:3,2)
  enddo
end subroutine
subroutine analytic_location(source,alpha)
  type(particle_kinetic_leapfrog), intent(inout) :: source
  real*8, intent(in) :: alpha
  integer :: e
  real*8 :: r,z,st(2)
  source%i_elm=0
  do e=1,elements%n_elements
    call interp_RZ(nodes,elements,e,0.5d0,0.5d0,r,z)
    st=[(source%x(1)-alpha*source%x(2)**2-(r-alpha*z*z)+0.25d0)/0.5d0, &
        (source%x(2)-z+0.25d0)/0.5d0]
    if(all(st>=0.d0).and.all(st<=1.d0)) then
      source%i_elm=e;source%st=st
      return
    endif
  enddo
end subroutine
subroutine false_loss_recovery_cases()
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state,before
  type(particle_kinetic_leapfrog) :: source,expected
  type(particle_slot_registry) :: pool,pool_before
  class(particle_base), allocatable :: ions(:),electrons(:)
  integer :: which,j,steps,error,nloss,expected_loss
  real*8 :: alpha,wl,el
  alpha=0.2d0
  do which=1,5
    call recovery_mesh(alpha)
    call fixture(sim,passive,state,0)
    steps=1;if(which==4) steps=2
    source%x=[2.482d0,-0.4d0,7.3d0];source%v=[0.d0,8000.d0,0.d0]
    if(which==2) then
      source%x=[2.468d0,-0.8d0,7.3d0];source%v=[0.d0,16000.d0,3000.d0]
    else if(which==3.or.which==5) then
      source%x=[2.472d0,-0.4d0,7.3d0];source%v=[0.d0,8000.d0,3000.d0]
    else if(which==4) then
      source%x=[2.56d0,-0.6d0,7.3d0];source%v=[-640.d0,4000.d0,0.d0]
    endif
    source%q=0;source%weight=6.d10;source%i_life=41;source%t_birth=0.03125_4
    call analytic_location(source,alpha)
    call assert_true(source%i_elm>0,'analytic recovered-source location')
    expected=source
    do j=1,steps
      call push_gk_neutral_ballistic(expected,dt)
    enddo
    call analytic_location(expected,alpha)
    call assert_true(expected%i_elm>0,'analytic interior endpoint after all substeps')
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      if(owns(1)) p(1)=source
      if(which==5.and.owns(2)) then
        p(2)=source;p(2)%i_elm=wall_element;p(2)%st=[0.999d0,0.4d0]
        call interp_RZ(nodes,elements,wall_element,p(2)%st(1),p(2)%st(2),p(2)%x(1),p(2)%x(2))
        p(2)%v=[1000.d0,200.d0,300.d0];p(2)%weight=2.d10
      endif
    end select
    call initialize_particle_slot_registry(sim,[1,2,3],pool)
    call snapshot_slot_registry(pool,pool_before)
    allocate(ions,source=sim%groups(1)%particles);allocate(electrons,source=sim%groups(2)%particles)
    before=state
    call advance_gk_neutrals(sim,passive,dt,steps)
    call same_state(state,before,'passive recovery consumes no reflection RNG or diagnostics')
    call assert_equals(merge(1,0,owns(1)),passive%recovered_count,'one false loss recovered on its owner')
    expected_loss=0;if(which==5.and.owns(2)) expected_loss=1
    call get_gk_neutral_loss_totals(passive,nloss,wl,el)
    call assert_equals(expected_loss,nloss,'only physical wall crossings counted lost')
    call assert_equals(expected_loss,passive%verified_wall_count,'verified physical wall counter')
    call assert_equals(0,passive%unresolved_count,'no unresolved recovery')
    call assert_equals(expected_loss,passive%n_newly_lost_indices,'no recovered marker in event list')
    call assert_equals(real(expected_loss,8)*2.d10,wl,0.d0,'false loss contributes no lost weight')
    if(expected_loss==0) call assert_equals(0.d0,el,0.d0,'false loss contributes no lost energy')
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      if(owns(1)) then
        call assert_equals(expected%i_elm,p(1)%i_elm,'verified final element, including remaining substeps')
        call assert_true(norm2(expected%st-p(1)%st)<1.d-11,'analytic endpoint coordinates')
        call assert_true(bits(expected%x,p(1)%x).and.bits(expected%v,p(1)%v),'ballistic x/v exactly unchanged')
        call assert_equals(expected%weight,p(1)%weight,0.d0,'recovery preserves weight exactly')
        call assert_true(sum(expected%v**2)==sum(p(1)%v**2),'recovery preserves endpoint energy exactly')
        call assert_true(p(1)%q==source%q.and.p(1)%i_life==source%i_life.and. &
          p(1)%t_birth==source%t_birth,'recovery preserves charge, life and birth')
        call assert_true(passive%max_traversed>=2,'recovery crosses interior interfaces')
      endif
      if(expected_loss==1) then
        call assert_equals(2,passive%newly_lost_indices(1),'mixed event belongs only to genuine wall marker')
        call assert_true(p(2)%i_elm<0,'genuine loss remains negative before reflection')
        call assert_equals(0.5d0*p(2)%weight*sim%groups(3)%mass*ATOMIC_MASS_UNIT*sum(p(2)%v**2), &
                           el,abs(el)*1.d-15,'physical lost energy')
      endif
    end select
    call assert_slot_registry_equal(pool,pool_before,'recovery preserves complete live registry')
    call validate_particle_slot_registry(sim,pool,[1,2,3],error)
    call assert_equals(0,error,'registry membership remains correct before reflection')
    call checked_reflect(sim,passive,state,error,.false.)
    if(which/=5) call same_state(state,before,'no reflection or RNG draw after false loss')
    call assert_equals(expected_loss,state%counts(2),'reflection sees only genuine local events')
    call assert_slot_registry_equal(pool,pool_before,'mixed recovery/reflection preserves registry')
    call same_particles(sim%groups(1)%particles,ions,'recovery protects charged ions')
    call same_particles(sim%groups(2)%particles,electrons,'recovery protects charged electrons')
    call report_gk_neutral_diagnostics(sim,passive)
    call restart_roundtrip(sim)
    deallocate(ions,electrons)
  enddo
  call mesh()
end subroutine
subroutine interior_endpoint_geometry_cases()
  type(particle_kinetic_leapfrog) :: source,endpoint
  type(side_correction) :: hit
  integer :: which,start,neighbour
  real*8 :: alpha,expected_lambda
  do which=1,8
    alpha=0.d0;if(which==2) alpha=0.2d0
    call recovery_mesh(alpha)
    source%x=[1.25d0,0.1d0,7.3d0];source%v=[-10000.d0,500.d0,10000.d0]
    if(which==2) then ! Both endpoints inside, but a physical exit precedes re-entry.
      source%x=[2.52d0,-0.8d0,7.3d0];source%v=[0.d0,16000.d0,0.d0]
    else if(which==3.or.which==4) then
      source%x=[1.25d0,-0.25d0,7.3d0];source%v=[10000.d0,10000.d0,0.d0]
      if(which==4) source%v(2)=10000.d0*(1.d0+1.d-5)
    else if(which==5.or.which==6) then
      source%x=[1.25d0,-0.1d0,7.3d0];source%v=[0.d0,(0.1d0-1.d-8)/dt,0.d0]
      if(which==6) source%v(2)=0.1d0/dt
    endif
    call analytic_location(source,alpha);start=source%i_elm
    if(which==7) then
      neighbour=elements%element(start)%neighbours(4)
      elements%element(start)%neighbours(4)=elements%n_elements+1
    endif
    if(which==8) source%v(1)=nan_value
    endpoint=source;call push_gk_neutral_ballistic(endpoint,dt)
    call verify_neutral_trajectory(nodes,elements,start,source%st,source%x,source%v,dt,hit)
    select case(which)
    case(1)
      call assert_equals(INTERIOR_ENDPOINT,hit%status,'radial revisit ends inside')
      call assert_equals(start,hit%element,'legitimate chronological return to original cell')
      call assert_equals(3,hit%n_visits,'two crossings and three visits including revisited cell')
      call assert_true(hit%trace_fractions(2)>hit%trace_fractions(1),'revisit makes chronological progress')
    case(2)
      expected_lambda=(0.8d0-sqrt(0.02d0/alpha))/1.6d0
      call assert_equals(WALL_CROSSING,hit%status,'interior endpoint does not override earlier physical exit')
      call assert_equals(expected_lambda,hit%lambda,1.d-11,'first wall crossing before re-entry')
      call assert_true(endpoint%x(1)<2.5d0+alpha*endpoint%x(2)**2,'exit/re-entry endpoint is inside')
    case(3,7,8)
      call assert_equals(UNRESOLVED,hit%status,'ambiguous corner, invalid connectivity or nonfinite state rejects')
      call assert_true(len_trim(hit%reason)>0,'unresolved geometry has a reason')
    case(4,5,6)
      call assert_equals(INTERIOR_ENDPOINT,hit%status,'unambiguous near corner or near/exact face endpoint')
      call assert_true(hit%max_residual<1.d-10,'independent endpoint reconstruction residual')
    end select
  enddo
  call mesh()
end subroutine

! Fault injection affects only the locator's spatial index, not the FE mesh.
! With zero local iterations the unchanged global locator sees no candidates.
! This deliberately exercises interior/unresolved zero returns without a production hook.
subroutine empty_locator_index()
  use mod_element_rtree, only: populate_element_rtree
  type(type_element_list) :: empty
  empty%n_elements=0
  call populate_element_rtree(nodes,empty)
end subroutine

subroutine zero_locator_cases()
  use phys_module, only: find_RZ_nearby_iter
  use mod_element_rtree, only: populate_element_rtree
  use mod_find_rz_nearby, only: find_RZ_nearby
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state,before
  type(particle_slot_registry) :: pool,pool_before
  type(particle_kinetic_leapfrog) :: source,endpoint,raw
  type(side_correction) :: hit
  type(pcg32_rng) :: rng
  integer :: which,j,k,n,steps,saved_iter,error,nloss,want
  real*8 :: alpha,wl,el,draw(2),expected_v(3)
  logical :: wall,reject_case
  saved_iter=find_RZ_nearby_iter
  do which=1,5
    alpha=0.d0;if(which==3) alpha=0.2d0
    call recovery_mesh(alpha)
    call fixture(sim,passive,state,0)
    source%x=[1.25d0,-0.11d0,7.3d0];source%v=[20000.d0,500.d0,1000.d0]
    if(which==2) source%v=[1000.d0,500.d0,1000.d0]
    if(which==3) then
      source%x=[2.52d0,-0.8d0,7.3d0];source%v=[0.d0,16000.d0,0.d0]
    endif
    if(which==4) then
      source%x=[1.25d0,0.1d0,7.3d0];source%v=[-10000.d0,500.d0,10000.d0]
    endif
    source%q=0;source%weight=6.d10;source%i_life=41;source%t_birth=0.03125_4
    call analytic_location(source,alpha)
    endpoint=source;call push_gk_neutral_ballistic(endpoint,dt)
    call verify_neutral_trajectory(nodes,elements,source%i_elm,source%st,source%x,source%v,dt,hit)
    wall=which==1.or.which==3.or.which==5
    call assert_equals(merge(WALL_CROSSING,INTERIOR_ENDPOINT,wall),hit%status,'zero fixture geometry classification')
    if(which==4) then
      call assert_equals(3,hit%n_visits,'zero interior trajectory revisits starting element')
      call assert_equals(source%i_elm,hit%element,'zero revisit returns to its start')
    endif
    steps=1;if(which==2) steps=2
    n=1;if(which==5) n=2
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      p(6)%i_elm=0
      do j=1,n
        if(owns(j)) p(j)=source
      enddo
    end select
    call initialize_particle_slot_registry(sim,[1,2,3],pool)
    call snapshot_slot_registry(pool,pool_before)
    find_RZ_nearby_iter=0
    ! Case 1 is a natural outside-mesh fallback; the others inject a missing index.
    if(which/=1) call empty_locator_index()
    raw=endpoint
    call find_RZ_nearby(nodes,elements,source%x(1),source%x(2),source%st(1),source%st(2),source%i_elm, &
      raw%x(1),raw%x(2),raw%st(1),raw%st(2),raw%i_elm,error,phi=raw%x(3))
    call assert_equals(0,raw%i_elm,'real locator returns zero')
    call assert_equals(99,error,'real global locator exhausts candidates')
    before=state
    call advance_gk_neutrals(sim,passive,dt,steps)
    call same_state(state,before,'zero verification consumes no RNG')
    call populate_element_rtree(nodes,elements)
    find_RZ_nearby_iter=saved_iter
    want=count([(owns(j),j=1,n)])
    call assert_equals(want*steps,passive%zero_count,'all zero substeps classified')
    call assert_equals(merge(want,0,wall),passive%zero_wall_count,'zero wall count')
    call assert_equals(merge(0,want*steps,wall),passive%zero_interior_count,'zero interior count')
    call assert_equals(0,passive%zero_unresolved_count,'no unresolved zero fixture')
    call get_gk_neutral_loss_totals(passive,nloss,wl,el)
    call assert_equals(merge(want,0,wall),nloss,'zero physical loss count')
    call assert_equals(nloss,passive%n_newly_lost_indices,'zero exactly one event per wall marker')
    call assert_equals(nloss*source%weight,wl,0.d0,'zero physical lost weight')
    call assert_true(abs(el-nloss*0.5d0*source%weight*sim%groups(3)%mass*ATOMIC_MASS_UNIT*sum(endpoint%v**2)) &
      <=128*epsilon(1.d0)*max(el,tiny(1.d0)),'zero physical lost energy')
    if(steps==2) call push_gk_neutral_ballistic(endpoint,dt)
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      do j=1,n
        if(.not.owns(j)) cycle
        call assert_true(bits(endpoint%x,p(j)%x).and.bits(endpoint%v,p(j)%v),'zero keeps exact ballistic x/v')
        if(wall) then
          call assert_equals(-hit%element,p(j)%i_elm,'zero authoritative first physical wall')
          call assert_true(bits(hit%st,p(j)%st),'zero authoritative wall coordinates')
          k=findloc(passive%newly_lost_indices(:passive%n_newly_lost_indices),j,dim=1)
          call assert_true(bits(source%x,passive%newly_lost_events(k)%x_old).and. &
            bits(source%st,passive%newly_lost_events(k)%st_old).and. &
            bits(source%v,passive%newly_lost_events(k)%v_old),'zero original event provenance')
        else
          call analytic_location(endpoint,alpha)
          call assert_equals(endpoint%i_elm,p(j)%i_elm,'zero recovered final element after remaining substeps')
          call assert_true(norm2(endpoint%st-p(j)%st)<1.d-11,'zero recovered endpoint coordinates')
        endif
      enddo
      if(which==5.and.owns(2)) p(2)%q=1
    end select
    call assert_slot_registry_equal(pool,pool_before,'zero passive path preserves live pool')
    reject_case=which==5
    rng=state%rng
    call checked_reflect(sim,passive,state,error,reject_case)
    if(.not.reject_case) then
      select type(p=>sim%groups(3)%particles)
      type is(particle_kinetic_leapfrog)
        if(wall.and.owns(1)) then
          call rng%next(draw);expected_v=norm2(endpoint%v)*sample_cosine(draw,hit%normal)
          call assert_true(norm2(expected_v-p(1)%v)<2.d-10,'zero reflection consumes exactly two uniforms')
          call assert_true(dot_product(p(1)%v,hit%normal)>0.d0,'zero reflected velocity inward')
          call assert_true(norm2(p(1)%x(:2)-hit%wall)<1.d-13,'zero reflection reconstructs verified wall')
          call assert_true(abs(sum(p(1)%v**2)/sum(endpoint%v**2)-1.d0)<128*epsilon(1.d0),'zero energy conserved')
        endif
        if(owns(1)) call assert_true(p(1)%weight==source%weight.and.p(1)%q==0.and. &
          p(1)%i_life==source%i_life.and.p(1)%t_birth==source%t_birth,'zero exact marker identity')
      end select
      call assert_true(all(transfer(rng,[0_1])==transfer(state%rng,[0_1])),'zero exact RNG stream')
      before=state
      call checked_reflect(sim,passive,state,error,.false.)
      call same_state(state,before,'zero event cannot reflect twice')
      call restart_roundtrip(sim)
    endif
    call assert_slot_registry_equal(pool,pool_before,'zero reflection preserves live pool')
  enddo
  call mesh()
end subroutine

subroutine starting_coordinate_roundoff_cases()
  type(particle_kinetic_leapfrog) :: source,original
  type(side_correction) :: hit
  real*8 :: r,z
  integer :: which,expected,j
  do which=1,12
    call mesh()
    source%i_elm=1;source%st=[0.524271007149745d0,-3.469446951953614d-18]
    call interp_RZ(nodes,elements,1,source%st(1),0.d0,r,z)
    source%x=[r,z,0.7d0];source%v=[50.d0,500.d0,0.d0]
    expected=INTERIOR_ENDPOINT
    if(which==2.or.which==4.or.which==6) then
      source%i_elm=wall_element;source%st=[1.d0+epsilon(1.d0),0.4d0]
      call interp_RZ(nodes,elements,wall_element,1.d0,0.4d0,source%x(1),source%x(2))
      source%v=[-500.d0,50.d0,0.d0]
    endif
    select case(which)
    case(3,4)
      source%v=-source%v;expected=WALL_CROSSING
    case(5)
      source%st(2)=-1.d-10;expected=UNRESOLVED
    case(6)
      source%st(1)=1.d0+1.d-10;expected=UNRESOLVED
    case(7)
      source%x(2)=source%x(2)+1.d-8;expected=UNRESOLVED
    case(8)
      do j=1,nodes%n_nodes
        nodes%node(j)%x(:,:,2)=-nodes%node(j)%x(:,:,2)
      enddo
      source%x(2)=-source%x(2);expected=UNRESOLVED
    case(9,10)
      source%v=[500.d0,0.d0,0.d0]
      if(which==10) source%v(2)=1.d-12
      expected=UNRESOLVED
    case(11)
      source%st(2)=-16.d0*epsilon(1.d0)
    case(12)
      source%st(2)=-64.d0*epsilon(1.d0);expected=UNRESOLVED
    end select
    original=source
    call verify_neutral_trajectory(nodes,elements,source%i_elm,source%st,source%x,source%v,dt,hit)
    call assert_equals(expected,hit%status,'precision-scaled starting coordinate classification')
    call assert_true(bits(source%st,original%st).and.bits(source%x,original%x).and. &
      bits(source%v,original%v),'roundoff clips only the private verification copy')
    if(expected/=UNRESOLVED) then
      call assert_true(hit%start_clamped,'validated roundoff clipping reported')
      call interp_RZ(nodes,elements,hit%element,hit%st(1),hit%st(2),r,z)
      call assert_true(norm2([r,z]-hit%wall)<1.d-13,'roundoff result reconstructs FE geometry')
      if(expected==WALL_CROSSING) call assert_equals(0.d0,hit%lambda,1.d-12,'immediate outward physical wall departure')
    else
      call assert_true(len_trim(hit%reason)>0,'invalid or tangent start remains controlled')
    endif
  enddo
  call mesh()
end subroutine

! Exactly at the limit succeeds; the next cell remains a controlled failure.
! These are real, reciprocal, positive-Jacobian meshes; no malformed connectivity
! or tolerance change is used to manufacture the limit rejection.
subroutine traversal_limit_cases()
  type(particle_kinetic_leapfrog) :: source
  type(side_correction) :: hit
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state
  integer :: visits,nr,e,error,which
  real*8 :: r,z,wall_fraction
  call assert_equals(128,max_trace_elements,'explicit authorized traversal limit')
  do which=1,3
    visits=65
    if(which==2) visits=128
    if(which==3) visits=129
    nr=visits+1
    call mesh(nr)
    e=wall_element-(nr-2)
    source%i_elm=e;source%st=[0.5d0,0.4d0]
    call interp_RZ(nodes,elements,e,source%st(1),source%st(2),r,z)
    source%x=[r,z,0.7d0];source%v=[22000.d0,100.d0,0.d0]
    source%q=0;source%weight=1.d10
    wall_fraction=(2.5d0-r)/(dt*source%v(1))
    call verify_neutral_trajectory(nodes,elements,e,source%st,source%x,source%v,dt,hit)
    call assert_equals(min(visits,128),hit%n_visits,'bounded traversal visit count')
    call assert_true(all(hit%trace_fractions(2:hit%n_visits)>hit%trace_fractions(1:hit%n_visits-1)), &
      'every recorded exit fraction strictly increases')
    if(which<3) then
      call assert_equals(WALL_CROSSING,hit%status,'long valid first-wall crossing')
      call assert_equals(wall_element,hit%element,'long traversal reaches correct wall element')
      call assert_equals(wall_fraction,hit%lambda,1.d-12,'analytic long-traversal wall fraction')
      call fixture(sim,passive,state,0)
      select type(p=>sim%groups(3)%particles)
      type is(particle_kinetic_leapfrog)
        p(6)%i_elm=0
        if(owns(1)) p(1)=source
      end select
      call advance_gk_neutrals(sim,passive,dt,1)
      call assert_equals(merge(visits,0,owns(1)),passive%max_successful_traversed,'successful maximum at 65 and 128')
      call assert_true(sum(passive%visit_histogram(:,1))==merge(1_8,0_8,owns(1)), &
        'one classified trajectory in successful histogram')
      call assert_true(sum(passive%visit_histogram(:,2))==0_8,'no unresolved traversal histogram entry')
      call checked_reflect(sim,passive,state,error,.false.)
      call report_gk_neutral_diagnostics(sim,passive)
    else
      call assert_equals(UNRESOLVED,hit%status,'129th-element wall exceeds the new bound')
      call assert_true(index(hit%reason,'element traversal limit exceeded (128 elements)')==1, &
        'explicit new-limit diagnostic')
      call assert_equals(0,hit%side,'no fabricated wall side after exhaustion')
      call assert_equals(0.d0,hit%lambda,0.d0,'no fabricated wall intersection after exhaustion')
      call assert_true(hit%trace_fractions(128)<wall_fraction,'real wall lies beyond the last visited element')
    endif
  enddo
  call mesh()
end subroutine

! Independent analytic curved map; references come from scalar bisection of the
! original ballistic radius, never the production Z-eliminated polynomial.
subroutine small_vertical_velocity_cases()
  type(particle_kinetic_leapfrog) :: source,impact
  type(side_correction) :: hit,repeat
  real*8, parameter :: vz(4)=[0.d0,7.649253735366487d0,-7.649253735366487d0,500.d0]
  real*8 :: z,r,a,b,l,limit,alpha
  integer :: j,k,it
  do k=1,size(vz)
    call mesh();alpha=0.2d0
    do j=1,nodes%n_nodes
      z=nodes%node(j)%x(1,1,2)
      nodes%node(j)%x(1,1,1)=nodes%node(j)%x(1,1,1)+alpha*z*z
      nodes%node(j)%x(1,2:3,1)=nodes%node(j)%x(1,2:3,1)+2*alpha*z*nodes%node(j)%x(1,2:3,2)
    enddo
    source%i_elm=wall_element;source%st=[0.5d0,0.4d0]
    call interp_RZ(nodes,elements,wall_element,source%st(1),source%st(2),r,z)
    source%x=[r,z,0.7d0];source%v=[5000.d0,vz(k),4000.d0]
    call verify_neutral_trajectory(nodes,elements,source%i_elm,source%st,source%x,source%v,dt,hit)
    call assert_equals(WALL_CROSSING,hit%status,'zero/small signed/ordinary vertical velocity crosses curved wall')
    a=0.d0;b=1.d0
    do it=1,64
      l=(a+b)*0.5d0;z=source%x(2)+l*dt*source%v(2)
      r=sqrt((source%x(1)+l*dt*source%v(1))**2+(l*dt*source%v(3))**2);limit=2.5d0+alpha*z*z
      if(r<limit)then
        a=l
      else
        b=l
      endif
    enddo
    call assert_equals((a+b)*0.5d0,hit%lambda,2.d-10,'independent ballistic curved-wall fraction')
    impact=source;call push_gk_neutral_ballistic(impact,dt*hit%lambda)
    call assert_true(norm2(impact%x(:2)-hit%wall)<2.d-10,'small vertical velocity original 2D residual')
    call assert_true(hit%vn<0.d0.and.hit%refinement_counts(2)==0,'unambiguous outward crossing')
    call verify_neutral_trajectory(nodes,elements,source%i_elm,source%st,source%x,source%v,dt,repeat)
    call assert_true(all(hit%st==repeat%st).and.hit%lambda==repeat%lambda,'repeat verifier has no mutable geometry state')
    if(k==1.or.k==4) call assert_equals(0,hit%refinement_counts(1),'already accurate ordinary path does not refine')
  enddo
  call mesh()
end subroutine

#ifdef GK_NEUTRAL_REFINEMENT_TEST
! Exercise private safeguards with deliberately perturbed seeds. This test-only
! build uses the identical production helper; references use independent REAL(16)
! analytic geometry and never polynomial candidates or Newton iteration.
subroutine refinement_safeguard_cases()
  type(trajectory_origin) :: origin
  real*8 :: alpha,z,r,u,lambda,residual,ur,lr,uref,lref,st(2),wall(2),tangent(2),normal(3)
  real*8 :: zlo,zhi,rate,root_gap
  real(kind=16) :: a,b,lq,zq,rq,aa,bb,cc,disc,rootsq(2),v0(3),x0(3),tq,aq
  integer :: which,j,side,it,iterations
  logical :: ok,valid
  character(len=160) :: reason
  do which=1,13
    call mesh();alpha=0.2d0
    if(which==9.or.which==10) alpha=0.d0
    do j=1,nodes%n_nodes
      z=nodes%node(j)%x(1,1,2)
      nodes%node(j)%x(1,1,1)=nodes%node(j)%x(1,1,1)+alpha*z*z
      nodes%node(j)%x(1,2:3,1)=nodes%node(j)%x(1,2:3,1)+2*alpha*z*nodes%node(j)%x(1,2:3,2)
    enddo
    origin%element=wall_element;origin%st_old=[0.5d0,0.4d0];origin%dt=dt
    call interp_RZ(nodes,elements,wall_element,0.5d0,0.4d0,r,z)
    origin%x_old=[r,z,0.7d0];origin%v_old=[5000.d0,0.d0,4000.d0];side=2
    if(which==2) origin%v_old(2)=1.d-8
    if(which==3) origin%v_old(2)=-1.d-8
    if(which==4.or.which==8) origin%v_old(2)=500.d0
    if(which==5) origin%v_old=[-10000.d0,0.d0,40000.d0]
    if(which==6.or.which==7) then
      origin%v_old=[-20000.d0,0.d0,30000.d0];side=4
    endif
    tq=real(dt,16);v0=real(origin%v_old,16);x0=real(origin%x_old,16);aq=real(alpha,16)
    call side_geometry(nodes,elements,wall_element,side,0.d0,st,wall,tangent,normal,valid);zlo=wall(2)
    call side_geometry(nodes,elements,wall_element,side,1.d0,st,wall,tangent,normal,valid);zhi=wall(2)
    root_gap=1.d0
    if(which==6.or.which==7) then
      aa=tq*tq*(v0(1)**2+v0(3)**2);bb=2*x0(1)*tq*v0(1)
      rq=2._16+aq*x0(2)**2;cc=x0(1)**2-rq**2;disc=bb*bb-4*aa*cc
      rootsq=[(-bb-sqrt(disc))/(2*aa),(-bb+sqrt(disc))/(2*aa)]
      lq=rootsq(which-5);root_gap=real(rootsq(2)-rootsq(1),8)
    else
      a=0._16;b=1._16
      do it=1,120
        lq=(a+b)/2;zq=x0(2)+lq*tq*v0(2)
        rq=sqrt((x0(1)+lq*tq*v0(1))**2+(lq*tq*v0(3))**2)
        if(rq<2.5_16+aq*zq*zq)then
          a=lq
        else
          b=lq
        endif
      enddo
      lq=(a+b)/2
    endif
    zq=x0(2)+lq*tq*v0(2)
    lref=real(lq,8);uref=real((zq-real(zlo,16))/real(zhi-zlo,16),8)
    u=uref+1.d-5;lambda=lref+1.d-4;ur=0.25d0;lr=min(0.25d0,0.25d0*root_gap)
    if(which==8)then
      u=uref+0.08d0;lambda=lref+0.12d0
    endif
    if(which==9.or.which==10)then
      origin%v_old=[0.d0,500.d0,0.d0]
      if(which==10) origin%v_old(1)=1.d-12
    endif
    if(which==11)ur=1.d-7
    if(which==12)u=nan_value
    if(which==13)then
      lambda=0.5d0;origin%v_old=[-5000.d0,0.d0,0.d0]
    endif
    call refine_side_intersection(nodes,elements,wall_element,side,origin,dt,u,lambda,ur,lr, &
                                  2.d-12,iterations,residual,ok,reason)
    if(which<=8)then
      call assert_true(ok,'local safeguarded refinement converges on independently identified root')
      call assert_equals(uref,u,5.d-12,'refined side coordinate agrees with REAL(16) reference')
      call assert_equals(lref,lambda,5.d-12,'refined fraction agrees with REAL(16) reference')
      call assert_true(residual<=2.d-12,'original two-coordinate residual passes unchanged criterion')
      call assert_true(iterations>0.and.iterations<=8,'refinement iterations bounded')
      if(which==8)call assert_true(iterations>=2,'nonlinear seed requires several safeguarded iterations')
    else
      call assert_true(.not.ok,'singular/ill-conditioned/nonfinite/branch-jumping seed rejects')
      call assert_true(len_trim(reason)>0,'failed refinement retains explicit reason')
      if(which==9.or.which==10)call assert_true(index(reason,'Jacobian')>0,'poor Jacobian identified')
      if(which==11)call assert_true(index(reason,'root identity')>0,'candidate cannot jump its isolation neighbourhood')
    endif
    write(*,'(A,3I8,5ES26.17,1X,L1,1X,A)')'REFINEMENT UNIT ',which,side,iterations,u,lambda, &
      residual,u-uref,lambda-lref,ok,trim(reason)
  enddo
  call mesh()
end subroutine
#endif

subroutine run_reflection_tests()
  call mesh()
#ifdef GK_NEUTRAL_REFINEMENT_TEST
  call run_test_case(refinement_safeguard_cases,'private refinement safeguards and REAL(16) references')
#endif
  call run_test_case(small_vertical_velocity_cases,'curved wall with zero, small signed and ordinary vZ')
  call run_test_case(traversal_limit_cases,'65 and 128 visits succeed, 129 visits remain unresolved')
  call run_test_case(zero_locator_cases,'zero locator classification, reflection, substeps and MPI rollback')
  call run_test_case(starting_coordinate_roundoff_cases,'FE-validated starting-coordinate roundoff')
  call run_test_case(noops_and_absorption,'disabled/absent/empty neutral reflection')
  call run_test_case(reflections,'diffuse identity, energy, replay, repeated impacts and restart')
  call run_test_case(grazing_and_zero,'grazing tolerance and zero-speed RNG')
  call run_test_case(multiple_side_and_initialization,'boundary-side resolution and initialization validation')
  call run_test_case(rejections,'global rollback of every particle/pool/RNG/diagnostic/event')
  call run_test_case(selected_wall_cases,'selected curved-wall intersections, first exit and transactional rollback')
  call run_test_case(starting_element_rejections,'starting-element global rollback and invalid topology')
  call run_test_case(false_loss_recovery_cases,'passive false loss recovery, mixed events and restart')
  call run_test_case(interior_endpoint_geometry_cases,'interior endpoints, revisits, exit/re-entry and unresolved geometry')
  call run_test_case(traversal_cases,'bounded traversal, curved interior faces, corners, loops and global rollback')
end subroutine
subroutine run_reflection_invalid(name)
  use phys_module, only: find_RZ_nearby_iter
  character(len=*), intent(in) :: name
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state
  type(particle_kinetic_leapfrog) :: source
  integer :: error
  call mesh()
  if(name=='traversal_limit') then
    call mesh(max_trace_elements+2)
    call fixture(sim,passive,state,0)
    source%i_elm=wall_element-max_trace_elements;source%st=[0.5d0,0.4d0]
    call interp_RZ(nodes,elements,source%i_elm,source%st(1),source%st(2),source%x(1),source%x(2))
    source%x(3)=0.7d0;source%v=[22000.d0,100.d0,0.d0];source%q=0;source%weight=1.d10
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      p(6)%i_elm=0
      if(rank==min(1,tasks-1)) p(1)=source
    end select
    call advance_gk_neutrals(sim,passive,dt,1)
    call MPI_BARRIER(MPI_COMM_WORLD,error)
    call MPI_ABORT(MPI_COMM_WORLD,99,error)
  endif
  if(name=='zero_unresolved') then
    call fixture(sim,passive,state,0)
    source%i_elm=wall_element;source%st=[0.5d0,-1.d-6]
    call interp_RZ(nodes,elements,wall_element,0.5d0,0.d0,source%x(1),source%x(2))
    source%x(3)=0.d0;source%v=[10000.d0,100.d0,0.d0];source%q=0;source%weight=1.d10
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      p(6)%i_elm=0
      if(rank==min(1,tasks-1)) p(1)=source
    end select
    find_RZ_nearby_iter=0
    call empty_locator_index()
    call advance_gk_neutrals(sim,passive,dt,1)
    call MPI_BARRIER(MPI_COMM_WORLD,error) ! Let a failing rank reach its intended rejection first.
    call MPI_ABORT(MPI_COMM_WORLD,99,error)
  endif
  if(index(name,'recovery_')==1) then
    call recovery_mesh(0.2d0)
    call fixture(sim,passive,state,0)
    source%x=[2.482d0,-0.4d0,7.3d0];source%v=[0.d0,8000.d0,0.d0]
    source%q=0;source%weight=6.d10
    call analytic_location(source,0.2d0)
    if(rank==min(1,tasks-1)) then
      if(name=='recovery_connectivity') elements%element(source%i_elm)%neighbours(3)=source%i_elm
      if(name=='recovery_nonfinite') source%v(1)=nan_value
    endif
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      p(1)=source
    end select
    call advance_gk_neutrals(sim,passive,dt,1)
    call MPI_BARRIER(MPI_COMM_WORLD,error) ! Let a failing rank reach its intended rejection first.
    call MPI_ABORT(MPI_COMM_WORLD,99,error)
  endif
  if(index(name,'wall_')==1) then
    select case(name)
    case('wall_no_cross'); error=8
    case('wall_interior'); error=9
    case('wall_corner'); error=10
    case('wall_grazing'); error=11
    case('wall_metadata'); error=12
    case('wall_overflow'); error=16
    case default;call MPI_ABORT(MPI_COMM_WORLD,98,error)
    end select
    call selected_fixture(sim,passive,state,error,3)
    call corrupt_selected(sim,passive,state,error)
    call reflect_newly_lost_gk_neutrals(sim,passive,state,sim%time,dt)
    call MPI_BARRIER(MPI_COMM_WORLD,error) ! Let a failing rank reach its intended rejection first.
    call MPI_ABORT(MPI_COMM_WORLD,99,error)
  endif
  call fixture(sim,passive,state,3)
  if(name=='type') then
    if(rank==min(1,tasks-1)) then
      deallocate(sim%groups(3)%particles)
      allocate(particle_gc_vpar::sim%groups(3)%particles(nslots))
    endif
    call initialize_gk_neutral_reflection(sim,3,.true.,state)
  else
    call advance_gk_neutrals(sim,passive,dt,1)
    call corrupt(sim,passive,state,name)
    call reflect_newly_lost_gk_neutrals(sim,passive,state,sim%time,dt)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,error)
  call MPI_ABORT(MPI_COMM_WORLD,99,error)
end subroutine
end module mod_gk_neutral_reflection_test
