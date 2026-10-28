module mod_gk_neutral_reflection_test
use mpi
use fruit
use data_structure
use, intrinsic :: ieee_arithmetic
use mod_particle_sim, only: particle_sim,configure_particle_groups
use mod_particle_types
use mod_gk_neutrals
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
subroutine mesh()
  integer :: error,e,s,n1,n2,k
  logical :: sides(4)
  real*8 :: r,z
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,tasks,error)
  nan_value=ieee_value(0.d0,ieee_quiet_nan)
  call init_node_list(nodes,n_nodes_max,nodes%n_dof,n_var)
  call default_square_grid(rank,tasks,5,5,nodes,elements,error,0.5d0,2.5d0,-1.d0,1.d0)
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
  integer :: count_before,valid
  allocate(ions,source=sim%groups(1)%particles); allocate(electrons,source=sim%groups(2)%particles)
  allocate(neutrals,source=sim%groups(3)%particles)
  call initialize_particle_slot_registry(sim,[1,2,3],pool)
  call snapshot_slot_registry(pool,before)
  old=state; events=passive%newly_lost_indices; count_before=passive%n_newly_lost_indices
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
  old=state
  call checked_reflect(sim,passive,state,error,.false.)
  call same_state(state,old,'empty events no-op')
end subroutine
subroutine reflections()
  type(particle_sim) :: sim,replay
  type(gk_neutral_state) :: passive,again
  type(gk_neutral_reflection_state) :: state,other,before
  type(pcg32_rng) :: expected_rng
  class(particle_base), allocatable :: hit(:),saved(:)
  real*8 :: r,z,normal(3),u(2),expected(3),expected_energy
  integer :: j,error,n
  call fixture(sim,passive,state,3)
  call fixture(replay,again,other,3)
  call advance_gk_neutrals(sim,passive,dt,3)
  call advance_gk_neutrals(replay,again,dt,3)
  allocate(hit,source=sim%groups(3)%particles)
  n=count([(owns(j),j=1,3)])
  call assert_equals(n,passive%n_newly_lost_indices,'one event per physical current loss')
  call assert_true(all(passive%newly_lost_indices(:n)<=3),'free/old-negative/interior protected')
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
subroutine run_reflection_tests()
  call mesh()
  call run_test_case(noops_and_absorption,'disabled/absent/empty neutral reflection')
  call run_test_case(reflections,'diffuse identity, energy, replay, repeated impacts and restart')
  call run_test_case(grazing_and_zero,'grazing tolerance and zero-speed RNG')
  call run_test_case(multiple_side_and_initialization,'boundary-side resolution and initialization validation')
  call run_test_case(rejections,'global rollback of every particle/pool/RNG/diagnostic/event')
end subroutine
subroutine run_reflection_invalid(name)
  character(len=*), intent(in) :: name
  type(particle_sim) :: sim
  type(gk_neutral_state) :: passive
  type(gk_neutral_reflection_state) :: state
  integer :: error
  call mesh()
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
  call MPI_ABORT(MPI_COMM_WORLD,99,error)
end subroutine
end module mod_gk_neutral_reflection_test
