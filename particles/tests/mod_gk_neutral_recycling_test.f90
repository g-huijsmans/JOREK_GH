module mod_gk_neutral_recycling_test
use mpi
use fruit
use data_structure
use mod_particle_sim, only: particle_sim, configure_particle_groups
use mod_particle_types
use mod_gk_neutral_recycling
use mod_particle_slot_pool
use mod_particle_slot_pool_test, only: assert_slot_registry_equal, snapshot_slot_registry
use mod_gk_neutrals
use mod_pcg32_rng, only: pcg32_rng
use mod_fields_linear, only: jorek_fields_interp_linear
use mod_projection_helpers_test_tools, only: default_square_grid
use mod_neighbours, only: update_neighbours
use mod_boundary, only: wall_normal_vector
use mod_interp, only: interp_RZ
use constants, only: ATOMIC_MASS_UNIT, EL_CHG
implicit none
private
public :: run_recycling_tests, run_recycling_invalid
integer, parameter :: ion_group=1, neutral_group=3
type(particle_slot_registry) :: registry
integer :: rank,tasks,boundary_element,interior_element,corner_element,boundary_side,corner_side
real*8 :: wall_st(2),corner_st(2)
type(type_node_list), target :: nodes
type(type_element_list), target :: elements
contains

subroutine checked_report(sim,state)
  type(particle_sim), intent(in) :: sim
  type(gk_neutral_recycling_state), intent(inout) :: state
  type(particle_slot_registry) :: before
  call snapshot_slot_registry(registry,before)
  call report_gk_neutral_recycling(sim,state)
  call assert_slot_registry_equal(registry,before,'report leaves registry unchanged')
end subroutine


subroutine checked_recycle(sim,state,time,dt,ierr)
  type(particle_sim), intent(inout) :: sim
  type(gk_neutral_recycling_state), intent(inout) :: state
  real*8, intent(in) :: time,dt
  integer, optional, intent(out) :: ierr
  type(particle_slot_registry) :: before
  type(gk_neutral_recycling_state) :: old
  type(pcg32_rng) :: rng_before,rng_after
  real*8 :: u_before(2),u_after(2)
  integer :: error,g,n,old_free(3)
  if(.not.state%enabled) then
    call recycle_lost_gk_ions(sim,state,registry,time,dt,ierr)
    return
  endif
  call validate_particle_slot_registry(sim,registry,[1,2,3],error)
  call assert_equals(0,error,'pool agrees before recycling')
  call snapshot_slot_registry(registry,before); old=state; old_free=registry%pool%n_free
  n=count(sim%groups(1)%particles%i_elm<0)
  call recycle_lost_gk_ions(sim,state,registry,time,dt,error)
  if(error/=0) then
    call assert_slot_registry_equal(registry,before,'recycling rollback')
    rng_before=old%rng; rng_after=state%rng
    call rng_before%next(u_before); call rng_after%next(u_after)
    call assert_true(all(u_before==u_after),'recycling failure preserves RNG')
    call assert_true(all(old%totals==state%totals).and.all(old%counts(:4)==state%counts(:4)), &
                     'recycling failure preserves diagnostics except failed attempts')
  else
    call assert_true(all(registry%pool%n_free==old_free+[n,0,-n]),'recycling pool transfers')
  endif
  call validate_particle_slot_registry(sim,registry,[1,2,3],g)
  call assert_equals(0,g,'pool agrees after recycling')
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

subroutine setup_mesh()
  integer :: ierr,e,s,nside,n1,n2,side
  logical :: sides(4)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,tasks,ierr)
  call init_node_list(nodes,n_nodes_max,nodes%n_dof,n_var)
  call default_square_grid(rank,tasks,5,5,nodes,elements,ierr,0.5d0,2.5d0,-1.d0,1.d0)
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
  boundary_element=0; interior_element=0; corner_element=0
  do e=1,elements%n_elements
    do s=1,4
      n1=elements%element(e)%vertex(s); n2=elements%element(e)%vertex(mod(s,4)+1)
      sides(s)=nodes%node(n1)%boundary>0.and.nodes%node(n2)%boundary>0
    enddo
    nside=count(sides)
    if(nside==0) interior_element=e
    side=findloc(sides,.true.,dim=1)
    if(nside==1) then
      boundary_element=e; boundary_side=side; wall_st=side_st(side)
    endif
    if(nside==2) then
      corner_element=e; corner_side=side; corner_st=side_st(side)
    endif
  enddo
end subroutine

function side_st(side) result(st)
  integer, intent(in) :: side
  real*8 :: st(2)
  select case(side)
  case(1); st=[0.4d0,0.d0]
  case(2); st=[1.d0,0.4d0]
  case(3); st=[0.4d0,1.d0]
  case(4); st=[0.d0,0.4d0]
  end select
end function

logical function owns(j)
  integer, intent(in) :: j
  owns=rank==mod(j-1,tasks)
end function

subroutine fixture(sim,number)
  use phys_module, only: use_manual_random_seed,manual_seed
  type(particle_sim), intent(out) :: sim
  integer, intent(in) :: number
  integer :: j
  use_manual_random_seed=.true.; manual_seed=12345
  sim%my_id=rank; sim%n_mpi=tasks; sim%time=0.125d0
  allocate(jorek_fields_interp_linear::sim%fields)
  sim%fields%node_list=>nodes; sim%fields%element_list=>elements
  allocate(sim%groups(3))
  sim%groups%mass=2.01410178d0
  sim%groups%coupling_scheme='non'
  sim%groups%ad%suffix=''
  sim%groups(1)%id='ion'; sim%groups(2)%id='ele'; sim%groups(3)%id='neu'
  allocate(particle_gc_vpar::sim%groups(1)%particles(3),sim%groups(2)%particles(1))
  allocate(particle_kinetic_leapfrog::sim%groups(3)%particles(5))
  select type(p=>sim%groups(1)%particles)
  type is(particle_gc_vpar)
    do j=1,number
      if(.not.owns(j)) cycle
      p(j)%i_elm=-boundary_element; p(j)%st=wall_st
      p(j)%x=[9.d0,9.d0,0.2d0*j] ! deliberately not the wall R,Z
      p(j)%q=1; p(j)%weight=1.d10*j
      p(j)%vpar=100.d0*j; p(j)%mu=20.d0*j; p(j)%B_norm=2.d0
    enddo
  end select
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    p%weight=0.d0
    p(1)%i_elm=interior_element; p(1)%x=[1.5d0,0.d0,0.1d0]; p(1)%st=0.5d0
    p(2)%i_elm=-boundary_element; p(2)%x=[8.d0,8.d0,0.2d0]
    p(3:5)%i_life=7
    p(3:5)%q=1 ! inactive charge is overwritten on creation
  end select
  call initialize_particle_slot_registry(sim,[1,2,3],registry)
end subroutine

subroutine run_recycling_tests()
  call setup_mesh()
  call run_test_case(single_recycle,'single product, inward motion, replay and no double recycling')
  call run_test_case(multiple_restart,'multiple unequal markers and MPI-safe restart')
  call run_test_case(transaction_checks,'capacity, adjacent boundary nodes, corner and invalid wall state')
  call run_test_case(cross_rank_rollback,'global recycling rollback of particles, pools and RNG')
end subroutine

subroutine check_neutral(a,b,label)
  type(particle_kinetic_leapfrog), intent(in) :: a,b
  character(len=*), intent(in) :: label
  call assert_equals(0.d0,maxval(abs(a%x-b%x)),0.d0,label//' x')
  call assert_equals(0.d0,maxval(abs(a%v-b%v)),0.d0,label//' v')
  call assert_equals(0.d0,maxval(abs(a%st-b%st)),0.d0,label//' st')
  call assert_equals(a%i_elm,b%i_elm,label//' i_elm')
  call assert_equals(int(a%q),int(b%q),label//' q')
  call assert_equals(a%weight,b%weight,0.d0,label//' weight')
  call assert_equals(a%i_life,b%i_life,label//' life')
  call assert_equals(dble(a%t_birth),dble(b%t_birth),0.d0,label//' birth')
end subroutine

subroutine single_recycle()
  type(particle_sim) :: sim,replay,empty
  type(gk_neutral_recycling_state) :: state,again,disabled
  type(gk_neutral_state) :: passive
  type(particle_kinetic_leapfrog) :: saved(5)
  real*8 :: r,z,normal(3),energy
  integer :: j
  allocate(empty%groups(2))
  call initialize_gk_neutral_recycling(empty,ion_group,neutral_group,disabled)
  call checked_recycle(empty,disabled,0.d0,1.d-7)
  call checked_report(empty,disabled)
  call assert_true(.not.disabled%enabled,'charged-only recycling no-op')
  call fixture(sim,1)
  call initialize_gk_neutral_recycling(sim,ion_group,neutral_group,state)
  call checked_recycle(sim,state,0.25d0,1.d-7)
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    saved=p
    if(owns(1)) then
      call interp_RZ(nodes,elements,boundary_element,wall_st(1),wall_st(2),r,z)
      normal=wall_normal_vector(nodes,elements,boundary_element,wall_st(1),wall_st(2))
      call assert_equals(boundary_element,p(3)%i_elm,'positive wall element')
      call assert_equals(0.d0,maxval(abs(p(3)%st-wall_st)),0.d0,'retained wall st')
      call assert_equals(0.d0,maxval(abs(p(3)%x-[r,z,0.2d0])),0.d0,'wall interpolation and phi')
      call assert_equals(0,int(p(3)%q),'neutral charge')
      call assert_equals(1.d10,p(3)%weight,0.d0,'incident weight')
      call assert_true(dot_product(p(3)%v,normal)>0.d0,'inward velocity')
      energy=sum(p(3)%v**2)/(100.d0**2+80.d0)
      call assert_equals(1.d0,energy,2.d-14,'incident kinetic energy')
      call assert_equals(8,p(3)%i_life,'free-slot lifetime increment')
      call assert_equals(0.25d0,dble(p(3)%t_birth),0.d0,'explicit physical birth time')
      call assert_equals(0,sim%groups(1)%particles(1)%i_elm,'source inactive')
    endif
  end select
  call checked_report(sim,state)
  if(rank==0) call assert_equals(1,state%global_counts(1),'one global source')
  call checked_recycle(sim,state,0.3d0,1.d-7)
  call assert_true(all(state%counts==0).and.all(state%totals==0.d0),'no double recycling increments')
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,5
      call check_neutral(saved(j),p(j),'second call unchanged')
    enddo
  end select
  call checked_report(sim,state)
  call fixture(replay,1)
  call initialize_gk_neutral_recycling(replay,ion_group,neutral_group,again)
  call checked_recycle(replay,again,0.25d0,1.d-7)
  select type(p=>replay%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,5
      call check_neutral(saved(j),p(j),'fixed-seed replay')
    enddo
    ! Remove the unrelated active sentinel before testing boundary launch.
    p(1)%i_elm=0
  end select
  call initialize_gk_neutrals(replay,neutral_group,passive)
  call advance_gk_neutrals(replay,passive,1.d-7,1)
  if(owns(1)) then
    call assert_true(replay%groups(3)%particles(3)%i_elm>0,'boundary product remains inside on next advance')
    call assert_true(any(replay%groups(3)%particles(3)%x/=saved(3)%x),'product moves next step')
  endif
end subroutine

subroutine multiple_restart()
  use hdf5, only: H5open_f
  use mod_particle_io, only: write_simulation_hdf5,read_simulation_hdf5
  use phys_module, only: n_part_groups,part_groups_in_use,part_group_configs
  use mod_particle_group_id, only: matching_part_config_indices
  type(particle_sim) :: sim,restored
  type(gk_neutral_recycling_state) :: state
  type(particle_kinetic_leapfrog) :: saved(5)
  real*8 :: expected_energy,weight,energy
  integer :: j,g,ierr,unit
  character(len=160) :: filename
  call fixture(sim,3)
  call initialize_gk_neutral_recycling(sim,ion_group,neutral_group,state)
  call checked_recycle(sim,state,0.25d0,1.d-7)
  call checked_report(sim,state)
  expected_energy=0.d0
  do j=1,3
    expected_energy=expected_energy+0.5d0*2.01410178d0*ATOMIC_MASS_UNIT*1.d10*j*(100.d0**2*j*j+80.d0*j)
  enddo
  if(rank==0) then
    call assert_equals(3,state%global_counts(1),'global incident count')
    call assert_equals(3,state%global_counts(2),'global created count')
    call assert_equals(6.d10,state%global_totals(1),0.d0,'global incident weight')
    call assert_equals(6.d10,state%global_totals(3),0.d0,'global created weight')
    call assert_equals(1.d0,state%global_totals(2)/expected_energy,2.d-14,'global incident energy')
    call assert_equals(1.d0,state%global_totals(4)/expected_energy,2.d-14,'global created energy')
    call assert_equals(1.d0,state%global_totals(5)/(EL_CHG*6.d10),2.d-14,'wall charge')
  endif
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    saved=p
    call assert_equals(interior_element,p(1)%i_elm,'active slot protected')
    call assert_equals(-boundary_element,p(2)%i_elm,'negative slot protected')
    call assert_equals(8.d0,p(2)%x(1),0.d0,'negative slot data protected')
  end select
  n_part_groups=3; part_groups_in_use(1:3)=sim%groups%id
  call unique_filename(filename)
  call H5open_f(ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call write_simulation_hdf5(sim,trim(filename))
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  n_part_groups=3; part_groups_in_use(1:3)=sim%groups%id
  do g=1,3
    matching_part_config_indices(g)=g
    part_group_configs(g)%id=sim%groups(g)%id
    part_group_configs(g)%mass=sim%groups(g)%mass
    part_group_configs(g)%Z=1
    part_group_configs(g)%coupling_scheme='non'
    part_group_configs(g)%n_particles=size(sim%groups(g)%particles)*tasks
    part_group_configs(g)%type='particle_gc_vpar'
  enddo
  part_group_configs(3)%type='particle_kinetic_leapfrog'
  restored%my_id=rank; restored%n_mpi=tasks
  call restored%allocate_groups(3)
  call configure_particle_groups(restored)
  call read_simulation_hdf5(restored,trim(filename),test_in=.true.)
  call assert_equals('neu',restored%groups(3)%id,'restart group ID')
  call check_restart_pool(restored)
  select type(p=>restored%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,5
      call check_neutral(saved(j),p(j),'restart')
    enddo
  class default
    call assert_true(.false.,'restart dynamic type')
  end select
  call assert_true(all(restored%groups(1)%particles%i_elm==0),'processed ions remain inactive after restart')
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if(rank==0) then
    open(newunit=unit,file=trim(filename),status='old',iostat=ierr)
    call assert_equals(0,ierr,'open restart for cleanup')
    if(ierr==0) close(unit,status='delete')
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine

subroutine transaction_checks()
  type(particle_sim) :: sim
  type(gk_neutral_recycling_state) :: state
  type(particle_kinetic_leapfrog) :: saved(5)
  integer :: ierr,j,k,n1,n3,old1,old3
  call assert_true(interior_element>0.and.boundary_element>0.and.corner_element>0,'mesh boundary fixtures')
  do k=1,5
    call fixture(sim,3)
    call initialize_gk_neutral_recycling(sim,ion_group,neutral_group,state)
    if(k==1) sim%groups(3)%particles(3:5)%i_elm=-boundary_element
    if(k==2) then
      do j=1,3
        if(owns(j)) sim%groups(1)%particles(j)%i_elm=-interior_element
      enddo
    endif
    if(k==3) then
      ! Two opposite nodes are boundary-marked; there is still no boundary side.
      n1=elements%element(interior_element)%vertex(1); n3=elements%element(interior_element)%vertex(3)
      old1=nodes%node(n1)%boundary; old3=nodes%node(n3)%boundary
      nodes%node(n1)%boundary=1; nodes%node(n3)%boundary=1
      do j=1,3
        if(owns(j)) sim%groups(1)%particles(j)%i_elm=-interior_element
      enddo
    endif
    if(k==4) then
      do j=1,3
        if(.not.owns(j)) cycle
        sim%groups(1)%particles(j)%i_elm=-corner_element
        sim%groups(1)%particles(j)%st=corner_st
      enddo
    endif
    if(k==5) then
      do j=1,3
        if(owns(j)) sim%groups(1)%particles(j)%st=0.5d0
      enddo
    endif
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      saved=p
    end select
    call initialize_particle_slot_registry(sim,[1,2,3],registry) ! finalized test fixture
    call checked_recycle(sim,state,0.25d0,1.d-7,ierr)
    if(k==3) then
      nodes%node(n1)%boundary=old1; nodes%node(n3)%boundary=old3
    endif
    if(k==4) then
      call assert_equals(0,ierr,'unambiguous multiple-side element accepted')
      call assert_equals(count([(owns(j),j=1,3)]),state%counts(4),'multiple-side diagnostic')
    else
      call assert_true(ierr/=0,'invalid transaction rejected')
      call assert_equals(count([(owns(j),j=1,3)]),count(sim%groups(1)%particles%i_elm<0),'sources unchanged on failure')
      select type(p=>sim%groups(3)%particles)
      type is(particle_kinetic_leapfrog)
        do j=1,5
          call check_neutral(saved(j),p(j),'transaction rollback')
        enddo
      end select
    endif
  enddo
end subroutine

subroutine cross_rank_rollback()
  type(particle_sim) :: sim
  type(gk_neutral_recycling_state) :: state
  type(particle_gc_vpar) :: ions(3),electron(1)
  type(particle_kinetic_leapfrog) :: neutrals(5)
  integer :: j,error
  call fixture(sim,3)
  call initialize_gk_neutral_recycling(sim,ion_group,neutral_group,state)
  if(rank==min(1,tasks-1)) sim%groups(3)%particles(3:5)%i_elm=-boundary_element
  call initialize_particle_slot_registry(sim,[1,2,3],registry)
  select type(p=>sim%groups(1)%particles)
  type is(particle_gc_vpar)
    ions=p
  end select
  select type(p=>sim%groups(2)%particles)
  type is(particle_gc_vpar)
    electron=p
  end select
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    neutrals=p
  end select
  call checked_recycle(sim,state,0.25d0,1.d-7,error)
  call assert_true(error/=0,'rank-one failure rejects on every rank')
  select type(p=>sim%groups(1)%particles)
  type is(particle_gc_vpar)
    do j=1,size(p)
      call check_gc(p(j),ions(j))
    enddo
  end select
  select type(p=>sim%groups(2)%particles)
  type is(particle_gc_vpar)
    call check_gc(p(1),electron(1))
  end select
  select type(p=>sim%groups(3)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,size(p)
      call check_neutral(p(j),neutrals(j),'cross-rank rollback')
    enddo
  end select
end subroutine

subroutine check_gc(a,b)
  type(particle_gc_vpar), intent(in) :: a,b
  call assert_true(all(a%x==b%x).and.all(a%st==b%st).and.a%i_elm==b%i_elm.and. &
       a%q==b%q.and.a%weight==b%weight.and.a%vpar==b%vpar.and.a%mu==b%mu.and. &
       a%B_norm==b%B_norm.and.a%i_life==b%i_life.and.a%t_birth==b%t_birth,'complete GC rollback')
end subroutine

subroutine run_recycling_invalid(mode)
  character(len=*), intent(in) :: mode
  type(particle_sim) :: sim
  type(gk_neutral_recycling_state) :: state
  integer :: j
  call setup_mesh()
  call fixture(sim,3)
  select case(trim(mode))
  case('ion_type')
    deallocate(sim%groups(1)%particles)
    allocate(particle_kinetic_leapfrog::sim%groups(1)%particles(1))
  case('neutral_type')
    deallocate(sim%groups(3)%particles)
    allocate(particle_gc_vpar::sim%groups(3)%particles(1))
  case('ion_charge')
    select type(p=>sim%groups(1)%particles)
    type is(particle_gc_vpar)
      p%q=2
    end select
  case('neutral_charge')
    select type(p=>sim%groups(3)%particles)
    type is(particle_kinetic_leapfrog)
      p(1)%q=1
    end select
  case('mass')
    sim%groups(3)%mass=1.d0
  case('capacity')
    sim%groups(3)%particles(3:5)%i_elm=-boundary_element
  case('interior')
    do j=1,3
      if(owns(j)) sim%groups(1)%particles(j)%i_elm=-interior_element
    enddo
  case('energy')
    select type(p=>sim%groups(1)%particles)
    type is(particle_gc_vpar)
      p%mu=-1.d20
    end select
  case('corner')
    do j=1,3
      if(.not.owns(j)) cycle
      sim%groups(1)%particles(j)%i_elm=-corner_element
      ! This top-right corner fixture has sides 2 and 3.
      sim%groups(1)%particles(j)%st=1.d0
    enddo
  end select
  call initialize_gk_neutral_recycling(sim,ion_group,neutral_group,state)
  call initialize_particle_slot_registry(sim,[1,2,3],registry)
  call recycle_lost_gk_ions(sim,state,registry,0.25d0,1.d-7)
end subroutine

subroutine unique_filename(filename)
  use, intrinsic :: iso_c_binding, only: c_int
  character(len=*), intent(out) :: filename
  integer :: ierr
  integer*8 :: tick
  logical :: exists
  interface
    function fgetpid() bind(C)
      import c_int
      integer(c_int) :: fgetpid
    end function
  end interface
  if(rank==0) then
    do
      call system_clock(count=tick)
      write(filename,'(A,I0,A,I0,A)') 'recycling_',fgetpid(),'_',tick,'.h5'
      inquire(file=trim(filename),exist=exists)
      if(.not.exists) exit
    enddo
  endif
  call MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
end subroutine
end module mod_gk_neutral_recycling_test
