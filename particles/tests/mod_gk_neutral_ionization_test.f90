module mod_gk_neutral_ionization_test
use mpi
use fruit
use, intrinsic :: ieee_arithmetic
use, intrinsic :: iso_fortran_env, only: int64
use data_structure
use mod_particle_sim, only: particle_sim,configure_particle_groups
use mod_particle_types
use mod_gk_neutral_ionization
use mod_particle_slot_pool
use mod_particle_slot_pool_test, only: assert_slot_registry_equal, snapshot_slot_registry
use mod_fields_linear, only: jorek_fields_interp_linear
use mod_projection_helpers_test_tools, only: default_square_grid
use mod_neighbours, only: update_neighbours
use mod_interp, only: interp_RZ
use constants, only: ATOMIC_MASS_UNIT,EL_CHG
implicit none
private
public :: run_ionization_tests,run_ionization_invalid
integer, parameter :: ion_group=1,electron_group=2,neutral_group=3
type(particle_slot_registry) :: registry
integer :: rank,tasks,interior
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
end module
