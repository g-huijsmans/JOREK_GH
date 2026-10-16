module mod_gk_neutrals_mpi_test
use fruit
use mpi
use mod_gk_neutrals
use mod_particle_types
use mod_particle_sim, only: particle_sim, configure_particle_groups
use constants, only: ATOMIC_MASS_UNIT
implicit none
private
public :: run_fruit_gk_neutrals_mpi
integer :: test_rank, test_tasks
integer, parameter :: neutral_group = 3
contains
subroutine run_fruit_gk_neutrals_mpi(rank, n_tasks, ifail)
  integer, intent(in) :: rank, n_tasks
  integer, intent(inout) :: ifail
  test_rank=rank
  test_tasks=n_tasks
  call run_test_case(test_cartesian_line, 'neutral Cartesian straight line')
  call run_test_case(test_mesh_and_restart, 'neutral mesh loss, groups and restart')
end subroutine

subroutine test_cartesian_line()
  type(particle_kinetic_leapfrog) :: p, initial
  real*8 :: x0(3), v0(3), x(3), v(3), exact(3), errx, errv, energy0, energy
  integer :: i
  if (test_rank /= 0) return
  initial%x = [1.2d0, -0.1d0, 0.7d0]
  initial%v = [-0.3d0, 0.2d0, 0.4d0]
  initial%q = 0
  initial%weight = 7.d13
  call copy_particle_kinetic_leapfrog(initial,p)
  call cartesian(p,x0,v0)
  energy0 = 0.5d0*p%weight*2.01410178d0*ATOMIC_MASS_UNIT*dot_product(p%v,p%v)
  errx = 0.d0
  errv = 0.d0
  do i=1,1000
    call push_gk_neutral_ballistic(p,0.01d0)
    call cartesian(p,x,v)
    exact = x0 + i*0.01d0*v0
    errx = max(errx,maxval(abs(x-exact)))
    errv = max(errv,maxval(abs(v-v0)))
  enddo
  energy = 0.5d0*p%weight*2.01410178d0*ATOMIC_MASS_UNIT*dot_product(p%v,p%v)
  call assert_equals(0,int(p%q),'zero charge')
  call assert_equals(initial%weight,p%weight,0.d0,'unchanged weight')
  call assert_equals(0.d0,errx,2.d-12,'Cartesian position')
  call assert_equals(0.d0,errv,2.d-12,'Cartesian velocity')
  call assert_equals(sqrt(exact(1)**2+exact(3)**2),p%x(1),2.d-12,'R')
  call assert_equals(exact(2),p%x(2),2.d-12,'Z')
  call assert_equals(0.d0,sin(p%x(3)-atan2(exact(3),exact(1))),2.d-12,'phi modulo 2pi')
  call assert_equals(norm2(initial%v),norm2(p%v),2.d-12,'speed')
  call assert_equals(1.d0,energy/energy0,2.d-12,'relative energy')
  write(*,'(A,3ES16.8)') 'NEUTRAL analytic max position/velocity/relative-energy error: ', &
       errx,errv,abs(energy/energy0-1.d0)
  p%x = [1.d0,0.d0,0.d0]
  p%v = [-1.d0,0.d0,0.d0]
  call push_gk_neutral_ballistic(p,1.d0)
  call push_gk_neutral_ballistic(p,1.d0)
  call cartesian(p,x,v)
  call assert_equals(-1.d0,x(1),2.d-12,'axis crossing')
end subroutine

subroutine cartesian(p,x,v)
  type(particle_kinetic_leapfrog), intent(in) :: p
  real*8, intent(out) :: x(3),v(3)
  real*8 :: c,s
  c=cos(p%x(3)); s=sin(p%x(3))
  x=[p%x(1)*c,p%x(2),p%x(1)*s]
  v=[p%v(1)*c-p%v(3)*s,p%v(2),p%v(1)*s+p%v(3)*c]
end subroutine

subroutine test_mesh_and_restart()
  use data_structure
  use hdf5, only: H5open_f
  use mod_fields_linear, only: jorek_fields_interp_linear
  use mod_projection_helpers_test_tools, only: default_square_grid
  use mod_neighbours, only: update_neighbours
  use mod_interp, only: interp_RZ
  use mod_particle_io, only: write_simulation_hdf5, read_simulation_hdf5
  use phys_module, only: n_part_groups, part_groups_in_use, part_group_configs
  use mod_particle_group_id, only: matching_part_config_indices
  type(type_node_list), target :: nodes
  type(type_element_list), target :: elements
  type(particle_sim) :: sim, restored, empty
  type(gk_neutral_state) :: state, no_neutrals
  type(particle_kinetic_leapfrog) :: initial, saved, lost, snapshot(3)
  real*8 :: dummy1,dummy2,r,z,rs,rt,zs,zt,x0(3),v0(3),x(3),v(3)
  integer :: ierr,g,j,unit
  character(len=160) :: restart_file
  character(len=80) :: label
  call init_node_list(nodes,n_nodes_max,nodes%n_dof,n_var)
  call default_square_grid(test_rank,test_tasks,3,3,nodes,elements,ierr,0.5d0,2.5d0,-1.d0,1.d0)
  call update_neighbours(nodes,elements,force_rtree_initialize=.true.)
  sim%my_id=test_rank
  sim%n_mpi=test_tasks
  restored%my_id=test_rank
  restored%n_mpi=test_tasks
  allocate(jorek_fields_interp_linear::sim%fields)
  sim%fields%node_list => nodes
  sim%fields%element_list => elements
  allocate(sim%groups(3), empty%groups(2))
  call initialize_gk_neutrals(empty,neutral_group,no_neutrals)
  call advance_gk_neutrals(empty,no_neutrals,0.01d0,3)
  call report_gk_neutral_diagnostics(empty,no_neutrals)
  call assert_equals(2,size(empty%groups),'absent neutral group is a no-op')
  allocate(particle_gc_vpar::sim%groups(1)%particles(1))
  sim%groups(1)%id='ion'
  sim%groups(1)%particles(1)%x=[1.d0,0.d0,0.d0]
  if (owns_marker(1)) sim%groups(1)%particles(1)%i_elm=1
  ! An unrelated full-orbit charged group also uses non: only the supplied index is neutral.
  allocate(particle_kinetic_leapfrog::sim%groups(2)%particles(1))
  sim%groups(2)%id='ele'
  select type(p=>sim%groups(2)%particles)
  type is(particle_kinetic_leapfrog)
    p(1)%x=[1.d0,0.d0,0.d0]
    p(1)%v=[0.1d0,0.2d0,0.3d0]
    p(1)%q=-1
    if (owns_marker(1)) p(1)%i_elm=1
  end select
  g=neutral_group
  sim%groups(g)%coupling_scheme='non'
  sim%groups(g)%mass=2.01410178d0
  sim%groups(g)%n_particles=3.d0*test_tasks
  sim%groups(g)%Z=1
  sim%groups(g)%ad%suffix=''
  write(sim%groups(g)%id,'(A,I2.2)') 'n',g
  allocate(particle_kinetic_leapfrog::sim%groups(g)%particles(3))
  select type(p=>sim%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
    ! Only the owner initializes each physical marker; other slots stay free.
    p%weight=0.d0
    if (owns_marker(1)) then
      p(1)%x=[1.2d0,0.1d0,0.7d0]
      p(1)%v=[0.3d0,0.2d0,0.4d0]
      p(1)%weight=7.d13
      p(1)%i_life=8
      p(1)%t_birth=0.125
    endif
    if (owns_marker(2)) then
      p(2)%x=[2.49d0,0.1d0,0.7d0]
      p(2)%v=[2.d0,0.d0,0.d0]
      p(2)%weight=3.d13
    endif
    p(3)%i_elm=0
    p(3)%q=1 ! inactive entries are deliberately ignored
    do j=1,2
      if (.not. owns_marker(j)) cycle
      call find_RZ(nodes,elements,p(j)%x(1),p(j)%x(2),dummy1,dummy2, &
                   p(j)%i_elm,p(j)%st(1),p(j)%st(2),ierr)
      call assert_equals(0,ierr,'initial location')
    enddo
    call copy_particle_kinetic_leapfrog(p(1),initial)
  end select
  call initialize_gk_neutrals(sim,neutral_group,state)
  call advance_gk_neutrals(sim,state,0.01d0,3)
  sim%time=0.03d0
  call report_gk_neutral_diagnostics(sim,state)
  if (owns_marker(1)) call cartesian(initial,x0,v0)
  select type(p=>sim%groups(neutral_group)%particles)
  type is(particle_kinetic_leapfrog)
    if (owns_marker(1)) then
      call cartesian(p(1),x,v)
      call assert_equals(0.d0,maxval(abs(x-x0-0.03d0*v0)),2.d-12,'three mesh substeps')
      call assert_equals(0.d0,maxval(abs(v-v0)),2.d-12,'mesh Cartesian velocity')
      call assert_true(p(1)%i_elm>0,'inside mesh')
      call assert_true(all(p(1)%st>=0.d0).and.all(p(1)%st<=1.d0),'valid st')
      call interp_RZ(nodes,elements,p(1)%i_elm,p(1)%st(1),p(1)%st(2),r,rs,rt,z,zs,zt)
      call assert_equals(p(1)%x(1),r,1.d-10,'R reconstructed from st')
      call assert_equals(p(1)%x(2),z,1.d-10,'Z reconstructed from st')
      call copy_particle_kinetic_leapfrog(p(1),saved)
    endif
    if (owns_marker(2)) call assert_true(p(2)%i_elm<0,'negative boundary element')
    call copy_particle_kinetic_leapfrog(p(2),lost)
  end select
  call assert_equals(1.d0,sim%groups(1)%particles(1)%x(1),0.d0,'charged group untouched')
  call check_global_fixture(sim)
  ! Poison fields, leaving the independently stored mesh geometry untouched.
  do j=1,nodes%n_nodes
    nodes%node(j)%values=huge(1.d0)/100.d0
    nodes%node(j)%deltas=-huge(1.d0)/100.d0
  enddo
  select type(p=>sim%groups(neutral_group)%particles)
  type is(particle_kinetic_leapfrog)
    if (owns_marker(1)) call copy_particle_kinetic_leapfrog(initial,p(1))
  end select
  call advance_gk_neutrals(sim,state,0.01d0,3)
  call report_gk_neutral_diagnostics(sim,state)
  select type(p=>sim%groups(neutral_group)%particles)
  type is(particle_kinetic_leapfrog)
    if (owns_marker(1)) &
      call assert_equals(0.d0,maxval(abs(p(1)%x-saved%x)),0.d0,'independent of fields')
  end select
  select type(p=>sim%groups(neutral_group)%particles)
  type is(particle_kinetic_leapfrog)
    if (owns_marker(2)) &
      call assert_equals(0.d0,maxval(abs(p(2)%x-lost%x)),0.d0,'lost marker stops advancing')
    call copy_particle_kinetic_leapfrog(p(1),saved)
  end select
  call assert_equals(1.d0,sim%groups(2)%particles(1)%x(1),0.d0,'unrelated full-orbit group x')
  select type(p=>sim%groups(2)%particles)
  type is(particle_kinetic_leapfrog)
    call assert_equals(0.d0,maxval(abs(p(1)%v-[0.1d0,0.2d0,0.3d0])),0.d0,'unrelated group velocity')
    call assert_equals(-1,int(p(1)%q),'unrelated group charge')
  end select
  n_part_groups=3
  part_groups_in_use(1:3)=sim%groups%id
  g=neutral_group
  select type(p=>sim%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,size(p)
      call copy_particle_kinetic_leapfrog(p(j),snapshot(j))
    enddo
  end select
  call unique_restart_filename(restart_file)
  ! HDF5 Fortran constants must be initialized before the writer selects its access mode.
  call H5open_f(ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call write_simulation_hdf5(sim,trim(restart_file))
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! root has closed the gathered file before any rank reads

  ! The common reader expects groups configured from the input before import.
  do g=1,3
    matching_part_config_indices(g)=g
    part_group_configs(g)%id=sim%groups(g)%id
    part_group_configs(g)%Z=sim%groups(g)%Z
    part_group_configs(g)%mass=sim%groups(g)%mass
    part_group_configs(g)%coupling_scheme=sim%groups(g)%coupling_scheme
    part_group_configs(g)%n_particles=size(sim%groups(g)%particles)*test_tasks
    part_group_configs(g)%type='particle_kinetic_leapfrog'
  enddo
  part_group_configs(1)%type='particle_gc_vpar'
  call restored%allocate_groups(3)
  call configure_particle_groups(restored)
  call read_simulation_hdf5(restored,trim(restart_file),test_in=.true.)
  g=neutral_group
  call assert_equals(sim%groups(g)%id,restored%groups(g)%id,'restart group ID')
  select type(p=>restored%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
    call assert_equals(3,size(p),'restart local slots')
    do j=1,size(p)
      call copy_particle_kinetic_leapfrog(snapshot(j),saved)
      write(label,'(A,I0,A,I0,A,I0)') ' rank=',test_rank,' group=',g,' slot=',j
      call assert_equals(0.d0,maxval(abs(p(j)%x-saved%x)),0.d0,'restart x'//trim(label))
      call assert_equals(0.d0,maxval(abs(p(j)%v-saved%v)),0.d0,'restart v'//trim(label))
      call assert_equals(0.d0,maxval(abs(p(j)%st-saved%st)),0.d0,'restart st'//trim(label))
      call assert_equals(saved%weight,p(j)%weight,0.d0,'restart weight'//trim(label))
      call assert_equals(int(saved%q),int(p(j)%q),'restart q'//trim(label))
      call assert_equals(saved%i_elm,p(j)%i_elm,'restart i_elm'//trim(label))
      call assert_equals(saved%i_life,p(j)%i_life,'restart i_life'//trim(label))
      call assert_equals(dble(saved%t_birth),dble(p(j)%t_birth),0.d0,'restart t_birth'//trim(label))
    enddo
  class default
    call assert_true(.false.,'restart particle type')
  end select
  call check_global_fixture(restored)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! all readers have closed the file
  if (test_rank == 0) then
    open(newunit=unit,file=trim(restart_file),status='old',iostat=ierr)
    call assert_equals(0,ierr,'open temporary restart for cleanup')
    if (ierr == 0) then
      close(unit,status='delete',iostat=ierr)
      call assert_equals(0,ierr,'delete temporary restart')
    endif
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
end subroutine

logical function owns_marker(global_index)
  integer, intent(in) :: global_index
  owns_marker = test_rank == mod(global_index-1,test_tasks)
end function

! Independent particle-data reference: the physical population never scales with MPI size.
subroutine check_global_fixture(sim)
  type(particle_sim), intent(in) :: sim
  integer :: g,j,ierr,counts(2),global_counts(2)
  real*8 :: totals(4),global_totals(4),mass
  g=neutral_group
  counts=0
  totals=0.d0
  mass=2.01410178d0*ATOMIC_MASS_UNIT
  select type(p=>sim%groups(g)%particles)
  type is(particle_kinetic_leapfrog)
    do j=1,size(p)
      if(p(j)%i_elm>0) then
        counts(1)=counts(1)+1
        totals(1)=totals(1)+p(j)%weight
        totals(2)=totals(2)+0.5d0*mass*p(j)%weight*dot_product(p(j)%v,p(j)%v)
      else if(p(j)%i_elm<0) then
        counts(2)=counts(2)+1
        totals(3)=totals(3)+p(j)%weight
        totals(4)=totals(4)+0.5d0*mass*p(j)%weight*dot_product(p(j)%v,p(j)%v)
      endif
    enddo
  end select
  call MPI_ALLREDUCE(counts,global_counts,2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(totals,global_totals,4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(test_rank==0) then
    call assert_equals(1,global_counts(1),'global active markers')
    call assert_equals(1,global_counts(2),'global lost markers')
    call assert_equals(7.d13,global_totals(1),0.d0,'global active weight')
    call assert_equals(3.d13,global_totals(3),0.d0,'global lost weight')
    call assert_equals(1.d0,global_totals(2)/(0.5d0*mass*7.d13*0.29d0),2.d-12,'global active energy')
    call assert_equals(1.d0,global_totals(4)/(0.5d0*mass*3.d13*4.d0),2.d-12,'global lost energy')
  endif
end subroutine

subroutine unique_restart_filename(filename)
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
  if(test_rank==0) then
    do
      call system_clock(count=tick)
      write(filename,'(A,I0,A,I0,A)') 'gk_neutrals_roundtrip_',fgetpid(),'_',tick,'.h5'
      inquire(file=trim(filename),exist=exists)
      if(.not.exists) exit
    enddo
  endif
  call MPI_BCAST(filename,len(filename),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
end subroutine
end module mod_gk_neutrals_mpi_test
