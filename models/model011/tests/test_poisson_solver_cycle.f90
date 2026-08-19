!> Developer regression for the persistent model011 harmonic Poisson solve.
program test_poisson_solver_cycle
  use mpi_mod
  use data_structure,      only: type_RHS, init_node_list, nbthreads
  use mod_parameters,      only: n_tor, n_degrees, n_eq_var, var_psi, var_u, var_rho, var_T
  use phys_module,         only: use_mumps, use_pastix, use_strumpack, F0, central_mass
  use basis_at_gaussian,   only: initialise_basis
  use nodes_elements,      only: node_list, element_list, bnd_node_list, bnd_elm_list, aux_node_list
  use equil_info,          only: ES
  use mod_simulation_data, only: type_MHD_SIM
  use mod_poisson_solver,  only: poisson_solve_action
  use tr_module,           only: tr_meminit

  implicit none

  type(type_MHD_SIM), target :: mhd_sim
  type(poisson_solve_action) :: poisson
  type(type_RHS)             :: rhs
  real*8, allocatable        :: phi_first(:), phi_second(:)
  real*8, allocatable        :: phi_old(:,:,:), psi_before(:,:,:), rho_before(:,:,:)
  real*8, allocatable        :: temperature_before(:,:,:), rhs_before_store(:)
  integer                    :: rank, n_tasks, ierr, inode, idof, index, entry
  integer                    :: row_component, column_component
  real*8                     :: changed, repeated, scale_error, cross_harmonic, cross_harmonic_global
  real*8                     :: store_error, addition_separation, delta_error, other_field_error

  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, ierr)
  if (n_tasks.ne.2) then
    if (rank.eq.0) write(*,*) 'This test requires exactly two MPI ranks.'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  call tr_meminit(rank, n_tasks)
  call preset_parameters()
  call det_modes()
  call initialise_basis()
  nbthreads = 1

  ! Exercise MUMPS explicitly; the production selector remains namelist-driven.
  use_mumps = .true.
  use_pastix = .false.
  use_strumpack = .false.
  F0 = 10.d0
  central_mass = 2.d0

  call init_node_list(node_list, 4, 0, n_eq_var)
  node_list%n_nodes = 4
  node_list%n_values = n_eq_var
  element_list%n_elements = 1
  element_list%element(1)%vertex = (/ 1, 2, 3, 4 /)
  element_list%element(1)%neighbours = 0
  element_list%element(1)%size = 1.d0
  element_list%element(1)%father = 0
  element_list%element(1)%n_sons = 0
  element_list%element(1)%n_gen = 0
  element_list%element(1)%sons = 0
  element_list%element(1)%contain_node = 0
  element_list%element(1)%nref = 0
  bnd_node_list%n_bnd_nodes = 0
  bnd_elm_list%n_bnd_elements = 0
  nullify(aux_node_list)

  ! Unit square in bicubic Hermite data: value, d/ds, d/dt, d2/dsdt.
  node_list%node(1)%x = 0.d0
  node_list%node(2)%x = 0.d0
  node_list%node(3)%x = 0.d0
  node_list%node(4)%x = 0.d0
  node_list%node(1)%x(1,1,:) = (/ 1.d0, 0.d0 /)
  node_list%node(2)%x(1,1,:) = (/ 2.d0, 0.d0 /)
  node_list%node(3)%x(1,1,:) = (/ 2.d0, 1.d0 /)
  node_list%node(4)%x(1,1,:) = (/ 1.d0, 1.d0 /)
  node_list%node(:)%boundary = 1
  node_list%node(:)%boundary_index = 0
  node_list%node(:)%axis_node = .false.
  node_list%node(:)%axis_dof = 0
  node_list%node(:)%parent_elem = 0
  node_list%node(:)%ref_lambda = 0.d0
  node_list%node(:)%ref_mu = 0.d0
  node_list%node(:)%constrained = .false.

  index = 0
  do inode = 1, 4
    node_list%node(inode)%parents = 0
    node_list%node(inode)%values = 0.d0
    node_list%node(inode)%deltas = 0.d0
    node_list%node(inode)%values(1,1,var_psi) = 0.d0
    node_list%node(inode)%values(1,1,var_rho) = 1.d0
    node_list%node(inode)%values(:,:,var_u) = 10.d0 + real(inode,8)
    node_list%node(inode)%deltas(:,:,var_u) = -3.d0-real(inode,8)
    node_list%node(inode)%values(:,:,var_T) = 20.d0 + real(inode,8)
    do idof = 1, n_degrees
      index = index + 1
      node_list%node(inode)%index(idof) = index
    enddo
  enddo

  ES%xpoint = .false.
  ES%xcase = 0
  ES%R_axis = 1.5d0
  ES%Z_axis = 0.5d0
  ES%psi_axis = 0.d0
  ES%psi_bnd = 1.d0
  ES%R_xpoint = 0.d0
  ES%Z_xpoint = 0.d0
  ES%psi_xpoint = 0.d0

  mhd_sim%my_id = rank
  mhd_sim%n_mpi = n_tasks
  mhd_sim%n_tor = n_tor
  mhd_sim%freeboundary = .false.
  mhd_sim%restart = .false.
  mhd_sim%sr_n_tor = 0
  mhd_sim%node_list => node_list
  mhd_sim%element_list => element_list
  mhd_sim%bnd_node_list => bnd_node_list
  mhd_sim%bnd_elm_list => bnd_elm_list
  mhd_sim%es => ES

  call poisson%setup(mhd_sim, MPI_COMM_WORLD)
  call poisson%construct_matrix()
  if (.not.poisson%matrix_ready .or. poisson%factorized) &
    error stop 'Invalid state after Poisson matrix construction.'

  cross_harmonic = 0.d0
  do entry = 1, poisson%solver%pc%mat%nnz
    row_component = mod(poisson%solver%pc%mat%irn(entry)-1, poisson%solver%pc%mode_set_n)
    column_component = mod(poisson%solver%pc%mat%jcn(entry)-1, poisson%solver%pc%mode_set_n)
    if (row_component.ne.column_component) &
      cross_harmonic = max(cross_harmonic, abs(poisson%solver%pc%mat%val(entry)))
  enddo
  call MPI_AllReduce(cross_harmonic, cross_harmonic_global, 1, MPI_DOUBLE_PRECISION, &
                     MPI_MAX, MPI_COMM_WORLD, ierr)
  if (cross_harmonic_global.gt.1.d-13) error stop 'Poisson matrix couples harmonics.'

  allocate(rhs%val(poisson%rhs_global%n), phi_first(poisson%rhs_global%n), &
           phi_second(poisson%rhs_global%n), rhs_before_store(poisson%rhs_global%n))
  allocate(phi_old(n_tor,n_degrees,node_list%n_nodes))
  allocate(psi_before(n_tor,n_degrees,node_list%n_nodes))
  allocate(rho_before(n_tor,n_degrees,node_list%n_nodes))
  allocate(temperature_before(n_tor,n_degrees,node_list%n_nodes))
  do inode = 1, node_list%n_nodes
    phi_old(:,:,inode) = node_list%node(inode)%values(:,:,var_u)
    psi_before(:,:,inode) = node_list%node(inode)%values(:,:,var_psi)
    rho_before(:,:,inode) = node_list%node(inode)%values(:,:,var_rho)
    temperature_before(:,:,inode) = node_list%node(inode)%values(:,:,var_T)
  enddo
  rhs%n = poisson%rhs_global%n
  do index = 1, rhs%n
    rhs%val(index) = sin(0.17d0*real(index,8))
  enddo
  call poisson%set_rhs(rhs)
  call poisson%solve()
  if (.not.poisson%factorized) error stop 'First Poisson solve did not retain factors.'
  call poisson%gather()
  phi_first = poisson%phi_global%val
  rhs_before_store = poisson%rhs_global%val
  call poisson%store_phi()

  store_error = 0.d0
  addition_separation = huge(1.d0)
  delta_error = 0.d0
  other_field_error = 0.d0
  do inode = 1, node_list%n_nodes
    do idof = 1, n_degrees
      index = node_list%node(inode)%index(idof)
      do entry = 1, n_tor
        store_error = max(store_error, abs(node_list%node(inode)%values(entry,idof,var_u) - &
             phi_first(n_tor*(index-1)+entry)))
        addition_separation = min(addition_separation, &
             abs(node_list%node(inode)%values(entry,idof,var_u) - &
             (phi_old(entry,idof,inode)+phi_first(n_tor*(index-1)+entry))))
      enddo
    enddo
    delta_error = max(delta_error, maxval(abs(node_list%node(inode)%deltas(:,:,var_u))))
    other_field_error = max(other_field_error, &
         maxval(abs(node_list%node(inode)%values(:,:,var_psi)-psi_before(:,:,inode))))
    other_field_error = max(other_field_error, &
         maxval(abs(node_list%node(inode)%values(:,:,var_rho)-rho_before(:,:,inode))))
    other_field_error = max(other_field_error, &
         maxval(abs(node_list%node(inode)%values(:,:,var_T)-temperature_before(:,:,inode))))
  enddo
  if (store_error.gt.1.d-12 .or. addition_separation.le.1.d0 .or. &
      delta_error.gt.1.d-14 .or. other_field_error.gt.1.d-14) &
    error stop 'Absolute Poisson storage validation failed.'
  if (.not.poisson%factorized .or. .not.poisson%solution_gathered) &
    error stop 'Poisson storage changed solver state.'
  if (maxval(abs(poisson%rhs_global%val-rhs_before_store)).gt.0.d0) &
    error stop 'Poisson storage changed the global RHS.'

  rhs%val = 2.d0*rhs%val
  call poisson%set_rhs(rhs)
  call poisson%solve()
  if (.not.poisson%factorized) error stop 'Second Poisson solve lost retained factors.'
  call poisson%gather()
  phi_second = poisson%phi_global%val

  changed = maxval(abs(phi_second-phi_first))
  scale_error = maxval(abs(phi_second-2.d0*phi_first))

  rhs%val = 0.5d0*rhs%val
  call poisson%set_rhs(rhs)
  call poisson%solve()
  call poisson%gather()
  repeated = maxval(abs(poisson%phi_global%val-phi_first))

  if (changed.le.1.d-16 .or. scale_error.gt.1.d-10 .or. repeated.gt.1.d-10) then
    if (rank.eq.0) write(*,'(A,3ES14.5)') 'Poisson cycle FAILED: ', changed, scale_error, repeated
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
  endif
  if (rank.eq.0) then
    write(*,'(A,ES14.5)') 'changed RHS difference: ', changed
    write(*,'(A,ES14.5)') 'linear scaling error:   ', scale_error
    write(*,'(A,ES14.5)') 'repeated RHS error:     ', repeated
    write(*,'(A,ES14.5)') 'absolute storage error: ', store_error
    write(*,'(A,ES14.5)') 'untouched fields error: ', other_field_error
    write(*,*) 'Poisson solver cycle PASSED'
  endif

  call poisson%finalize()
  deallocate(rhs%val, phi_first, phi_second, rhs_before_store)
  deallocate(phi_old, psi_before, rho_before, temperature_before)
  call MPI_Finalize(ierr)
end program test_poisson_solver_cycle
