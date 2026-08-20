!> Developer regression for the persistent model011 harmonic Poisson solve.
program test_poisson_solver_cycle
  use mpi_mod
  use data_structure,      only: type_RHS, init_node_list, nbthreads
  use mod_parameters,      only: n_tor, n_degrees, n_eq_var, var_psi, var_u, var_rho, var_T
  use phys_module,         only: use_mumps, use_pastix, use_strumpack, F0, central_mass, central_density, &
                                 filter_par, mode, mode_type, autodistribute_modes, n_mode_families, &
                                 modes_per_family, mode_families_modes
  use basis_at_gaussian,   only: initialise_basis
  use nodes_elements,      only: node_list, element_list, bnd_node_list, bnd_elm_list, aux_node_list
  use equil_info,          only: ES
  use mod_simulation_data, only: type_MHD_SIM
  use mod_poisson_solver,  only: poisson_solve_action
  use mod_poisson_rhs,     only: assemble_poisson_rhs, assemble_direct_poisson_rhs
  use mod_uncoupled_projection, only: assemble_projection_rhs
  use tr_module,           only: tr_meminit

  implicit none

  type(type_MHD_SIM), target :: mhd_sim
  type(poisson_solve_action) :: poisson
  type(type_RHS)             :: rhs, deposition_load, direct_rhs
  real*8, allocatable        :: element_deposition(:,:,:,:)
  real*8, allocatable        :: phi_first(:), phi_second(:)
  real*8, allocatable        :: phi_old(:,:,:), psi_before(:,:,:), rho_before(:,:,:)
  real*8, allocatable        :: temperature_before(:,:,:), rhs_before_store(:)
  real*8, allocatable        :: dense_matrix(:,:), dense_work(:,:), dense_rhs(:), dense_solution(:)
  real*8, allocatable        :: coupled_rhs_global(:)
  integer, allocatable       :: pivots(:)
  integer                    :: rank, n_tasks, ierr, inode, idof, index, entry
  integer                    :: row_component, column_component, cos_slot, sin_slot, spatial_dof, info
  integer                    :: family, member, mapping_family
  real*8                     :: changed, repeated, scale_error, cross_harmonic, cross_harmonic_global
  real*8                     :: store_error, addition_separation, delta_error, other_field_error
  real*8                     :: direct_error, direct_relative_error, density_norm
  real*8                     :: kcs, ksc, kcs_global, ksc_global, coupled_error, coupled_error_global
  real*8                     :: induced, induced_global, reference_scale, reference_scale_global, rebuild_error

  external dgesv

  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, ierr)
  if (n_tasks.ne.2) then
    if (rank.eq.0) write(*,*) 'This test requires exactly two MPI ranks.'
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  call tr_meminit(rank, n_tasks)
  call preset_parameters()
  ! The normal model011 test build has n_tor=3.  A temporary n_tor=7 test
  ! build exercises two families that each contain multiple physical n.
  if (n_tor.ge.7) then
    autodistribute_modes=.false.
    n_mode_families=2
    modes_per_family=0
    mode_families_modes=0
    modes_per_family(1:2)=(/3,4/)
    mode_families_modes(1,1:3)=(/1,2,3/)
    mode_families_modes(2,1:4)=(/4,5,6,7/)
  endif
  call det_modes()
  call initialise_basis()
  nbthreads = 1

  ! Exercise MUMPS explicitly; the production selector remains namelist-driven.
  use_mumps = .true.
  use_pastix = .false.
  use_strumpack = .false.
  F0 = 10.d0
  central_mass = 2.d0
  central_density = 2.5d0
  filter_par = 0.35d0

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
    node_list%node(inode)%values(1,1,var_psi) = 0.08d0*real(inode-1,8)
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
  if (rank.eq.0) then
    write(*,'(A)') 'Poisson physical-harmonic family mapping:'
    write(*,'(A,I0,A,I0)') 'physical n=0 global cos=1 family=',1
    do index=2,n_tor,2
      mapping_family=0
      do family=1,poisson%solver%pc%n_mode_families
        do member=1,poisson%solver%pc%modes_per_family(family)
          if (poisson%solver%pc%mode_families_modes(family,member).eq.index) mapping_family=family
        enddo
      enddo
      write(*,'(A,I0,A,I0,A,I0,A,I0)') 'physical n=',mode(index),' global cos=',index, &
        ' global sin=',index+1,' family=',mapping_family
    enddo
  endif
  call poisson%construct_matrix()
  if (.not.poisson%matrix_ready .or. poisson%factorized) &
    error stop 'Invalid state after Poisson matrix construction.'

  cross_harmonic = 0.d0
  kcs = 0.d0
  ksc = 0.d0
  cos_slot = 0
  sin_slot = 0
  do index = 1, poisson%solver%pc%mode_set_n
    if (mode_type(poisson%solver%pc%mode_set(index)).eq.'cos' .and. &
        mode(poisson%solver%pc%mode_set(index)).gt.0) cos_slot = index
    if (mode_type(poisson%solver%pc%mode_set(index)).eq.'sin') sin_slot = index
  enddo
  do entry = 1, poisson%solver%pc%mat%nnz
    row_component = mod(poisson%solver%pc%mat%irn(entry)-1, poisson%solver%pc%mode_set_n)
    column_component = mod(poisson%solver%pc%mat%jcn(entry)-1, poisson%solver%pc%mode_set_n)
    if (mode(poisson%solver%pc%mode_set(row_component+1)).ne. &
        mode(poisson%solver%pc%mode_set(column_component+1))) &
      cross_harmonic = max(cross_harmonic, abs(poisson%solver%pc%mat%val(entry)))
    if (row_component+1.eq.cos_slot .and. column_component+1.eq.sin_slot) &
      kcs = max(kcs,abs(poisson%solver%pc%mat%val(entry)))
    if (row_component+1.eq.sin_slot .and. column_component+1.eq.cos_slot) &
      ksc = max(ksc,abs(poisson%solver%pc%mat%val(entry)))
  enddo
  call MPI_AllReduce(cross_harmonic, cross_harmonic_global, 1, MPI_DOUBLE_PRECISION, &
                     MPI_MAX, MPI_COMM_WORLD, ierr)
  if (cross_harmonic_global.gt.1.d-13) error stop 'Poisson matrix couples distinct mode numbers.'
  call MPI_AllReduce(kcs,kcs_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(ksc,ksc_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  if (kcs_global.le.1.d-14 .or. ksc_global.le.1.d-14) &
    error stop 'Poisson sparse matrix is missing cosine/sine cross blocks.'

  ! Build an independent dense representation of each family matrix.  On two
  ! ranks each family communicator has one member, so this also checks the
  ! exact sparse local-to-family ordering used by the direct backend.
  allocate(dense_matrix(poisson%solver%pc%rhs%n,poisson%solver%pc%rhs%n))
  allocate(dense_work(poisson%solver%pc%rhs%n,poisson%solver%pc%rhs%n))
  allocate(dense_rhs(poisson%solver%pc%rhs%n),dense_solution(poisson%solver%pc%rhs%n))
  allocate(pivots(poisson%solver%pc%rhs%n),coupled_rhs_global(poisson%rhs_global%n))
  dense_matrix = 0.d0
  do entry=1,poisson%solver%pc%mat%nnz
    dense_matrix(poisson%solver%pc%mat%irn(entry),poisson%solver%pc%mat%jcn(entry)) = &
      dense_matrix(poisson%solver%pc%mat%irn(entry),poisson%solver%pc%mat%jcn(entry)) + &
      poisson%solver%pc%mat%val(entry)
  enddo

  if (cos_slot.gt.0 .and. sin_slot.gt.0) then
    write(*,'(A,I0,A,I0,A,I0,A,I0)') 'physical n=',mode(poisson%solver%pc%mode_set(cos_slot)), &
      ' global cos=',poisson%solver%pc%mode_set(cos_slot),' global sin=', &
      poisson%solver%pc%mode_set(sin_slot),' family=',poisson%solver%pc%family_id
  endif

  call check_coupled_rhs(cos_slot,sin_slot,'cosine-only')
  if (.not.poisson%factorized) error stop 'Coupled first solve did not retain factors.'
  call check_coupled_rhs(sin_slot,cos_slot,'sine-only')

  allocate(phi_first(poisson%rhs_global%n), &
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
  call assemble_poisson_rhs(node_list,element_list,poisson%solver%pc%local_elms, &
       poisson%solver%pc%n_local_elms,node_list,var_T,rhs)

  ! For this one-element mesh, construct the exact element deposition d=M*c
  ! represented by the reference projected-field Poisson load.  Only rank zero
  ! owns the synthetic particles; the reusable helper performs the MPI sum.
  density_norm=central_density*1.d20
  allocate(element_deposition(n_degrees,4,1,n_tor))
  element_deposition=0.d0
  if (rank.eq.0) then
    do inode=1,4
      do idof=1,n_degrees
        index=node_list%node(inode)%index(idof)
        element_deposition(idof,inode,1,:)=density_norm* &
             rhs%val(n_tor*(index-1)+1:n_tor*index)
      enddo
    enddo
  endif
  call assemble_projection_rhs(node_list,element_list,element_deposition,deposition_load,MPI_COMM_WORLD)
  call assemble_direct_poisson_rhs(deposition_load,direct_rhs)
  direct_error=maxval(abs(direct_rhs%val-rhs%val))
  direct_relative_error=maxval(abs(direct_rhs%val-rhs%val)/ &
       max(max(abs(rhs%val),abs(direct_rhs%val)),tiny(1.d0)))
  if (direct_error.gt.2.d-13 .or. direct_relative_error.gt.2.d-13) &
    error stop 'Direct and projection-based Poisson RHS differ.'
  call poisson%set_rhs(rhs)
  call poisson%solve()
  if (.not.poisson%factorized) error stop 'Poisson solve did not retain factors.'
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

  do inode = 1, node_list%n_nodes
    node_list%node(inode)%values(:,:,var_T)=2.d0*temperature_before(:,:,inode)
  enddo
  call assemble_poisson_rhs(node_list,element_list,poisson%solver%pc%local_elms, &
       poisson%solver%pc%n_local_elms,node_list,var_T,rhs)
  call poisson%set_rhs(rhs)
  call poisson%solve()
  if (.not.poisson%factorized) error stop 'Second Poisson solve lost retained factors.'
  call poisson%gather()
  phi_second = poisson%phi_global%val

  changed = maxval(abs(phi_second-phi_first))
  scale_error = maxval(abs(phi_second-2.d0*phi_first))

  do inode = 1, node_list%n_nodes
    node_list%node(inode)%values(:,:,var_T)=temperature_before(:,:,inode)
  enddo
  call assemble_poisson_rhs(node_list,element_list,poisson%solver%pc%local_elms, &
       poisson%solver%pc%n_local_elms,node_list,var_T,rhs)
  call poisson%set_rhs(rhs)
  call poisson%solve()
  call poisson%gather()
  repeated = maxval(abs(poisson%phi_global%val-phi_first))

  ! A matrix rebuild must invalidate and then replace the retained numerical
  ! factors; the same RHS must still produce the same coupled solution.
  call poisson%construct_matrix()
  if (poisson%factorized) error stop 'Poisson matrix rebuild retained stale factors.'
  call poisson%set_rhs(rhs)
  call poisson%solve()
  call poisson%gather()
  rebuild_error=maxval(abs(poisson%phi_global%val-phi_first))

  if (changed.le.1.d-16 .or. scale_error.gt.1.d-10 .or. repeated.gt.1.d-10 .or. rebuild_error.gt.1.d-10) then
    if (rank.eq.0) write(*,'(A,4ES14.5)') 'Poisson cycle FAILED: ', changed, scale_error, repeated, rebuild_error
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
  endif
  if (rank.eq.0) then
    write(*,'(A,ES14.5)') 'changed RHS difference: ', changed
    write(*,'(A,ES14.5)') 'linear scaling error:   ', scale_error
    write(*,'(A,ES14.5)') 'repeated RHS error:     ', repeated
    write(*,'(A,ES14.5)') 'matrix rebuild error:   ', rebuild_error
    write(*,'(A,ES14.5)') 'maximum Kcs entry:      ', kcs_global
    write(*,'(A,ES14.5)') 'maximum Ksc entry:      ', ksc_global
    write(*,'(A,ES14.5)') 'absolute storage error: ', store_error
    write(*,'(A,ES14.5)') 'untouched fields error: ', other_field_error
    write(*,'(A,ES14.5)') 'direct RHS abs error:    ', direct_error
    write(*,'(A,ES14.5)') 'direct RHS rel error:    ', direct_relative_error
    write(*,*) 'Poisson solver cycle PASSED'
  endif

  call poisson%finalize()
  deallocate(rhs%val, phi_first, phi_second, rhs_before_store)
  deallocate(deposition_load%val,direct_rhs%val,element_deposition)
  deallocate(phi_old, psi_before, rho_before, temperature_before)
  deallocate(dense_matrix,dense_work,dense_rhs,dense_solution,pivots,coupled_rhs_global)
  call MPI_Finalize(ierr)

contains

  subroutine check_coupled_rhs(source_slot,response_slot,label)
    integer, intent(in)          :: source_slot,response_slot
    character(len=*), intent(in) :: label

    dense_rhs=0.d0
    coupled_rhs_global=0.d0
    if (source_slot.gt.0 .and. response_slot.gt.0) then
      do spatial_dof=1,poisson%solver%pc%rhs%n/poisson%solver%pc%mode_set_n
        index=source_slot+(spatial_dof-1)*poisson%solver%pc%mode_set_n
        dense_rhs(index)=sin(0.37d0*real(spatial_dof,8))+0.25d0*cos(0.19d0*real(spatial_dof,8))
        coupled_rhs_global(poisson%solver%pc%row_index(index))=dense_rhs(index)
      enddo
    endif
    call MPI_AllReduce(MPI_IN_PLACE,coupled_rhs_global,size(coupled_rhs_global), &
                       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    dense_work=dense_matrix
    dense_solution=dense_rhs
    call dgesv(size(dense_solution),1,dense_work,size(dense_solution),pivots, &
               dense_solution,size(dense_solution),info)
    if (info.ne.0) error stop 'Independent dense Poisson reference solve failed.'

    rhs%n=poisson%rhs_global%n
    if (.not.associated(rhs%val)) allocate(rhs%val(rhs%n))
    rhs%val=coupled_rhs_global
    call poisson%set_rhs(rhs)
    call poisson%solve()
    call poisson%gather()

    coupled_error=0.d0
    induced=0.d0
    reference_scale=1.d0
    if (source_slot.gt.0 .and. response_slot.gt.0) then
      coupled_error=maxval(abs(poisson%phi_global%val(poisson%solver%pc%row_index)-dense_solution))
      reference_scale=max(1.d0,maxval(abs(dense_solution)))
      do spatial_dof=1,poisson%solver%pc%rhs%n/poisson%solver%pc%mode_set_n
        index=response_slot+(spatial_dof-1)*poisson%solver%pc%mode_set_n
        induced=max(induced,abs(dense_solution(index)))
      enddo
    endif
    call MPI_AllReduce(coupled_error,coupled_error_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(induced,induced_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(reference_scale,reference_scale_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (coupled_error_global.gt.2.d-10*reference_scale_global .or. induced_global.le.1.d-14) &
      error stop 'Coupled Poisson solve disagrees with dense reference or did not induce the paired component.'
    if (rank.eq.0) write(*,'(A,A,A,2ES14.5)') trim(label),' coupled solve: ', &
      'error, induced=',coupled_error_global,induced_global
  end subroutine check_coupled_rhs
end program test_poisson_solver_cycle
