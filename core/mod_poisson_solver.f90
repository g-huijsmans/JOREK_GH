!> Stateful infrastructure for a direct harmonic Poisson solve.
module mod_poisson_solver
  use mod_event,           only: action, event
  use mod_particle_sim,    only: particle_sim
  use mod_simulation_data, only: type_MHD_SIM
  use mod_sparse_data,     only: type_SP_SOLVER
  use data_structure,      only: type_RHS

  implicit none

  private
  public :: poisson_solve_action

  type, extends(action) :: poisson_solve_action
    type(type_SP_SOLVER)        :: solver
    type(type_RHS)              :: rhs_global
    type(type_RHS)              :: phi_global
    type(type_MHD_SIM), pointer :: mhd_sim => null()

    logical :: setup_done  = .false.
    logical :: matrix_ready = .false.
    logical :: factorized   = .false.
    logical :: solution_gathered = .false.
  contains
    procedure :: setup    => setup_poisson_solver
    procedure :: construct_matrix => construct_poisson_matrix
    procedure :: set_rhs  => set_poisson_rhs
    procedure :: solve    => solve_poisson_system
    procedure :: gather   => gather_poisson_solution
    procedure :: store_phi => store_poisson_solution
    procedure :: residual_norm => poisson_residual_norm
    procedure :: finalize => finalize_poisson_solver
    procedure :: do       => do_poisson_solve
  end type poisson_solve_action

contains

  !> Set up mode-family communicators, harmonic structure, and global vectors.
  subroutine setup_poisson_solver(this, mhd_sim, comm)
    use mod_preconditioner, only: initialize_preconditioner
#ifdef DIRECT_CONSTRUCTION
    use mod_direct_construction, only: setup_pc_structure
#endif

    class(poisson_solve_action), intent(inout) :: this
    type(type_MHD_SIM), intent(inout), target  :: mhd_sim
    integer, intent(in)                        :: comm
    integer                                    :: n_global

    if (this%setup_done) return

#ifndef DIRECT_CONSTRUCTION
    write(*,*) 'ERROR: The dedicated Poisson solver requires DIRECT_CONSTRUCTION.'
    error stop
#else
    this%mhd_sim => mhd_sim

    call this%solver%setup()
    call initialize_preconditioner(this%solver%pc, comm)
    call validate_poisson_mode_families(this%solver%pc)
    call setup_pc_structure(this%solver%pc, mhd_sim)

    n_global = mhd_sim%node_list%n_dof
    allocate(this%rhs_global%val(n_global))
    allocate(this%phi_global%val(n_global))
    this%rhs_global%n = n_global
    this%phi_global%n = n_global
    this%rhs_global%val = 0.d0
    this%phi_global%val = 0.d0

    this%matrix_ready = .false.
    this%factorized = .false.
    this%solution_gathered = .false.
    this%setup_done = .true.
#endif
  end subroutine setup_poisson_solver


  !> Check the model011 real-Fourier requirements without changing the generic
  !! mode-family distribution.  The direct-construction implementation uses a
  !! contiguous component interval for each family; projection preconditioners
  !! remain free to use the more general overlapping/noncontiguous metadata.
  subroutine validate_poisson_mode_families(pc)
    use data_structure, only: type_PRECOND
    use mod_parameters, only: n_tor
    use phys_module,    only: mode, mode_type

    type(type_PRECOND), intent(in) :: pc
    integer                       :: component, family, member, previous
    integer                       :: cos_family, sin_family
    integer                       :: membership(n_tor)

    membership=0
    do family=1,pc%n_mode_families
      previous=0
      do member=1,pc%modes_per_family(family)
        component=pc%mode_families_modes(family,member)
        if (component.lt.1 .or. component.gt.n_tor) &
          error stop 'Poisson mode family contains an invalid global Fourier component.'
        membership(component)=membership(component)+1
        if (member.gt.1 .and. component.ne.previous+1) &
          error stop 'Direct Poisson construction requires contiguous components within a mode family.'
        previous=component
      enddo
    enddo

    if (any(membership.ne.1)) &
      error stop 'Poisson mode families must contain every global Fourier component exactly once.'
    if (mode(1).ne.0 .or. mode_type(1).ne.'cos') &
      error stop 'Unexpected JOREK n=0 real-Fourier component convention.'

    do component=2,n_tor,2
      if (component+1.gt.n_tor) &
        error stop 'Poisson real-Fourier storage ends with an incomplete cosine/sine pair.'
      if (mode_type(component).ne.'cos' .or. mode_type(component+1).ne.'sin' .or. &
          mode(component).ne.mode(component+1)) &
        error stop 'Unexpected JOREK nonzero real-Fourier component convention.'
      cos_family=0
      sin_family=0
      do family=1,pc%n_mode_families
        if (any(pc%mode_families_modes(family,1:pc%modes_per_family(family)).eq.component)) &
          cos_family=family
        if (any(pc%mode_families_modes(family,1:pc%modes_per_family(family)).eq.component+1)) &
          sin_family=family
      enddo
      if (cos_family.ne.sin_family) &
        error stop 'Poisson cosine/sine components for one physical harmonic are split across mode families.'
    enddo
  end subroutine validate_poisson_mode_families


  !> Assemble the model-selected harmonic matrix into the persistent PC storage.
  subroutine construct_poisson_matrix(this)
    use construct_matrix_mod, only: assemble_harmonic_matrix => construct_matrix
    use mod_sparse,           only: invalidate_pc_factorization

    class(poisson_solve_action), intent(inout) :: this

    if (.not.this%setup_done) then
      write(*,*) 'ERROR: Poisson matrix construction requested before setup.'
      error stop
    endif

    if (this%factorized) call invalidate_pc_factorization(this%solver)

    call assemble_harmonic_matrix(this%mhd_sim, this%solver%pc%local_elms, &
         this%solver%pc%n_local_elms, this%solver%pc%mat, this%solver%pc%rhs, &
         harmonic_matrix=.true.)

    this%matrix_ready = .true.
    this%factorized = .false.
    this%solution_gathered = .false.
  end subroutine construct_poisson_matrix


  !> Copy a global JOREK RHS and redistribute it to this rank's mode family.
  subroutine set_poisson_rhs(this, rhs)
    use mod_preconditioner, only: update_pc_rhs

    class(poisson_solve_action), intent(inout) :: this
    type(type_RHS), intent(in)                 :: rhs

    if (.not.this%setup_done) then
      write(*,*) 'ERROR: Poisson solver RHS set before setup.'
      error stop
    endif
    if (.not.associated(rhs%val)) then
      write(*,*) 'ERROR: Poisson solver received an unallocated global RHS.'
      error stop
    endif
    if (rhs%n.ne.this%rhs_global%n) then
      write(*,*) 'ERROR: Invalid Poisson RHS size:', rhs%n, ' expected:', this%rhs_global%n
      error stop
    endif

    this%rhs_global%val(1:rhs%n) = rhs%val(1:rhs%n)
    call update_pc_rhs(this%solver%pc, this%rhs_global)
    this%solution_gathered = .false.
  end subroutine set_poisson_rhs


  !> Solve the mode-family system, factorizing only on the first call.
  subroutine solve_poisson_system(this)
    use mod_sparse, only: solve_pc_direct

    class(poisson_solve_action), intent(inout) :: this
    logical                                   :: solve_only

    if (.not.this%setup_done) then
      write(*,*) 'ERROR: Poisson solve requested before setup.'
      error stop
    endif
    if (.not.this%matrix_ready) then
      write(*,*) 'ERROR: Poisson solve requested before matrix assembly.'
      error stop
    endif

    solve_only = this%factorized
    call solve_pc_direct(this%solver, solve_only, -1)
    this%factorized = .true.
    this%solution_gathered = .false.
  end subroutine solve_poisson_system


  !> Gather the in-place mode-family solution into global JOREK ordering.
  subroutine gather_poisson_solution(this)
    use mod_preconditioner, only: gather_solution

    class(poisson_solve_action), intent(inout) :: this

    if (.not.this%factorized) then
      write(*,*) 'ERROR: Poisson solution gather requested before a solve.'
      error stop
    endif
    call gather_solution(this%solver%pc, this%phi_global)
    this%solution_gathered = .true.
  end subroutine gather_poisson_solution


  !> Evaluate the residual of the first retained nonzero physical harmonic
  !! using the assembled mode-family matrix and the RHS supplied to the solve.
  !! Both real-Fourier components of the selected harmonic are included.
  subroutine poisson_residual_norm(this, harmonic, absolute_norm, relative_norm, available)
    use mpi_mod
    use mod_parameters, only: n_tor
    use phys_module, only: mode

    class(poisson_solve_action), intent(in) :: this
    integer, intent(out)                    :: harmonic
    real*8, intent(out)                     :: absolute_norm, relative_norm
    logical, intent(out)                    :: available

    real*8, allocatable :: ax(:), solver_x(:)
    real*8              :: residual_sq_local, residual_sq_global
    real*8              :: rhs_sq_local, rhs_sq_global, residual
    integer             :: component, slot, i, ierr, available_local, available_global
    logical             :: owns_rows, selected_family

    harmonic = 0
    do component=1,n_tor
      if (mode(component).ne.0) then
        harmonic = mode(component)
        exit
      endif
    enddo

    available_local = 0
    if (this%solution_gathered .and. harmonic.ne.0 .and. &
        associated(this%solver%pc%mat%irn) .and. &
        associated(this%solver%pc%mat%jcn) .and. &
        associated(this%solver%pc%mat%val)) available_local = 1
    call MPI_AllReduce(available_local,available_global,1,MPI_INTEGER,MPI_MIN, &
                       this%solver%pc%comm,ierr)
    available = available_global.eq.1
    if (.not.available) then
      absolute_norm = -1.d0
      relative_norm = -1.d0
      return
    endif

    allocate(ax(this%solver%pc%rhs%n),solver_x(this%solver%pc%rhs%n))
    ax = 0.d0
    solver_x = this%solver%pc%rhs%val
    if (this%solver%pc%mat%scaled) then
      if (.not.associated(this%solver%pc%mat%column_scaling)) then
        available_local = 0
      endif
    endif
    call MPI_AllReduce(available_local,available_global,1,MPI_INTEGER,MPI_MIN, &
                       this%solver%pc%comm,ierr)
    available = available_global.eq.1
    if (.not.available) then
      absolute_norm = -1.d0
      relative_norm = -1.d0
      deallocate(ax,solver_x)
      return
    endif
    if (this%solver%pc%mat%scaled) then
      solver_x = solver_x*this%solver%pc%mat%column_scaling
    endif

    selected_family = any(mode(this%solver%pc%mode_set).eq.harmonic)
    if (selected_family) then
      do i=1,this%solver%pc%mat%nnz
        ax(this%solver%pc%mat%irn(i)) = ax(this%solver%pc%mat%irn(i)) + &
             this%solver%pc%mat%val(i)*solver_x(this%solver%pc%mat%jcn(i))
      enddo
    endif

    ! A distributed matrix contributes its owned rows on every family rank;
    ! a replicated matrix contributes once, from the family master.
    owns_rows = this%solver%pc%mat%row_distributed .or. this%solver%pc%my_id_n.eq.0
    residual_sq_local = 0.d0
    rhs_sq_local = 0.d0
    if (selected_family .and. owns_rows) then
      do i=1,this%solver%pc%rhs%n
        slot = mod(i-1,this%solver%pc%mode_set_n)+1
        if (mode(this%solver%pc%mode_set(slot)).ne.harmonic) cycle
        residual = ax(i)-this%rhs_global%val(this%solver%pc%row_index(i))
        residual_sq_local = residual_sq_local+residual*residual
        rhs_sq_local = rhs_sq_local + &
             this%rhs_global%val(this%solver%pc%row_index(i))**2
      enddo
    endif

    call MPI_AllReduce(residual_sq_local,residual_sq_global,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                       this%solver%pc%comm,ierr)
    call MPI_AllReduce(rhs_sq_local,rhs_sq_global,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                       this%solver%pc%comm,ierr)
    absolute_norm = sqrt(residual_sq_global)
    relative_norm = absolute_norm/max(sqrt(rhs_sq_global),tiny(1.d0))
    deallocate(ax,solver_x)
  end subroutine poisson_residual_norm


  !> Store the gathered equation solution as absolute physical node values.
  !! This deliberately does not use update_values: phi is not a Newton
  !! increment, n=0 must not be frozen, and its delta storage must be cleared.
  subroutine store_poisson_solution(this)
    use mod_parameters,    only: n_tor, n_var, n_degrees, n_vertex_max, var_index
    use phys_module,       only: treat_axis
    use data_structure,    only: type_element_list, type_node_list
    use mod_basisfunctions, only: basisfunctions
    use mod_axis_treatment, only: new_to_old_dofs_on_the_axis

    class(poisson_solve_action), intent(inout) :: this

    type(type_element_list), pointer :: element_list
    type(type_node_list), pointer    :: node_list
    real*8                           :: H(4,n_degrees), H_s(4,n_degrees)
    real*8                           :: H_t(4,n_degrees), H_st(4,n_degrees)
    real*8                           :: new_dofs(4), old_dofs(4)
    real*8                           :: phi, phi_s, phi_t, phi_st
    real*8                           :: lambda, mu
    integer                          :: parent(2), vertices(n_vertex_max)
    integer                          :: inode, idof, ivar, itor, index_node, index
    integer                          :: ielm, vertex, parent_dof, ivar_store

    if (.not.this%solution_gathered) then
      write(*,*) 'ERROR: Poisson storage requested before gathering the solution.'
      error stop
    endif

    node_list => this%mhd_sim%node_list
    element_list => this%mhd_sim%element_list

    ! Absolute Poisson storage has no increment semantics.  Only equation
    ! fields selected by var_index are touched; imported background fields
    ! outside that mapping are preserved.
    do ivar = 1, n_var
      ivar_store = var_index(ivar)
      do inode = 1, node_list%n_nodes
        node_list%node(inode)%deltas(:,:,ivar_store) = 0.d0
      enddo
    enddo

    ! Store independent nodes first.  On-axis equations use transformed
    ! degrees of freedom, so convert them back to physical Hermite storage.
    do inode = 1, node_list%n_nodes
      if (node_list%node(inode)%constrained) cycle

      do ivar = 1, n_var
        ivar_store = var_index(ivar)
        do itor = 1, n_tor
          if (treat_axis .and. node_list%node(inode)%axis_node) then
            do idof = 1, n_degrees
              index_node = node_list%node(inode)%index(idof)
              index = n_tor*n_var*(index_node-1) + n_tor*(ivar-1) + itor
              new_dofs(idof) = this%phi_global%val(index)
            enddo
            call new_to_old_dofs_on_the_axis(node_list, inode, new_dofs, old_dofs)
            node_list%node(inode)%values(itor,:,ivar_store) = old_dofs
          else
            do idof = 1, n_degrees
              index_node = node_list%node(inode)%index(idof)
              index = n_tor*n_var*(index_node-1) + n_tor*(ivar-1) + itor
              if (index.gt.0) node_list%node(inode)%values(itor,idof,ivar_store) = &
                   this%phi_global%val(index)
            enddo
          endif
        enddo
      enddo
    enddo

    ! Constrained nodes have no independent equation entries.  Reconstruct
    ! their value and Hermite derivatives from the two parent vertices using
    ! exactly the interpolation and scale factors used by update_values.
    do inode = 1, node_list%n_nodes
      if (.not.node_list%node(inode)%constrained) cycle

      lambda = node_list%node(inode)%ref_lambda
      mu = node_list%node(inode)%ref_mu
      ielm = node_list%node(inode)%parent_elem
      parent = node_list%node(inode)%parents
      vertices = element_list%element(ielm)%vertex
      call basisfunctions(lambda, mu, H, H_s, H_t, H_st)

      do ivar = 1, n_var
        ivar_store = var_index(ivar)
        do itor = 1, n_tor
          phi = 0.d0
          phi_s = 0.d0
          phi_t = 0.d0
          phi_st = 0.d0
          do vertex = 1, n_vertex_max
            if (vertices(vertex).ne.parent(1) .and. vertices(vertex).ne.parent(2)) cycle
            do parent_dof = 1, n_degrees
              phi = phi + node_list%node(vertices(vertex))%values(itor,parent_dof,ivar_store) * &
                   H(vertex,parent_dof)*element_list%element(ielm)%size(vertex,parent_dof)
              phi_s = phi_s + node_list%node(vertices(vertex))%values(itor,parent_dof,ivar_store) * &
                   H_s(vertex,parent_dof)*element_list%element(ielm)%size(vertex,parent_dof)
              phi_t = phi_t + node_list%node(vertices(vertex))%values(itor,parent_dof,ivar_store) * &
                   H_t(vertex,parent_dof)*element_list%element(ielm)%size(vertex,parent_dof)
              phi_st = phi_st + node_list%node(vertices(vertex))%values(itor,parent_dof,ivar_store) * &
                   H_st(vertex,parent_dof)*element_list%element(ielm)%size(vertex,parent_dof)
            enddo
          enddo
          node_list%node(inode)%values(itor,1,ivar_store) = phi
          node_list%node(inode)%values(itor,2,ivar_store) = phi_s/3.d0
          node_list%node(inode)%values(itor,3,ivar_store) = phi_t/3.d0
          node_list%node(inode)%values(itor,4,ivar_store) = phi_st/9.d0
        enddo
      enddo
    enddo
  end subroutine store_poisson_solution


  !> Action wrapper: perform only the already-configured solve and gather.
  subroutine do_poisson_solve(this, sim, ev)
    class(poisson_solve_action), intent(inout) :: this
    type(particle_sim), intent(inout)          :: sim
    type(event), intent(inout), optional       :: ev

    call this%solve()
    call this%gather()
  end subroutine do_poisson_solve


  !> Release backend state, preconditioner storage/communicators, and vectors.
  subroutine finalize_poisson_solver(this)
    use mod_preconditioner, only: reset_preconditioner

    class(poisson_solve_action), intent(inout) :: this

    call this%solver%finalize()
    call reset_preconditioner(this%solver%pc)

    if (associated(this%rhs_global%val)) deallocate(this%rhs_global%val)
    if (associated(this%phi_global%val)) deallocate(this%phi_global%val)
    this%rhs_global%val => null()
    this%phi_global%val => null()
    this%rhs_global%n = 0
    this%phi_global%n = 0
    this%mhd_sim => null()

    this%setup_done = .false.
    this%matrix_ready = .false.
    this%factorized = .false.
    this%solution_gathered = .false.
  end subroutine finalize_poisson_solver

end module mod_poisson_solver
