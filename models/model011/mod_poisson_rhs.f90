!> Physical model011 Poisson load assembly from a projected charge field.
module mod_poisson_rhs
  implicit none
  private
  public :: assemble_poisson_rhs, assemble_direct_poisson_rhs
contains

  !> Convert the already assembled particle deposition load directly to the
  !! physical model011 Poisson RHS.  The deposition contains charge sign,
  !! macroparticle weight and real-Fourier test-function sampling; only the
  !! legacy density normalization remains.
  subroutine assemble_direct_poisson_rhs(deposition_load, rhs_global)
    use phys_module, only: central_density
    use data_structure, only: type_RHS
    implicit none

    type(type_RHS), intent(in)    :: deposition_load
    type(type_RHS), intent(inout) :: rhs_global

    if (.not.associated(deposition_load%val)) error stop 'Unallocated direct Poisson deposition load.'
    if (central_density.le.0.d0) error stop 'Invalid model011 central density.'
    rhs_global%n = deposition_load%n
    if (associated(rhs_global%val)) then
      if (size(rhs_global%val).ne.rhs_global%n) deallocate(rhs_global%val)
    endif
    if (.not.associated(rhs_global%val)) allocate(rhs_global%val(rhs_global%n))
    rhs_global%val = deposition_load%val/(central_density*1.d20)
  end subroutine assemble_direct_poisson_rhs

  !> Assemble projected charge coefficients into normal global JOREK ordering.
  !! `charge_nodes` is the output of the particle projection, not its element
  !! deposition RHS. Charge signs and macroparticle weights are already in it.
  subroutine assemble_poisson_rhs(node_list, element_list, local_elms, n_local_elms, &
                                  charge_nodes, charge_variable, rhs_global)
    use mpi_mod
    use mod_parameters, only: n_vertex_max, n_degrees, n_tor, n_plane, n_var
    use phys_module, only: central_density, mode, mode_type
    use data_structure, only: type_node_list, type_element_list, type_RHS
    use gauss, only: n_gauss, wgauss
    use basis_at_gaussian, only: H, H_s, H_t
    use mod_poisson_element_kernel, only: accumulate_poisson_load_mass, apply_poisson_load_harmonics
    implicit none

    type(type_node_list), intent(in)    :: node_list, charge_nodes
    type(type_element_list), intent(in) :: element_list
    integer, intent(in)                 :: local_elms(:), n_local_elms, charge_variable
    type(type_RHS), intent(inout)       :: rhs_global
    integer, parameter                 :: basis_size=n_vertex_max*n_degrees
    real*8, allocatable                :: rhs_local(:)
    real*8                             :: value(basis_size), charge(basis_size,n_tor)
    real*8                             :: load_mass(basis_size,basis_size), element_rhs(basis_size,n_tor)
    real*8                             :: x_s, x_t, y_s, y_t, big_r, xjac, weight
    integer                            :: ife, ielm, i, j, a, ms, mt, inode, idof, itor
    integer                            :: index_node, index_global, ierr

    if (n_var.ne.1) error stop 'model011 Poisson RHS assembly requires n_var=1.'
    if (n_local_elms.gt.size(local_elms)) error stop 'Invalid local Poisson element list.'
    if (charge_variable.lt.1 .or. charge_variable.gt.charge_nodes%n_values) &
      error stop 'Invalid projected charge variable.'

    rhs_global%n = node_list%n_dof
    if (associated(rhs_global%val)) then
      if (size(rhs_global%val).ne.rhs_global%n) deallocate(rhs_global%val)
    endif
    if (.not.associated(rhs_global%val)) allocate(rhs_global%val(rhs_global%n))
    allocate(rhs_local(rhs_global%n))
    rhs_global%val = 0.d0
    rhs_local = 0.d0

    do ife = 1, n_local_elms
      ielm = local_elms(ife)
      load_mass = 0.d0
      charge = 0.d0

      do i = 1, n_vertex_max
        inode = element_list%element(ielm)%vertex(i)
        do j = 1, n_degrees
          a = (i-1)*n_degrees+j
          charge(a,:) = charge_nodes%node(inode)%values(:,j,charge_variable)
        enddo
      enddo

      do ms = 1, n_gauss
        do mt = 1, n_gauss
          x_s=0.d0; x_t=0.d0; y_s=0.d0; y_t=0.d0; big_r=0.d0
          do i = 1, n_vertex_max
            inode = element_list%element(ielm)%vertex(i)
            do j = 1, n_degrees
              a = (i-1)*n_degrees+j
              value(a) = element_list%element(ielm)%size(i,j)*H(i,j,ms,mt)
              big_r = big_r + node_list%node(inode)%x(1,j,1)*value(a)
              x_s = x_s + node_list%node(inode)%x(1,j,1)* &
                    element_list%element(ielm)%size(i,j)*H_s(i,j,ms,mt)
              x_t = x_t + node_list%node(inode)%x(1,j,1)* &
                    element_list%element(ielm)%size(i,j)*H_t(i,j,ms,mt)
              y_s = y_s + node_list%node(inode)%x(1,j,2)* &
                    element_list%element(ielm)%size(i,j)*H_s(i,j,ms,mt)
              y_t = y_t + node_list%node(inode)%x(1,j,2)* &
                    element_list%element(ielm)%size(i,j)*H_t(i,j,ms,mt)
            enddo
          enddo
          xjac = x_s*y_t-x_t*y_s
          weight = wgauss(ms)*wgauss(mt)
          call accumulate_poisson_load_mass(weight,big_r,xjac,value,load_mass)
        enddo
      enddo

      call apply_poisson_load_harmonics(load_mass,charge,mode,mode_type,n_plane, &
                                        central_density*1.d20,element_rhs)
      do i = 1, n_vertex_max
        inode = element_list%element(ielm)%vertex(i)
        do idof = 1, n_degrees
          index_node = node_list%node(inode)%index(idof)
          if (index_node.le.0) cycle
          a = (i-1)*n_degrees+idof
          do itor = 1, n_tor
            index_global = n_tor*n_var*(index_node-1)+itor
            rhs_local(index_global) = rhs_local(index_global)+element_rhs(a,itor)
          enddo
        enddo
      enddo
    enddo

    call MPI_AllReduce(rhs_local,rhs_global%val,rhs_global%n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate(rhs_local)
  end subroutine assemble_poisson_rhs
end module mod_poisson_rhs
