!> Fixed-potential boundary condition retained from the old model011.
module mod_boundary_conditions
  implicit none
contains
  subroutine boundary_conditions(my_id, node_list, element_list, bnd_node_list, local_elms, n_local_elms, &
                                 index_min, index_max, rhs_loc, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                 psi_bnd, R_xpoint, Z_xpoint, psi_xpoint, a_mat)
    use mod_parameters, only: n_order, n_vertex_max
    use data_structure, only: type_node_list, type_element_list, type_bnd_node_list, type_SP_MATRIX
    use phys_module, only: bcs, keep_n0_const
    use mod_assembly, only: boundary_conditions_add_one_entry
    use mod_node_indices, only: calculate_node_indices
    implicit none

    integer, intent(in)                    :: my_id, local_elms(*), n_local_elms, index_min, index_max
    type(type_node_list), intent(in)       :: node_list
    type(type_element_list), intent(in)    :: element_list
    type(type_bnd_node_list), intent(in)   :: bnd_node_list
    real*8, intent(inout)                  :: rhs_loc(*)
    logical, intent(in)                    :: xpoint2
    integer, intent(in)                    :: xcase2
    real*8, intent(in)                     :: R_axis, Z_axis, psi_axis, psi_bnd
    real*8, intent(in)                     :: R_xpoint(2), Z_xpoint(2), psi_xpoint(2)
    type(type_SP_MATRIX), intent(inout)    :: a_mat

    integer :: ie, edge, iv1, iv2, inode, idir, kk, ll, component, index_tmp, index_node
    integer :: node_indices((n_order+1)/2,(n_order+1)/2)
    real*8  :: penalty

    call calculate_node_indices(node_indices)

    ! The old model imposed homogeneous Dirichlet phi through a large diagonal
    ! penalty on every boundary type whose u flag is Dirichlet.  It constrained
    ! the nodal value and tangential Hermite derivatives for every real Fourier
    ! component, identically for n=0, cosine, and sine.
    do ie = 1, n_local_elms
      do edge = 1, n_vertex_max
        iv1 = edge
        iv2 = mod(edge,n_vertex_max)+1
        if (node_list%node(element_list%element(local_elms(ie))%vertex(iv1))%boundary.eq.0) cycle
        if (node_list%node(element_list%element(local_elms(ie))%vertex(iv2))%boundary.eq.0) cycle

        do idir = 1, 2
          if (idir.eq.1) then
            inode = element_list%element(local_elms(ie))%vertex(iv1)
          else
            inode = element_list%element(local_elms(ie))%vertex(iv2)
          endif
          if (.not.bcs(node_list%node(inode)%boundary)%dirichlet%u) cycle

          do component = a_mat%i_tor_min, a_mat%i_tor_max
            penalty = 1.d12
            if (keep_n0_const .and. component.eq.1) penalty = 1.d15
            do kk = 1, (n_order+1)/2
              do ll = 1, (n_order+1)/2
                if (mod(edge,2).eq.1 .and. ll.gt.1) cycle
                if (mod(edge,2).eq.0 .and. kk.gt.1) cycle
                index_tmp = node_indices(kk,ll)
                index_node = node_list%node(inode)%index(index_tmp)
                call boundary_conditions_add_one_entry(index_node, 1, component, index_node, 1, component, &
                                                       penalty, index_min, index_max, a_mat)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine boundary_conditions
end module mod_boundary_conditions
