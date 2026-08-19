!> No natural/open-boundary element term is present in the old model011 operator.
module mod_boundary_matrix_open
  implicit none
contains
  subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, &
                                  psi_axis, psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max)
    use mod_parameters, only: n_vertex_max
    use data_structure, only: type_element, type_node
    implicit none
    integer, intent(in)                   :: vertex(2), direction(2), xcase2, i_tor_min, i_tor_max
    type(type_element), intent(in)        :: element
    type(type_node), intent(in)           :: nodes(n_vertex_max)
    logical, intent(in)                   :: xpoint2
    real*8, intent(in)                    :: R_axis, Z_axis, psi_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
    real*8, allocatable, intent(inout)    :: ELM(:,:), RHS(:)
  end subroutine boundary_matrix_open
end module mod_boundary_matrix_open
